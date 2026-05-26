# ============================================================================
# Phase 1: ENVIRONMENT SETUP & CONFIGURATION
# ============================================================================
rm(list = ls(all.names = TRUE))
gc()
global_seed <- 624
set.seed(global_seed)
options(dplyr.summarise.inform = FALSE)
options(tensorflow.extract.one_based = TRUE)
Sys.setenv(PYTHONHASHSEED = as.character(global_seed))

source("Functions/Helper.R") # Helper functions for metrics, predictions, and data loading
source('Functions/wLC-NB.R') # Lee-Carter model fitting
source("Functions/CoherentLCDrift.R") # Coherent Lee-Carter κ2 forecast with AR(1)+drift
source("Functions/Models-QuantileLSTM.R") # Model definitions (CNN, LSTM, GNN)

# Install and load required packages
packages <- c(
    'dplyr',
    'sf',
    'ncdf4',
    'lubridate',
    'tidyverse',
    'tibble',
    'stats',
    'tidyr',
    'ISOweek',
    'ecmwfr',
    'httr',
    'utils',
    'sp',
    'spdep',
    'eurostat',
    'keras3',
    'keras',
    'abind',
    'ggplot2',
    'forecast',
    'data.table',
    'ggpubr',
    'strucchange',
    'foreach',
    'doSNOW',
    'snow',
    'progress'
)
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
    install.packages(new.packages)
}
suppressMessages(sapply(packages, require, character.only = TRUE))

################################################################################
# Hyperparameters
################################################################################
seq_length_cnn_gnn <- 4
seq_length_mort <- 2
batch_size <- 16

# CNN-LSTM specific training controls
epochs_final_cnn <- 70
patience_cnn <- 15

# MortFCNet specific training controls
epochs_final_mort <- 70
patience_mort <- 15

# GNN-LSTM specific training controls
epochs_final_gnn <- 200
patience_gnn <- 15

# Timing storage
model_timings <- list()

# Best hyperparameter choices (priority: GNN-LSTM > CNN-LSTM > MortFCNet)
# Preferred GNN-LSTM configuration
best_gnn <- list(
    gnn_hidden_units = c(16),
    lstm_units = 16,
    hidden_size_1 = 32,
    hidden_size_2 = 16,
    hidden_size_3 = 8,
    age_embedding_dim = 4,
    region_embedding_dim = 6,
    dropout_rate = 0.10,
    learning_rate = 1e-04
)

# Preferred CNN-LSTM configuration (second priority)
best_cnn <- list(
    cnn_filters = c(16),
    lstm_units = 16,
    hidden_size_1 = 32,
    hidden_size_2 = 16,
    hidden_size_3 = 8,
    region_embedding_dim = 6,
    age_embedding_dim = 4,
    dropout_rate = 0.1,
    learning_rate = 5e-05
)

# Preferred MortFCNet configuration (third priority)
best_mort <- list(
    lstm_units = 32,
    hidden_size_1 = 8,
    hidden_size_2 = 4,
    hidden_size_3 = 2,
    dropout_rate = 0.2,
    learning_rate = 5e-05
)

# Embedding toggles for ablation experiments
use_age_embedding_cnn <- FALSE
use_region_embedding_cnn <- TRUE
use_age_embedding_gnn <- FALSE
use_region_embedding_gnn <- TRUE

################################################################################
# TRAIN/TEST SPLIT CONFIGURATION
################################################################################
train_start <- 1990
train_end <- 2014
test_start <- 2015
test_end <- 2019

train_years <- train_start:train_end
test_years <- test_start:test_end

# With validation removed, use all training years for preprocessing stats.
fit_years <- train_years

# ==============================================================================
# PHASE 2: DATA PREPARATION
# ==============================================================================
# Weekly Mortality Model

# Prepare data
Df <- fread(
    file = 'C:/1 My Code/Climate-Mortality-Forecasting/Results/INSEE/Df_nuts2_final.txt'
)

if ("NUTS2_ID" %in% names(Df)) {
    names(Df)[names(Df) == "NUTS2_ID"] <- "region"
}
Df <- Df %>% rename(deaths = Deaths, expo = expo)

# Restrict to metropolitan NUTS2 regions and aggregate deaths/exposure into
# age groups 65, 70, 75, 80, 85, 90+ for weekly region-age analysis.
overseas_regions <- c(paste0('FRY', 1:5), 'FRM0')

DfM_all <- Df %>%
    dplyr::filter(
        region != 94,
        !region %in% overseas_regions,
        age >= 65,
        ISOYear <= test_end,
        ISOYear >= train_start,
        ISOWeek <= 52
    ) %>%
    dplyr::mutate(
        age = dplyr::case_when(
            age >= 90 ~ 90L,
            age >= 85 ~ 85L,
            age >= 80 ~ 80L,
            age >= 75 ~ 75L,
            age >= 70 ~ 70L,
            TRUE ~ 65L
        )
    ) %>%
    dplyr::group_by(ISODate, ISOYear, ISOWeek, age, region) %>%
    dplyr::summarise(
        deaths = sum(deaths, na.rm = TRUE),
        expo = sum(expo, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::arrange(age, ISODate)

# Training data only
DfM_train <- DfM_all %>%
    dplyr::filter(ISOYear >= train_start, ISOYear <= train_end)

# Test data
DfM_test <- DfM_all %>%
    dplyr::filter(ISOYear >= test_start, ISOYear <= test_end)

# Death counts and exposures - matrix format (TRAINING DATA ONLY)
dtxr <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM_train)
etxr <- stats::xtabs(expo ~ ISODate + age + region, data = DfM_train)

# Age range, year range, week range
xv <- as.integer(colnames(dtxr))
train_row_years <- lubridate::isoyear(as.Date(rownames(dtxr)))
isoyv <- resolve_year_index(
    train_row_years,
    train_years,
    "training LC year index"
)
isowv <- isoweek(as.Date(rownames(dtxr)))
R <- dim(dtxr)[3]

# Mark June-July 1997 as missing (data imputation artefact)
dates_dtxr <- as.Date(dimnames(dtxr)[[1]])
ind.rm <- which(
    dates_dtxr >= as.Date('1997-06-01') &
        dates_dtxr <= as.Date('1997-07-31')
)
dtxr[ind.rm, , ] <- NA

# Avoid numerical issues after the missing-week mask is applied.
dtxr[dtxr == 0] <- 0.001

# ==============================================================================
# PHASE 3: BASELINE LEE-CARTER MODELING & KAPPA2 FORECAST
# ==============================================================================
lca <- fit701M.nb(xv, isoyv, isowv, etxr, dtxr, xv * 0, 'ALL')

# Fill imputed weeks with model-fitted values
dtxr[ind.rm, , ] <- lca$mhat[ind.rm, , ] * etxr[ind.rm, , ]

################################################################################
# FORECAST KAPPA2 FOR TEST PERIOD — Li & Lee Coherent Decomposition
# κ_{t,r} = κ_t^(F) + u_{t,r}
#   κ_t^(F) follows RW with drift (common national factor)
#   u_{t,r} follows AR(1) (mean-reverting regional deviations)
################################################################################
kappa2_fitted <- lca$kappa2
n_years_test <- test_end - test_start + 1
n_years_forecast <- n_years_test

colnames(kappa2_fitted) <- dimnames(dtxr)[[3]]
rownames(kappa2_fitted) <- as.character(train_start:train_end)

coherent_res <- forecast_kappa2_coherent(
    kappa2_matrix = kappa2_fitted,
    n_ahead = n_years_forecast,
    weights = NULL
)
kappa2_forecast <- coherent_res$kappa2_forecast
kappa2_full <- rbind(kappa2_fitted, kappa2_forecast)
rownames(kappa2_full) <- train_start:(train_end + n_years_forecast)

# Re-estimate collective drift using strucchange break detection on the
# annual national mean fitted LC log-mortality.
kappa_fr_original <- coherent_res$kappa_fr

# Build yearly national mean fitted log-mortality from LC fitted rates
train_dates <- as.Date(dimnames(dtxr)[[1]])
train_date_years <- lubridate::isoyear(train_dates)
logmhat_lc <- log(pmax(lca$mhat, 1e-12))
annual_years <- train_start:train_end
annual_lc_mean_logmx <- vapply(
    annual_years,
    function(y) {
        idx <- which(train_date_years == y)
        if (length(idx) == 0) {
            return(NA_real_)
        }
        mean(logmhat_lc[idx, , ], na.rm = TRUE)
    },
    numeric(1)
)

valid_year_mask <- is.finite(annual_lc_mean_logmx)
annual_years_valid <- annual_years[valid_year_mask]
annual_lc_valid <- annual_lc_mean_logmx[valid_year_mask]
time_idx <- seq_along(annual_lc_valid)

# Bai-Perron style break detection on first differences of annual LC log-mortality
# (breaks are detected on log mxt dynamics, not on kappa directly).
max_breaks_cap <- 5L
min_seg_abs <- 3L
min_seg_frac <- 0.05
# Sensitivity override: allow an m>=1 solution when its BIC is close to the m=0 minimum.
bic_force_break_margin <- 10.0
# Complexity preference: prefer the largest m among near-optimal BIC values.
bic_more_breaks_margin <- 8.0

delta_lc <- diff(annual_lc_valid)
delta_years <- annual_years_valid[-1]
n_delta <- length(delta_lc)

min_seg_years <- max(min_seg_abs, ceiling(min_seg_frac * n_delta))
max_breaks_allowed <- max(0L, floor(n_delta / min_seg_years) - 1L)
max_breaks_allowed <- min(max_breaks_allowed, max_breaks_cap)

candidate_break_years <- integer(0)
bic_values <- rep(NA_real_, max_breaks_allowed + 1L)
names(bic_values) <- paste0("m=", 0:max_breaks_allowed)
m_star <- 0L
break_diagnostics <- vector("list", max_breaks_allowed + 1L)

if (n_delta >= 2L && max_breaks_allowed >= 0L) {
    for (m in 0:max_breaks_allowed) {
        opt_m <- find_optimal_break_idx(
            series_vec = delta_lc,
            m = m,
            min_seg = min_seg_years
        )
        break_idx_m <- opt_m$break_idx
        break_years_m <- if (length(break_idx_m) > 0) {
            sort(unique(delta_years[break_idx_m]))
        } else {
            integer(0)
        }
        rss_m <- opt_m$rss
        last_break_m <- if (length(break_years_m) > 0) {
            max(break_years_m)
        } else {
            NA_integer_
        }
        if (!is.finite(rss_m) || rss_m <= 0) {
            break_diagnostics[[m + 1L]] <- data.frame(
                m = m,
                bic = NA_real_,
                break_years = if (length(break_years_m) > 0) {
                    paste(break_years_m, collapse = ",")
                } else {
                    "none"
                },
                last_break_year = last_break_m,
                stringsAsFactors = FALSE
            )
            next
        }
        # Piecewise-constant drift model on delta(log mxt): p = m + 1 means.
        p_m <- m + 1L
        bic_values[m + 1L] <- n_delta *
            log(rss_m / n_delta) +
            p_m * log(n_delta)
        break_diagnostics[[m + 1L]] <- data.frame(
            m = m,
            bic = bic_values[m + 1L],
            break_years = if (length(break_years_m) > 0) {
                paste(break_years_m, collapse = ",")
            } else {
                "none"
            },
            last_break_year = last_break_m,
            stringsAsFactors = FALSE
        )
    }

    valid_bic <- which(is.finite(bic_values))
    if (length(valid_bic) > 0) {
        best_idx <- valid_bic[which.min(bic_values[valid_bic])]
        m_star <- best_idx - 1L
        best_bic <- bic_values[best_idx]

        close_idx <- valid_bic[
            bic_values[valid_bic] <= (best_bic + bic_more_breaks_margin)
        ]
        if (length(close_idx) > 0) {
            m_close <- close_idx - 1L
            m_pref <- max(m_close, na.rm = TRUE)
            if (is.finite(m_pref) && m_pref > m_star) {
                m_star <- as.integer(m_pref)
            }
        }

        if (m_star == 0L && is.finite(bic_values[1])) {
            nonzero_valid <- which(is.finite(bic_values[-1])) + 1L
            if (length(nonzero_valid) > 0) {
                best_nonzero_idx <- nonzero_valid[
                    which.min(bic_values[nonzero_valid])
                ]
                bic_gap <- bic_values[best_nonzero_idx] - bic_values[1]
                if (is.finite(bic_gap) && bic_gap <= bic_force_break_margin) {
                    m_star <- best_nonzero_idx - 1L
                }
            }
        }

        if (m_star > 0L) {
            opt_star <- find_optimal_break_idx(
                series_vec = delta_lc,
                m = m_star,
                min_seg = min_seg_years
            )
            if (length(opt_star$break_idx) > 0) {
                candidate_break_years <- sort(unique(
                    delta_years[opt_star$break_idx]
                ))
            } else {
                candidate_break_years <- integer(0)
            }
        } else {
            candidate_break_years <- integer(0)
        }
    }
}

has_break <- length(candidate_break_years) > 0

drift_recent <- NA
change_point_years <- integer(0)

# Adjust the breakpoint comparison if the breakpoint year should be excluded from the recent segment.
if (has_break) {
    change_point_years <- candidate_break_years
    last_break_year <- max(change_point_years)

    kappa_years <- as.integer(rownames(kappa2_fitted))
    seg_indices <- which(kappa_years >= last_break_year)
    if (length(seg_indices) < 2) {
        seg_indices <- seq_along(kappa_fr_original)
    }

    recent_data <- kappa_fr_original[seg_indices]
    recent_time <- seq_along(recent_data)
    lm_recent <- lm(recent_data ~ recent_time)
    drift_recent <- coef(lm_recent)["recent_time"]

    dir.create("Results/plots", showWarnings = FALSE, recursive = TRUE)
    png("Results/plots/structural_break_trend.png", width = 1000, height = 600)
    plot(
        annual_years_valid,
        annual_lc_valid,
        type = "l",
        lwd = 2,
        main = "Yearly National LC Mean Log-Mortality with Structural Breaks",
        xlab = "Year",
        ylab = "Mean log-mortality (LC fitted)"
    )
    grid()
    abline(v = change_point_years, col = "red", lty = 2)
    seg_factor <- cut(
        annual_years_valid,
        breaks = c(-Inf, change_point_years, Inf),
        labels = FALSE,
        right = TRUE
    )
    fitted_segments <- tryCatch(
        fitted(lm(annual_lc_valid ~ time_idx * as.factor(seg_factor))),
        error = function(e) NULL
    )
    if (!is.null(fitted_segments)) {
        lines(
            annual_years_valid,
            fitted_segments,
            col = "blue",
            lwd = 2,
            lty = 2
        )
    }
    legend(
        "topright",
        legend = c(
            "Annual LC mean log-mortality",
            "Breakpoint",
            "Fitted Segments"
        ),
        col = c("black", "red", "blue"),
        lty = c(1, 2, 2),
        lwd = 2
    )
    dev.off()
} else {
    drift_recent <- mean(diff(kappa_fr_original), na.rm = TRUE)
}

last_kappa_fr <- tail(kappa_fr_original, 1)
kappa_fr_forecast_recent <- last_kappa_fr +
    seq_len(n_years_forecast) * drift_recent

kappa2_forecast_recent <- sweep(
    coherent_res$u_forecast,
    1,
    kappa_fr_forecast_recent,
    "+"
)
colnames(kappa2_forecast_recent) <- colnames(kappa2_fitted)
rownames(kappa2_forecast_recent) <- as.character(
    (train_end + 1):(train_end + n_years_forecast)
)

kappa2_forecast <- kappa2_forecast_recent
kappa2_full <- rbind(kappa2_fitted, kappa2_forecast)
rownames(kappa2_full) <- as.character(
    train_start:(train_end + n_years_forecast)
)
forecast_years <- as.integer(rownames(kappa2_full))

coherent_res$kappa_fr_forecast <- kappa_fr_forecast_recent
coherent_res$drift_recent <- drift_recent

################################################################################
# COMPUTE BASELINE LOG MORTALITY
################################################################################
# Generate forecasts for test period
DfM_test <- DfM_test %>%
    mutate(
        age_idx = match(age, xv),
        region_idx = match(region, dimnames(dtxr)[[3]]),
        logmxt_baseline = mapply(
            predict_logmxt,
            year = ISOYear,
            week = ISOWeek,
            age_idx = age_idx,
            region_idx = region_idx,
            MoreArgs = list(
                beta1 = lca$beta1,
                beta2 = lca$beta2,
                beta3 = lca$beta3,
                kappa2_full = kappa2_full,
                kappa3 = lca$kappa3
            )
        )
    )

# Get dimensions
n_weeks_train <- nrow(dtxr)
n_ages <- length(xv)
n_regions <- R
region_names <- dimnames(dtxr)[[3]]
age_names <- as.character(xv)

# Create fitted log mortality for training period
logmxt_train <- array(
    NA,
    dim = c(n_weeks_train, n_ages, n_regions),
    dimnames = list(rownames(dtxr), age_names, region_names)
)

for (t in 1:n_weeks_train) {
    year_idx <- isoyv[t]
    week_idx <- isowv[t]
    for (a in 1:n_ages) {
        for (r in 1:n_regions) {
            logmxt_train[t, a, r] <- lca$beta1[a, r] +
                lca$beta2[a] * lca$kappa2[year_idx, r] +
                lca$beta3[a] * lca$kappa3[week_idx, r]
        }
    }
}

# Compute forecasted log mortality for test period
test_dates <- sort(unique(DfM_test$ISODate))
n_weeks_test <- length(test_dates)

logmxt_test <- array(
    NA,
    dim = c(n_weeks_test, n_ages, n_regions),
    dimnames = list(as.character(test_dates), age_names, region_names)
)

for (t in 1:n_weeks_test) {
    date_t <- test_dates[t]
    year_t <- isoyear(date_t)
    week_t <- isoweek(date_t)
    year_idx <- resolve_year_index(
        year_t,
        forecast_years,
        "forecast LC year index"
    )
    for (a in 1:n_ages) {
        for (r in 1:n_regions) {
            logmxt_test[t, a, r] <- lca$beta1[a, r] +
                lca$beta2[a] * kappa2_full[year_idx, r] +
                lca$beta3[a] * lca$kappa3[week_t, r]
        }
    }
}

# Combine into a single array
logmxt_combined <- abind::abind(logmxt_train, logmxt_test, along = 1)

# Create long-format data frame
logmxt_combined_long <- expand.grid(
    ISODate = c(as.Date(rownames(dtxr)), test_dates),
    age = xv,
    region = region_names,
    stringsAsFactors = FALSE
) %>%
    mutate(
        ISOYear = isoyear(ISODate),
        ISOWeek = isoweek(ISODate),
        logmxt = as.vector(logmxt_combined),
        mxt = exp(logmxt),
        period = ifelse(ISOYear <= train_end, "Fitted", "Forecast")
    ) %>%
    arrange(ISODate, age, region)

# Create wide-format matrices for train and test periods
baseline_logmxt_wide_train <- logmxt_combined_long %>%
    dplyr::filter(period == "Fitted") %>%
    dplyr::mutate(age = as.character(age)) %>%
    dplyr::select(ISODate, age, region, logmxt) %>%
    tidyr::pivot_wider(
        id_cols = ISODate,
        names_from = c(age, region),
        values_from = logmxt,
        names_sep = "_"
    ) %>%
    dplyr::arrange(ISODate)

baseline_logmxt_wide <- logmxt_combined_long %>%
    dplyr::filter(period == "Forecast") %>%
    dplyr::mutate(age = as.character(age)) %>%
    dplyr::select(ISODate, age, region, logmxt) %>%
    tidyr::pivot_wider(
        id_cols = ISODate,
        names_from = c(age, region),
        values_from = logmxt,
        names_sep = "_"
    ) %>%
    dplyr::arrange(ISODate)

names(baseline_logmxt_wide_train)[
    names(baseline_logmxt_wide_train) == "ISODate"
] <- "Week"
names(baseline_logmxt_wide)[names(baseline_logmxt_wide) == "ISODate"] <- "Week"

# Preserve raw baseline tables for MortFCNet before subsequent age-bias correction.
baseline_logmxt_wide_train_raw <- baseline_logmxt_wide_train
baseline_logmxt_wide_raw <- baseline_logmxt_wide

# Create observed log mortality wide format
observed_logmxt_wide <- DfM_all %>%
    dplyr::mutate(
        logmxt = log(deaths / expo),
        age = as.character(age)
    ) %>%
    dplyr::mutate(logmxt = ifelse(!is.finite(logmxt), NA, logmxt)) %>%
    dplyr::select(ISODate, age, region, logmxt) %>%
    tidyr::pivot_wider(
        id_cols = ISODate,
        names_from = c(age, region),
        values_from = logmxt,
        names_sep = "_"
    ) %>%
    dplyr::arrange(ISODate) %>%
    dplyr::rename(Week = ISODate)

# Create exposure (population-at-risk) wide format for weighted metrics.
exposure_wide <- DfM_all %>%
    dplyr::mutate(age = as.character(age)) %>%
    dplyr::select(ISODate, age, region, expo) %>%
    tidyr::pivot_wider(
        id_cols = ISODate,
        names_from = c(age, region),
        values_from = expo,
        names_sep = "_"
    ) %>%
    dplyr::arrange(ISODate) %>%
    dplyr::rename(Week = ISODate)

################################################################################
# ADJACENCY MATRIX FOR GNN
################################################################################
shapefile <- read_sf('Data/NUTS/NUTS_RG_20M_2021_3035.shp')
nuts.spec <- 2
ctry.spec <- 'FR'
shapef <- shapefile[
    shapefile$CNTR_CODE %in% ctry.spec & shapefile$LEVL_CODE == nuts.spec,
]

overseas <- c(paste0('FRY', 1:5), 'FRM0')
ind.rm <- unlist(sapply(overseas, function(x) {
    which(grepl(x, shapef$NUTS_ID, fixed = TRUE))
}))
shapef <- shapef[-c(ind.rm), ] %>%
    dplyr::arrange(NUTS_ID)

shapefile_regions <- sort(unique(shapef$NUTS_ID))

adjacency_file <- 'Data/NUTS/nuts2_adjacency_FR.rds'
adj_matrix_raw <- readRDS(adjacency_file)
adj_matrix_raw <- adj_matrix_raw[shapefile_regions, shapefile_regions]

# ==============================================================================
# PHASE 4: CLIMATE COVARIATE
# ==============================================================================
climate_df <- read_eobs_files(
    path = "Data/E-OBS",
    country = "FR",
    nuts = "NUTS2",
    combine = TRUE
)

col_clim <- "Region"

climate_df <- climate_df[
    !is.na(climate_df[[col_clim]]) &
        trimws(as.character(climate_df[[col_clim]])) != "FRM0",
]

names(climate_df)[names(climate_df) == col_clim] <- "Region"

climate_df <- climate_df %>%
    dplyr::mutate(
        Week = floor_date(Date, unit = "week", week_start = 1),
        Country_Code = substr(Region, 1, 2)
    )

################################################################################
# AGGREGATE DAILY to WEEKLY
################################################################################
measurement_cols <- c("tn", "tg", "tx", "rr", "hu", "fg")
min_vars <- intersect(measurement_cols, "tn")
max_vars <- intersect(measurement_cols, c("tx", "fg"))
sum_vars <- intersect(measurement_cols, "rr")
mean_vars <- setdiff(measurement_cols, c(min_vars, max_vars, sum_vars))

climate_weekly_region <- climate_df %>%
    dplyr::group_by(Week, Region, Country_Code) %>%
    dplyr::summarise(
        dplyr::across(tidyselect::any_of(min_vars), ~ min(.x, na.rm = TRUE)),
        dplyr::across(tidyselect::any_of(max_vars), ~ max(.x, na.rm = TRUE)),
        dplyr::across(tidyselect::any_of(sum_vars), ~ sum(.x, na.rm = TRUE)),
        dplyr::across(tidyselect::any_of(mean_vars), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
    ) %>%
    dplyr::arrange(Region, Week) %>%
    dplyr::filter(Week >= as.Date("1990-01-01"))

# Replace infinite values produced by aggregating empty groups
climate_weekly_region[measurement_cols][
    !is.finite(as.matrix(climate_weekly_region[measurement_cols]))
] <- 0

climate_weekly_region_z <- climate_weekly_region %>%
    dplyr::mutate(week_of_year = lubridate::isoweek(Week)) %>%
    dplyr::group_by(Region, week_of_year) %>%
    dplyr::mutate(
        dplyr::across(
            all_of(measurement_cols),
            ~ (.x -
                mean(
                    .x[lubridate::isoyear(Week) %in% fit_years],
                    na.rm = TRUE
                )) /
                sd(.x[lubridate::isoyear(Week) %in% fit_years], na.rm = TRUE)
        )
    ) %>%
    dplyr::ungroup()

# Replace infinite values introduced by zero variance in a training-week group.
climate_weekly_region_z[measurement_cols][
    !is.finite(as.matrix(climate_weekly_region_z[measurement_cols]))
] <- 0

################################################################################
# CALCULATE RESIDUALS
################################################################################
obs_mxt_train <- stats::xtabs(
    deaths ~ ISODate + age + region,
    data = DfM_train
) /
    stats::xtabs(expo ~ ISODate + age + region, data = DfM_train)
obs_mxt_train[obs_mxt_train <= 0] <- 1e-10
baseline_mxt_train <- exp(logmxt_train)
additive_res_train <- obs_mxt_train - baseline_mxt_train
additive_res_train[!is.finite(additive_res_train)] <- NA

obs_mxt_test <- stats::xtabs(deaths ~ ISODate + age + region, data = DfM_test) /
    stats::xtabs(expo ~ ISODate + age + region, data = DfM_test)
obs_mxt_test[obs_mxt_test <= 0] <- 1e-10
baseline_mxt_test <- exp(logmxt_test)
additive_res_test <- obs_mxt_test - baseline_mxt_test
additive_res_test[!is.finite(additive_res_test)] <- NA

additive_res <- abind::abind(additive_res_train, additive_res_test, along = 1)
age_bias <- apply(additive_res_train, 2, mean, na.rm = TRUE)
names(age_bias) <- age_names

for (a_idx in seq_along(age_names)) {
    additive_res_train[, a_idx, ] <- additive_res_train[, a_idx, ] -
        age_bias[a_idx]
    additive_res_test[, a_idx, ] <- additive_res_test[, a_idx, ] -
        age_bias[a_idx]
}
additive_res <- abind::abind(additive_res_train, additive_res_test, along = 1)

for (a_chr in age_names) {
    b <- age_bias[a_chr]
    cols_tr <- grep(
        paste0("^", a_chr, "_"),
        names(baseline_logmxt_wide_train),
        value = TRUE
    )
    cols_te <- grep(
        paste0("^", a_chr, "_"),
        names(baseline_logmxt_wide),
        value = TRUE
    )
    if (length(cols_tr) > 0) {
        baseline_logmxt_wide_train[, cols_tr] <- log(
            pmax(
                exp(baseline_logmxt_wide_train[, cols_tr, drop = FALSE]) + b,
                1e-12
            )
        )
    }
    if (length(cols_te) > 0) {
        baseline_logmxt_wide[, cols_te] <- log(
            pmax(
                exp(baseline_logmxt_wide[, cols_te, drop = FALSE]) + b,
                1e-12
            )
        )
    }
}

# ==============================================================================
# PHASE 5: DEEP LEARNING DATASET (TARGETS & SEQUENCES)
# ==============================================================================
residual_long <- expand.grid(
    Week = as.Date(dimnames(additive_res)[[1]]),
    AgeNumeric = dimnames(additive_res)[[2]],
    Region = dimnames(additive_res)[[3]],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
residual_long$residual <- as.vector(additive_res)

residual_wide_region <- residual_long %>%
    tidyr::pivot_wider(
        id_cols = c(Week, Region),
        names_from = AgeNumeric,
        values_from = residual,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Region, Week)

target_cols <- setdiff(names(residual_wide_region), c("Week", "Region"))

residual_wide_region_lag <- residual_wide_region %>%
    dplyr::arrange(Region, Week)
residual_wide_region_lag$mean_residual <- rowMeans(
    residual_wide_region_lag[, target_cols, drop = FALSE],
    na.rm = TRUE
)

lagged_res_wide <- residual_wide_region_lag %>%
    dplyr::group_by(Region) %>%
    dplyr::arrange(Week, .by_group = TRUE) %>%
    dplyr::mutate(
        mean_res_lag1 = dplyr::lag(mean_residual, 1L),
        mean_res_lag2 = dplyr::lag(mean_residual, 2L),
        mean_res_lag3 = dplyr::lag(mean_residual, 3L)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
        Week,
        Region,
        mean_res_lag1,
        mean_res_lag2,
        mean_res_lag3
    )

lagged_age_res_wide <- residual_wide_region_lag %>%
    dplyr::group_by(Region) %>%
    dplyr::arrange(Week, .by_group = TRUE) %>%
    dplyr::mutate(
        dplyr::across(
            all_of(target_cols),
            list(
                lag1 = ~ dplyr::lag(.x, 1L),
                lag2 = ~ dplyr::lag(.x, 2L),
                lag3 = ~ dplyr::lag(.x, 3L)
            ),
            .names = "age{.col}_{.fn}"
        )
    ) %>%
    dplyr::ungroup()

age_lag_feature_cols <- grep(
    "^age.*_lag[123]$",
    names(lagged_age_res_wide),
    value = TRUE
)

lagged_age_res_wide <- lagged_age_res_wide %>%
    dplyr::select(
        Week,
        Region,
        dplyr::all_of(age_lag_feature_cols)
    )

df_model <- dplyr::inner_join(
    climate_weekly_region_z,
    residual_wide_region,
    by = c("Week", "Region")
) %>%
    dplyr::left_join(lagged_res_wide, by = c("Week", "Region")) %>%
    dplyr::left_join(lagged_age_res_wide, by = c("Week", "Region")) %>%
    dplyr::arrange(Region, Week)

lag_feature_cols <- c(
    "mean_res_lag1",
    "mean_res_lag2",
    "mean_res_lag3",
    age_lag_feature_cols
)
lag_fit_rows <- lubridate::isoyear(df_model$Week) %in% fit_years
for (lag_col in lag_feature_cols) {
    lag_center <- mean(df_model[[lag_col]][lag_fit_rows], na.rm = TRUE)
    lag_scale <- sd(df_model[[lag_col]][lag_fit_rows], na.rm = TRUE)
    if (!is.finite(lag_scale) || lag_scale == 0) {
        lag_scale <- 1
    }
    df_model[[lag_col]] <- (df_model[[lag_col]] - lag_center) / lag_scale
    df_model[[lag_col]][!is.finite(df_model[[lag_col]])] <- 0
}

X_num <- as.matrix(
    df_model[, c(measurement_cols, lag_feature_cols), drop = FALSE]
)

# Add ISOYear as a standardised covariate to capture long-run mortality trend
year_vals <- lubridate::isoyear(df_model$Week)
year_center_val <- mean(year_vals[year_vals %in% fit_years])
year_scale_val <- sd(year_vals[year_vals %in% fit_years])
if (!is.finite(year_scale_val) || year_scale_val == 0) {
    year_scale_val <- 1
}
X_num <- cbind(
    X_num,
    ISOYear_z = (year_vals - year_center_val) / year_scale_val
)
# GNN gets the year signal through both the climate features and
# the explicit `year_numeric_input`.
X_num_gnn <- X_num[, colnames(X_num) != "ISOYear_z", drop = FALSE]
# MortFCNet should not receive the lagged mortality features.
# That includes both the region-mean lag features and the age-specific lag
# features added for CNN-LSTM and GNN-LSTM.
X_num_mort <- X_num[, !colnames(X_num) %in% lag_feature_cols, drop = FALSE]
n_climate_features_mort <- ncol(X_num_mort)

y_mat <- as.matrix(df_model[, target_cols, drop = FALSE])

################################################################################
# STANDARDISE TARGETS
################################################################################
df_model_years <- lubridate::isoyear(df_model$Week)
fit_rows_df <- df_model_years %in% fit_years

# Standardise using fit-period statistics only (exclude validation years)
y_center <- colMeans(y_mat[fit_rows_df, , drop = FALSE], na.rm = TRUE)
# Pool target scaling across regions within each age group so high-age
# mortality variance is not compressed by low-variance region columns.
age_numeric_cols <- as.integer(sub("_.*", "", target_cols))
y_scale_raw <- apply(
    y_mat[fit_rows_df, , drop = FALSE],
    2,
    stats::sd,
    na.rm = TRUE
)
y_scale_pooled <- tapply(y_scale_raw, age_numeric_cols, mean)
y_scale <- y_scale_pooled[as.character(age_numeric_cols)]
names(y_scale) <- target_cols
y_scale[y_scale == 0 | is.na(y_scale)] <- 1
y_mat_z <- scale(y_mat, center = y_center, scale = y_scale)

# Pre-center standardized targets so models learn climate-driven deviation
# rather than a static age-region offset.
col_bias_z <- colMeans(y_mat_z[fit_rows_df, , drop = FALSE], na.rm = TRUE)
for (j in seq_along(target_cols)) {
    y_mat_z[, j] <- y_mat_z[, j] - col_bias_z[j]
}

# Keep the same per-column centering bias in original residual units for
# inference-time reconstruction after unscaling model outputs.
col_bias_original <- col_bias_z * y_scale
names(col_bias_original) <- target_cols
age_group_bias_original <- col_bias_original

y_valid_mask <- is.finite(y_mat_z)
y_mat_z[!y_valid_mask] <- 0

# Mask out weeks that were imputed (so sequences don't include imputed rows)
imputed_week_dates <- dates_dtxr[ind.rm]
imputed_week_rows <- df_model$Week %in% imputed_week_dates
if (any(imputed_week_rows)) {
    y_valid_mask[imputed_week_rows, ] <- FALSE
    y_mat_z[imputed_week_rows, ] <- 0
}

# MortFCNet should learn a single residual per region-week so the climate
# effect is shared across all age groups.
mort_common_residual_raw <- rowMeans(y_mat, na.rm = TRUE)
mort_common_residual_raw[!is.finite(mort_common_residual_raw)] <- NA_real_
y_center_mort_common <- mean(
    mort_common_residual_raw[fit_rows_df],
    na.rm = TRUE
)
y_scale_mort_common <- sd(mort_common_residual_raw[fit_rows_df], na.rm = TRUE)
if (!is.finite(y_scale_mort_common) || y_scale_mort_common == 0) {
    y_scale_mort_common <- 1
}
y_mat_mort_common_z <- matrix(
    (mort_common_residual_raw - y_center_mort_common) / y_scale_mort_common,
    ncol = 1,
    dimnames = list(NULL, "mort_common")
)
mort_common_valid <- is.finite(mort_common_residual_raw)
if (any(imputed_week_rows)) {
    mort_common_valid[imputed_week_rows] <- FALSE
}
y_mat_mort_common_z[!mort_common_valid, 1] <- 0

baseline_logmxt_wide_train_mort <- baseline_logmxt_wide_train_raw
baseline_logmxt_wide_mort <- baseline_logmxt_wide_raw

y_mat_mort <- restore_age_bias_to_targets(y_mat, target_cols, age_bias)

y_center_mort <- colMeans(y_mat_mort[fit_rows_df, , drop = FALSE], na.rm = TRUE)
y_scale_raw_mort <- apply(
    y_mat_mort[fit_rows_df, , drop = FALSE],
    2,
    stats::sd,
    na.rm = TRUE
)
y_scale_pooled_mort <- tapply(y_scale_raw_mort, age_numeric_cols, mean)
y_scale_mort <- y_scale_pooled_mort[as.character(age_numeric_cols)]
names(y_scale_mort) <- target_cols
y_scale_mort[y_scale_mort == 0 | is.na(y_scale_mort)] <- 1
y_mat_mort_z <- scale(y_mat_mort, center = y_center_mort, scale = y_scale_mort)

# Recenter standardized targets for MortFCNet so its residuals remain zero-centered.
col_bias_z_mort <- colMeans(
    y_mat_mort_z[fit_rows_df, , drop = FALSE],
    na.rm = TRUE
)
for (j in seq_along(target_cols)) {
    y_mat_mort_z[, j] <- y_mat_mort_z[, j] - col_bias_z_mort[j]
}

col_bias_original_mort <- col_bias_z_mort * y_scale_mort
names(col_bias_original_mort) <- target_cols
age_group_bias_original_mort <- col_bias_original_mort

y_valid_mask_mort <- is.finite(y_mat_mort_z)
y_mat_mort_z[!y_valid_mask_mort] <- 0
if (any(imputed_week_rows)) {
    y_valid_mask_mort[imputed_week_rows, ] <- FALSE
    y_mat_mort_z[imputed_week_rows, ] <- 0
}

# Use the same age-bias-restored target baseline for all models so CNN/LSTM,
# GNN/LSTM and MortFCNet are trained on the same target signal.
baseline_logmxt_wide_train <- baseline_logmxt_wide_train_raw
baseline_logmxt_wide <- baseline_logmxt_wide_raw
y_scale <- y_scale_mort
y_center <- y_center_mort
col_bias_original <- col_bias_original_mort
age_group_bias_original <- age_group_bias_original_mort
y_valid_mask <- y_valid_mask_mort

################################################################################
# PREPARE SEQUENCES (CNN-LSTM + GNN-LSTM)
################################################################################
region_levels <- sort(unique(df_model$Region))
if (!setequal(region_levels, shapefile_regions)) {
    stop(
        "Region identifiers differ between the mortality/climate data and the adjacency matrix."
    )
}
regions <- region_levels
adj_matrix <- adj_matrix_raw[regions, regions]
adj_matrix_norm <- adj_matrix + diag(nrow(adj_matrix))
degree_matrix <- diag(1 / sqrt(rowSums(adj_matrix_norm)))
adj_matrix_norm <- degree_matrix %*% adj_matrix_norm %*% degree_matrix
gnn_adjacency_alpha <- 0.7
adj_matrix_norm_soft <- gnn_adjacency_alpha *
    adj_matrix_norm +
    (1 - gnn_adjacency_alpha) * diag(nrow(adj_matrix_norm))
n_climate_features <- ncol(X_num)

seqs_age_aware <- make_sequences_age_aware(
    df_key = df_model[, c("Week", "Region"), drop = FALSE],
    X_climate = X_num,
    y_mat = y_mat_mort_z,
    seq_length = seq_length_cnn_gnn,
    region_levels = region_levels,
    age_cols = target_cols
)

target_years_aa <- lubridate::isoyear(seqs_age_aware$target_week)
train_mask_aa <- target_years_aa %in% train_years
test_mask_aa <- target_years_aa %in% test_years

X_train_aa <- seqs_age_aware$X[train_mask_aa, , , drop = FALSE]
y_train_aa <- seqs_age_aware$y[train_mask_aa]
region_train_aa <- matrix(seqs_age_aware$region_idx[train_mask_aa], ncol = 1)
age_train_aa <- matrix(seqs_age_aware$age_idx[train_mask_aa], ncol = 1)
age_mid_train_aa <- matrix(seqs_age_aware$age_mid[train_mask_aa], ncol = 1)
X_test_aa <- seqs_age_aware$X[test_mask_aa, , , drop = FALSE]
region_test_aa <- matrix(seqs_age_aware$region_idx[test_mask_aa], ncol = 1)
age_test_aa <- matrix(seqs_age_aware$age_idx[test_mask_aa], ncol = 1)
age_mid_test_aa <- matrix(seqs_age_aware$age_mid[test_mask_aa], ncol = 1)

# Drop samples where the target is a structural zero (no mortality signal)
df_key_idx_aa <- paste(df_model$Region, format(df_model$Week), sep = "|")
seqs_key_aa <- paste(
    seqs_age_aware$Region,
    format(seqs_age_aware$target_week),
    sep = "|"
)
seqs_df_row_aa <- match(seqs_key_aa, df_key_idx_aa)
seqs_age_col_aa <- match(seqs_age_aware$Age, target_cols)
y_valid_seqs_aa <- y_valid_mask[cbind(seqs_df_row_aa, seqs_age_col_aa)]
y_valid_seqs_aa[is.na(y_valid_seqs_aa)] <- FALSE

keep_train <- which(y_valid_seqs_aa[train_mask_aa])
X_train_aa <- X_train_aa[keep_train, , , drop = FALSE]
y_train_aa <- y_train_aa[keep_train]
region_train_aa <- region_train_aa[keep_train, , drop = FALSE]
age_train_aa <- age_train_aa[keep_train, , drop = FALSE]
age_mid_train_aa <- age_mid_train_aa[keep_train, , drop = FALSE]
valid_train_mask_aa <- which(train_mask_aa)[keep_train]

# Train on the full training window with fixed epochs.
fit_idx <- seq_len(length(y_train_aa))

# Standardise numeric age midpoints on the full training set.
age_mid_center <- mean(age_mid_train_aa[fit_idx], na.rm = TRUE)
age_mid_scale <- sd(age_mid_train_aa[fit_idx], na.rm = TRUE)
if (!is.finite(age_mid_scale) || age_mid_scale == 0) {
    age_mid_scale <- 1
}
age_mid_train_aa <- (age_mid_train_aa - age_mid_center) / age_mid_scale
age_mid_test_aa <- (age_mid_test_aa - age_mid_center) / age_mid_scale

X_fit_aa <- X_train_aa
y_fit_aa <- y_train_aa
region_fit_aa <- region_train_aa
age_fit_aa <- age_train_aa
age_mid_fit_aa <- age_mid_train_aa

# Weight recent training samples more heavily and upweight peak residual weeks.
sample_weights_fit <- make_peak_aware_sample_weights(
    seqs_age_aware$target_week[train_mask_aa][keep_train],
    observed_residuals = y_train_aa,
    recency_strength = 0.0,
    peak_strength = 1,
    peak_quantile = 0.75
)

################################################################################
# MortFCNet SEQUENCES (dedicated shorter context window)
################################################################################
seqs_age_aware_mort <- make_sequences_age_aware(
    df_key = df_model[, c("Week", "Region"), drop = FALSE],
    X_climate = X_num_mort,
    y_mat = y_mat_mort_common_z,
    seq_length = seq_length_mort,
    region_levels = region_levels,
    age_cols = c("mort_common")
)

stopifnot(all(
    lubridate::isoyear(seqs_age_aware_mort$target_week) %in%
        c(train_years, test_years)
))

target_years_mort <- lubridate::isoyear(seqs_age_aware_mort$target_week)
train_mask_mort <- target_years_mort %in% train_years
test_mask_mort <- target_years_mort %in% test_years

X_train_mort <- seqs_age_aware_mort$X[train_mask_mort, , , drop = FALSE]
y_train_mort <- seqs_age_aware_mort$y[train_mask_mort]
X_test_mort <- seqs_age_aware_mort$X[test_mask_mort, , , drop = FALSE]

# Drop samples where the aggregated target is a structural zero
df_key_idx_mort <- paste(df_model$Region, format(df_model$Week), sep = "|")
seqs_key_mort <- paste(
    seqs_age_aware_mort$Region,
    format(seqs_age_aware_mort$target_week),
    sep = "|"
)
seqs_df_row_mort <- match(seqs_key_mort, df_key_idx_mort)
mort_common_valid_seq <- mort_common_valid[seqs_df_row_mort]
mort_common_valid_seq[is.na(mort_common_valid_seq)] <- FALSE

keep_train_mort <- which(mort_common_valid_seq[train_mask_mort])
X_train_mort <- X_train_mort[keep_train_mort, , , drop = FALSE]
y_train_mort <- y_train_mort[keep_train_mort]
valid_train_mask_mort <- which(train_mask_mort)[keep_train_mort]

# Train on the full training window with fixed epochs.

X_fit_mort <- X_train_mort
y_fit_mort <- y_train_mort

# MortFCNet uses unweighted training for the common-region residual target.

################################################################################
# GNN-LSTM SEQUENCES
################################################################################
seqs_gnn <- make_sequences_gnn_scalar(
    df_key = df_model[, c("Week", "Region"), drop = FALSE],
    X_mat = X_num_gnn,
    y_mat = y_mat_mort_z,
    seq_length = seq_length_cnn_gnn,
    regions_list = regions,
    age_cols = target_cols
)

target_years_gnn <- lubridate::isoyear(seqs_gnn$target_week)
train_mask_gnn <- target_years_gnn %in% train_years
test_mask_gnn <- target_years_gnn %in% test_years
X_train_gnn <- seqs_gnn$X[train_mask_gnn, , , drop = FALSE]
X_test_gnn <- seqs_gnn$X[test_mask_gnn, , , drop = FALSE]
y_train_gnn <- seqs_gnn$y[train_mask_gnn]
age_idx_train_gnn <- matrix(seqs_gnn$age_idx[train_mask_gnn], ncol = 1)
age_idx_test_gnn <- matrix(seqs_gnn$age_idx[test_mask_gnn], ncol = 1)
region_idx_train_gnn <- matrix(seqs_gnn$region_idx[train_mask_gnn], ncol = 1)
region_idx_test_gnn <- matrix(seqs_gnn$region_idx[test_mask_gnn], ncol = 1)
age_name_train_gnn <- seqs_gnn$Age[train_mask_gnn]
age_name_test_gnn <- seqs_gnn$Age[test_mask_gnn]
region_name_train_gnn <- seqs_gnn$Region[train_mask_gnn]
region_name_test_gnn <- seqs_gnn$Region[test_mask_gnn]
target_week_train_gnn <- seqs_gnn$target_week[train_mask_gnn]
target_week_test_gnn <- seqs_gnn$target_week[test_mask_gnn]

year_numeric_train_gnn <- matrix(
    (target_years_gnn[train_mask_gnn] - year_center_val) / year_scale_val,
    ncol = 1
)

df_key_idx_gnn <- paste(df_model$Region, format(df_model$Week), sep = "|")
seqs_key_gnn <- paste(
    seqs_gnn$Region,
    format(seqs_gnn$target_week),
    sep = "|"
)
seqs_df_row_gnn <- match(seqs_key_gnn, df_key_idx_gnn)
seqs_age_col_gnn <- match(seqs_gnn$Age, target_cols)
y_valid_seqs_gnn <- y_valid_mask[cbind(seqs_df_row_gnn, seqs_age_col_gnn)]
y_valid_seqs_gnn[is.na(y_valid_seqs_gnn)] <- FALSE

keep_train_gnn <- which(y_valid_seqs_gnn[train_mask_gnn])
X_train_gnn <- X_train_gnn[keep_train_gnn, , , drop = FALSE]
y_train_gnn <- y_train_gnn[keep_train_gnn]
age_idx_train_gnn <- age_idx_train_gnn[keep_train_gnn, , drop = FALSE]
region_idx_train_gnn <- region_idx_train_gnn[keep_train_gnn, , drop = FALSE]
age_name_train_gnn <- age_name_train_gnn[keep_train_gnn]
region_name_train_gnn <- region_name_train_gnn[keep_train_gnn]
target_week_train_gnn <- target_week_train_gnn[keep_train_gnn]
year_numeric_train_gnn <- year_numeric_train_gnn[
    keep_train_gnn,
    ,
    drop = FALSE
]

year_numeric_test_gnn <- matrix(
    (target_years_gnn[test_mask_gnn] - year_center_val) / year_scale_val,
    ncol = 1
)

all_results <- list()

# ==============================================================================
# PHASE 6: NEURAL NETWORK TRAINING (CNN-LSTM, MortFCNet, GNN-LSTM)
# ==============================================================================
# MODEL 1: CNN-LSTM
start_time_cnn <- Sys.time()

# Build args for cnn_lstm_age_aware.
cnn_args <- list(
    seq_length = seq_length_cnn_gnn,
    n_climate_features = n_climate_features,
    n_regions = length(region_levels),
    n_ages = length(xv),
    cnn_filters = best_cnn$cnn_filters,
    lstm_units = best_cnn$lstm_units,
    hidden_size_1 = best_cnn$hidden_size_1,
    hidden_size_2 = best_cnn$hidden_size_2,
    hidden_size_3 = best_cnn$hidden_size_3,
    region_embedding_dim = best_cnn$region_embedding_dim,
    age_embedding_dim = best_cnn$age_embedding_dim,
    dropout_rate = best_cnn$dropout_rate,
    cell_type = "lstm",
    use_age_embedding = use_age_embedding_cnn,
    use_region_embedding = use_region_embedding_cnn,
    use_age_numeric = TRUE
)

mdl_cnn <- do.call(cnn_lstm_age_aware, cnn_args)

opt_cnn <- keras::optimizer_adam(
    learning_rate = best_cnn$learning_rate,
    clipnorm = 0.5
)
mdl_cnn$compile(
    optimizer = opt_cnn,
    loss = "mse",
    metrics = list("mse")
)

x_fit_list <- select_model_inputs(
    mdl_cnn,
    list(
        climate_input = X_fit_aa,
        region_input = region_fit_aa,
        age_input = age_fit_aa,
        age_numeric_input = age_mid_fit_aa
    )
)

history_cnn <- keras::fit(
    mdl_cnn,
    x = x_fit_list,
    y = y_fit_aa,
    sample_weight = sample_weights_fit,
    shuffle = TRUE,
    epochs = epochs_final_cnn,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_cnn,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

end_time_cnn <- Sys.time()
time_cnn <- difftime(end_time_cnn, start_time_cnn, units = "mins")
epochs_cnn <- length(history_cnn$metrics$loss)

model_timings[["CNN-LSTM"]] <- list(
    time_minutes = as.numeric(time_cnn),
    epochs_to_converge = epochs_cnn,
    final_train_loss = tail(history_cnn$metrics$loss, 1)
)

################################################################################
# INTERCEPT CORRECTION - PER AGE GROUP
################################################################################
X_train_all_ic <- seqs_age_aware$X[valid_train_mask_aa, , , drop = FALSE]
region_train_all_ic <- matrix(
    seqs_age_aware$region_idx[valid_train_mask_aa],
    ncol = 1
)
age_train_all_ic <- matrix(
    seqs_age_aware$age_idx[valid_train_mask_aa],
    ncol = 1
)
age_train_all_ic_mid <- matrix(
    (seqs_age_aware$age_mid[valid_train_mask_aa] - age_mid_center) /
        age_mid_scale,
    ncol = 1
)

pred_ic_list <- select_model_inputs(
    mdl_cnn,
    list(
        climate_input = X_train_all_ic,
        region_input = region_train_all_ic,
        age_input = age_train_all_ic,
        age_numeric_input = age_train_all_ic_mid
    )
)
pred_train_cnn_z <- as.vector(predict(mdl_cnn, pred_ic_list))
age_idx_train_full <- seqs_age_aware$age_idx[valid_train_mask_aa]

intercept_fix_per_age <- col_bias_original

# Add back pre-centered age bias in original residual space.
age_idx_test_aa <- seqs_age_aware$age_idx[test_mask_aa]

col_names_train_full <- target_cols[age_idx_train_full]
col_names_test_aa <- target_cols[age_idx_test_aa]

stopifnot(all(col_names_train_full %in% target_cols))
stopifnot(all(col_names_test_aa %in% target_cols))

pred_test_list <- select_model_inputs(
    mdl_cnn,
    list(
        climate_input = X_test_aa,
        region_input = region_test_aa,
        age_input = age_test_aa,
        age_numeric_input = age_mid_test_aa
    )
)
pred_test_cnn_z <- as.vector(predict(mdl_cnn, pred_test_list))

pred_train_cnn <- pred_train_cnn_z *
    y_scale[col_names_train_full] +
    y_center[col_names_train_full] +
    intercept_fix_per_age[col_names_train_full]
pred_test_cnn <- pred_test_cnn_z *
    y_scale[col_names_test_aa] +
    y_center[col_names_test_aa] +
    intercept_fix_per_age[col_names_test_aa]

################################################################################
# BUILD WIDE PREDICTION TABLES — CNN-LSTM
################################################################################
test_pred_df_cnn <- data.frame(
    Week = seqs_age_aware$target_week[test_mask_aa],
    Region = seqs_age_aware$Region[test_mask_aa],
    Age = seqs_age_aware$Age[test_mask_aa],
    residual_pred = pred_test_cnn,
    stringsAsFactors = FALSE
)
test_pred_df_cnn$AgeNumeric <- safe_map_age(test_pred_df_cnn$Age)

test_residuals_wide_cnn <- test_pred_df_cnn %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(AgeNumeric, Region),
        values_from = residual_pred,
        values_fn = mean,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_test_wide_cnn <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide,
    residuals_wide = test_residuals_wide_cnn,
    week_col = "Week"
)

train_pred_df_cnn <- data.frame(
    Week = seqs_age_aware$target_week[valid_train_mask_aa],
    Region = seqs_age_aware$Region[valid_train_mask_aa],
    Age = seqs_age_aware$Age[valid_train_mask_aa],
    residual_pred = pred_train_cnn,
    stringsAsFactors = FALSE
)
train_pred_df_cnn$AgeNumeric <- safe_map_age(train_pred_df_cnn$Age)

train_residuals_wide_cnn <- train_pred_df_cnn %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(AgeNumeric, Region),
        values_from = residual_pred,
        values_fn = mean,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_train_wide_cnn <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide_train,
    residuals_wide = train_residuals_wide_cnn,
    week_col = "Week"
)

all_results[["CNN-LSTM"]] <- list(
    model_name = "CNN-LSTM",
    adjusted_train_wide = adjusted_train_wide_cnn,
    adjusted_test_wide = adjusted_test_wide_cnn,
    intercept_fix = intercept_fix_per_age
)

################################################################################
# MODEL 2: MortFCNet
################################################################################
start_time_mortfcnet <- Sys.time()

mdl_mortfcnet <- cnn_mort_fc_net(
    seq_length = seq_length_mort,
    input_size = n_climate_features_mort,
    lstm_units = best_mort$lstm_units,
    hidden_size_1 = best_mort$hidden_size_1,
    hidden_size_2 = best_mort$hidden_size_2,
    hidden_size_3 = best_mort$hidden_size_3,
    output_size = 1,
    cell_type = "gru",
    use_cnn = FALSE,
    use_dense = TRUE,
    use_age_embedding = FALSE,
    use_age_embedding_pre_lstm = FALSE,
    use_region_embedding = FALSE,
    dropout_rate = best_mort$dropout_rate
)

opt_mort <- keras::optimizer_adam(
    learning_rate = best_mort$learning_rate,
    clipnorm = 0.5
)
mdl_mortfcnet$compile(
    optimizer = opt_mort,
    loss = "mse",
    metrics = list("mse")
)

history_mortfcnet <- keras::fit(
    mdl_mortfcnet,
    x = X_fit_mort,
    y = y_fit_mort,
    shuffle = TRUE,
    epochs = epochs_final_mort,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_mort,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

end_time_mortfcnet <- Sys.time()
time_mortfcnet <- difftime(
    end_time_mortfcnet,
    start_time_mortfcnet,
    units = "mins"
)
epochs_mortfcnet <- length(history_mortfcnet$metrics$loss)

model_timings[["MortFCNet"]] <- list(
    time_minutes = as.numeric(time_mortfcnet),
    epochs_to_converge = epochs_mortfcnet,
    final_train_loss = tail(history_mortfcnet$metrics$loss, 1)
)

# Reconstruct MortFCNet predictions on the common residual scale.
X_train_all_ic_mort <- X_train_mort
pred_train_mort_z <- as.vector(predict(mdl_mortfcnet, X_train_all_ic_mort))

pred_test_mort_z <- as.vector(predict(mdl_mortfcnet, X_test_mort))

pred_test_mortfcnet <- pred_test_mort_z *
    y_scale_mort_common +
    y_center_mort_common
pred_train_mortfcnet <- pred_train_mort_z *
    y_scale_mort_common +
    y_center_mort_common

# Build wide tables — MortFCNet
test_pred_df_mortfcnet <- data.frame(
    Week = seqs_age_aware_mort$target_week[test_mask_mort],
    Region = seqs_age_aware_mort$Region[test_mask_mort],
    residual_pred = pred_test_mortfcnet,
    stringsAsFactors = FALSE
)

test_residuals_region_wide_mortfcnet <- test_pred_df_mortfcnet %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = Region,
        values_from = residual_pred,
        values_fn = mean,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_test_wide_mortfcnet <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide_mort,
    residuals_wide = broadcast_region_residuals_to_baseline(
        test_residuals_region_wide_mortfcnet,
        baseline_logmxt_wide_mort
    ),
    week_col = "Week"
)

train_pred_df_mortfcnet <- data.frame(
    Week = seqs_age_aware_mort$target_week[valid_train_mask_mort],
    Region = seqs_age_aware_mort$Region[valid_train_mask_mort],
    residual_pred = pred_train_mortfcnet,
    stringsAsFactors = FALSE
)

train_residuals_region_wide_mortfcnet <- train_pred_df_mortfcnet %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = Region,
        values_from = residual_pred,
        values_fn = mean,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_train_wide_mortfcnet <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide_train_mort,
    residuals_wide = broadcast_region_residuals_to_baseline(
        train_residuals_region_wide_mortfcnet,
        baseline_logmxt_wide_train_mort
    ),
    week_col = "Week"
)

all_results[["MortFCNet"]] <- list(
    model_name = "MortFCNet",
    adjusted_train_wide = adjusted_train_wide_mortfcnet,
    adjusted_test_wide = adjusted_test_wide_mortfcnet,
    intercept_fix = 0
)

################################################################################
# MODEL 3: GNN-LSTM
################################################################################
start_time_gnn <- Sys.time()

# Train on the full training window with fixed epochs.
fit_idx_gnn <- seq_len(length(y_train_gnn))

X_fit_gnn <- X_train_gnn
y_fit_gnn <- y_train_gnn
age_idx_fit_gnn <- age_idx_train_gnn
region_idx_fit_gnn <- region_idx_train_gnn

# Derive mid-point numeric age values for GNN samples from age labels.
age_mid_train_gnn <- matrix(
    as.numeric(safe_map_age(age_name_train_gnn)) + 2.5,
    ncol = 1
)
age_mid_test_gnn <- matrix(
    as.numeric(safe_map_age(age_name_test_gnn)) + 2.5,
    ncol = 1
)
# Standardise the GNN numeric age input on the full training set.
age_mid_gnn_center <- mean(age_mid_train_gnn[fit_idx_gnn], na.rm = TRUE)
age_mid_gnn_scale <- sd(age_mid_train_gnn[fit_idx_gnn], na.rm = TRUE)
if (!is.finite(age_mid_gnn_scale) || age_mid_gnn_scale == 0) {
    age_mid_gnn_scale <- 1
}
age_mid_train_gnn <- (age_mid_train_gnn - age_mid_gnn_center) /
    age_mid_gnn_scale
age_mid_test_gnn <- (age_mid_test_gnn - age_mid_gnn_center) / age_mid_gnn_scale
age_mid_fit_gnn <- age_mid_train_gnn

# Weight recent GNN training samples more heavily and upweight peak residual weeks.
sample_weights_fit_gnn <- make_peak_aware_sample_weights(
    target_week_train_gnn,
    observed_residuals = y_train_gnn,
    recency_strength = 0.0,
    peak_strength = 1,
    peak_quantile = 0.90
)

year_numeric_fit_gnn <- year_numeric_train_gnn

mdl_gnn <- gnn_lstm_scalar_output_region_preserving(
    seq_length = seq_length_cnn_gnn,
    n_regions = seqs_gnn$n_regions,
    features_per_region = seqs_gnn$features_per_region,
    adjacency_matrix = adj_matrix_norm_soft,
    n_ages = length(xv),
    gnn_hidden_units = best_gnn$gnn_hidden_units,
    gnn_activation = "relu",
    lstm_units = best_gnn$lstm_units,
    hidden_size_1 = best_gnn$hidden_size_1,
    hidden_size_2 = best_gnn$hidden_size_2,
    hidden_size_3 = best_gnn$hidden_size_3,
    age_embedding_dim = best_gnn$age_embedding_dim,
    region_embedding_dim = best_gnn$region_embedding_dim,
    dropout_rate = best_gnn$dropout_rate,
    use_age_embedding = use_age_embedding_gnn,
    use_region_embedding = use_region_embedding_gnn,
    use_age_numeric = TRUE,
    use_year_numeric = TRUE
)

opt_gnn <- keras::optimizer_adam(
    learning_rate = best_gnn$learning_rate,
    clipnorm = 0.5
)
mdl_gnn$compile(
    optimizer = opt_gnn,
    loss = "mse",
    metrics = list("mse")
)

x_fit_gnn_list <- select_model_inputs(
    mdl_gnn,
    list(
        climate_input = X_fit_gnn,
        age_input = age_idx_fit_gnn,
        region_input = region_idx_fit_gnn,
        age_numeric_input = age_mid_fit_gnn,
        year_numeric_input = year_numeric_fit_gnn
    )
)

# Compile and fit the GNN model.
history_gnn <- keras::fit(
    mdl_gnn,
    x = x_fit_gnn_list,
    y = y_fit_gnn,
    sample_weight = sample_weights_fit_gnn,
    shuffle = TRUE,
    epochs = epochs_final_gnn,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_gnn,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

# Recompile the GNN model after training to keep the object in a valid state.
mdl_gnn$compile(
    optimizer = opt_gnn,
    loss = "mse",
    metrics = list("mse")
)

end_time_gnn <- Sys.time()
time_gnn <- difftime(end_time_gnn, start_time_gnn, units = "mins")
epochs_gnn <- length(history_gnn$metrics$loss)

model_timings[["GNN-LSTM"]] <- list(
    time_minutes = as.numeric(time_gnn),
    epochs_to_converge = epochs_gnn,
    final_train_loss = tail(history_gnn$metrics$loss, 1)
)

# Intercept correction for GNN-LSTM
pred_train_gnn_z <- as.vector(predict(
    mdl_gnn,
    select_model_inputs(
        mdl_gnn,
        list(
            climate_input = X_train_gnn,
            age_input = age_idx_train_gnn,
            region_input = region_idx_train_gnn,
            age_numeric_input = age_mid_train_gnn,
            year_numeric_input = year_numeric_train_gnn
        )
    )
))
age_idx_train_gnn_vec <- age_idx_train_gnn[, 1]

intercept_fix_gnn <- col_bias_original

age_idx_test_gnn_vec <- age_idx_test_gnn[, 1]
col_names_train_gnn <- target_cols[age_idx_train_gnn_vec]
col_names_test_gnn <- target_cols[age_idx_test_gnn_vec]

stopifnot(all(col_names_train_gnn %in% target_cols))
stopifnot(all(col_names_test_gnn %in% target_cols))

pred_test_gnn_z <- as.vector(predict(
    mdl_gnn,
    select_model_inputs(
        mdl_gnn,
        list(
            climate_input = X_test_gnn,
            age_input = age_idx_test_gnn,
            region_input = region_idx_test_gnn,
            age_numeric_input = age_mid_test_gnn,
            year_numeric_input = year_numeric_test_gnn
        )
    )
))

pred_test_gnn <- pred_test_gnn_z *
    y_scale[col_names_test_gnn] +
    y_center[col_names_test_gnn] +
    intercept_fix_gnn[col_names_test_gnn]
pred_train_gnn <- pred_train_gnn_z *
    y_scale[col_names_train_gnn] +
    y_center[col_names_train_gnn] +
    intercept_fix_gnn[col_names_train_gnn]

# Build wide tables — GNN-LSTM
test_pred_df_gnn <- data.frame(
    Week = target_week_test_gnn,
    Age = age_name_test_gnn,
    Region = region_name_test_gnn,
    residual_pred = pred_test_gnn,
    check.names = FALSE
)

test_residuals_wide_gnn <- test_pred_df_gnn %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(Age, Region),
        values_from = residual_pred,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_test_wide_gnn <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide,
    residuals_wide = test_residuals_wide_gnn,
    week_col = "Week"
)

train_pred_df_gnn <- data.frame(
    Week = target_week_train_gnn,
    Age = age_name_train_gnn,
    Region = region_name_train_gnn,
    residual_pred = pred_train_gnn,
    check.names = FALSE
)

train_residuals_wide_gnn <- train_pred_df_gnn %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(Age, Region),
        values_from = residual_pred,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

adjusted_train_wide_gnn <- combine_logmxt_and_residuals_mxt(
    baseline_logmxt_wide = baseline_logmxt_wide_train,
    residuals_wide = train_residuals_wide_gnn,
    week_col = "Week"
)

all_results[["GNN-LSTM"]] <- list(
    model_name = "GNN-LSTM",
    adjusted_train_wide = adjusted_train_wide_gnn,
    adjusted_test_wide = adjusted_test_wide_gnn,
    intercept_fix = intercept_fix_gnn
)

# Sanity check: CNN/GNN intercept fixes should come from precomputed age bias.
cnn_bias_fix_ok <- isTRUE(all.equal(
    intercept_fix_per_age,
    age_group_bias_original,
    check.attributes = FALSE
))
gnn_bias_fix_ok <- isTRUE(all.equal(
    intercept_fix_gnn,
    age_group_bias_original,
    check.attributes = FALSE
))
# ==============================================================================
# PHASE 7: EVALUATION & METRICS COMPUTATION
# ==============================================================================
baseline_logmxt_wide_train$Week <- as.Date(baseline_logmxt_wide_train$Week)
baseline_logmxt_wide$Week <- as.Date(baseline_logmxt_wide$Week)
observed_logmxt_wide$Week <- as.Date(observed_logmxt_wide$Week)
exposure_wide$Week <- as.Date(exposure_wide$Week)

observed_train <- observed_logmxt_wide %>%
    dplyr::filter(
        lubridate::isoyear(Week) %in% train_years
    )

observed_test <- observed_logmxt_wide %>%
    dplyr::filter(
        lubridate::isoyear(Week) %in% test_years
    )

exposure_train <- exposure_wide %>%
    dplyr::filter(
        lubridate::isoyear(Week) %in% train_years
    )

exposure_test <- exposure_wide %>%
    dplyr::filter(
        lubridate::isoyear(Week) %in% test_years
    )

# Baseline metrics
metrics_baseline_train <- compute_metrics(
    baseline_logmxt_wide_train %>%
        dplyr::filter(
            lubridate::isoyear(Week) %in% train_years
        ),
    observed_train
)
metrics_baseline_test <- compute_metrics(
    baseline_logmxt_wide %>%
        dplyr::filter(
            lubridate::isoyear(Week) %in% test_years
        ),
    observed_test
)

# Build metrics table for all models
metrics_rows <- list()

# Baseline
metrics_rows[["BASELINE"]] <- data.frame(
    model = "BASELINE",
    mae_log_test = metrics_baseline_test$mae_log,
    mae_log_train = metrics_baseline_train$mae_log,
    mse_log_test = metrics_baseline_test$mse_log,
    mse_log_train = metrics_baseline_train$mse_log,
    mae_mxt_test = metrics_baseline_test$mae_mxt,
    mae_mxt_train = metrics_baseline_train$mae_mxt,
    mse_mxt_test = metrics_baseline_test$mse_mxt,
    mse_mxt_train = metrics_baseline_train$mse_mxt,
    time_minutes = NA,
    epochs = NA
)

# All models
for (model_name in names(all_results)) {
    res <- all_results[[model_name]]
    metrics_train <- compute_metrics(res$adjusted_train_wide, observed_train)
    metrics_test <- compute_metrics(res$adjusted_test_wide, observed_test)

    timing_info <- model_timings[[model_name]]
    # Guard against missing timing info (NULL fields break data.frame)
    timing_minutes_val <- if (!is.null(timing_info$time_minutes)) {
        timing_info$time_minutes
    } else {
        NA_real_
    }
    timing_epochs_val <- if (!is.null(timing_info$epochs_to_converge)) {
        timing_info$epochs_to_converge
    } else {
        NA_integer_
    }

    metrics_rows[[model_name]] <- data.frame(
        model = model_name,
        mae_log_test = metrics_test$mae_log,
        mae_log_train = metrics_train$mae_log,
        mse_log_test = metrics_test$mse_log,
        mse_log_train = metrics_train$mse_log,
        mae_mxt_test = metrics_test$mae_mxt,
        mae_mxt_train = metrics_train$mae_mxt,
        mse_mxt_test = metrics_test$mse_mxt,
        mse_mxt_train = metrics_train$mse_mxt,
        time_minutes = timing_minutes_val,
        epochs = timing_epochs_val
    )
}

metrics_table <- do.call(rbind, metrics_rows)
numeric_cols <- sapply(metrics_table, is.numeric)
metrics_table[, numeric_cols] <- round(metrics_table[, numeric_cols], 6)

################################################################################
# PER-AGE METRICS
################################################################################
age_levels <- as.character(xv)

models_data <- list()
models_data[["BASELINE"]] <- list(
    train = baseline_logmxt_wide_train,
    test = baseline_logmxt_wide
)

for (model_name in names(all_results)) {
    models_data[[model_name]] <- list(
        train = all_results[[model_name]]$adjusted_train_wide,
        test = all_results[[model_name]]$adjusted_test_wide
    )
}

per_age_rows <- list()
for (age in age_levels) {
    for (mname in names(models_data)) {
        age_cols_train <- get_age_metric_cols(models_data[[mname]]$train, age)
        age_cols_test <- get_age_metric_cols(models_data[[mname]]$test, age)
        pred_train_df <- models_data[[mname]]$train %>%
            dplyr::select(Week, dplyr::all_of(age_cols_train)) %>%
            dplyr::arrange(Week)
        pred_test_df <- models_data[[mname]]$test %>%
            dplyr::select(Week, dplyr::all_of(age_cols_test)) %>%
            dplyr::arrange(Week)
        obs_train_df <- observed_train %>%
            dplyr::select(
                Week,
                dplyr::all_of(get_age_metric_cols(observed_train, age))
            ) %>%
            dplyr::arrange(Week)
        obs_test_df <- observed_test %>%
            dplyr::select(
                Week,
                dplyr::all_of(get_age_metric_cols(observed_test, age))
            ) %>%
            dplyr::arrange(Week)

        metrics_train <- compute_metrics(pred_train_df, obs_train_df)
        metrics_test <- compute_metrics(pred_test_df, obs_test_df)

        per_age_rows[[paste(age, mname, sep = "_")]] <- data.frame(
            age = as.integer(age),
            model = mname,
            mae_log_test = metrics_test$mae_log,
            mae_log_train = metrics_train$mae_log,
            mse_log_test = metrics_test$mse_log,
            mse_log_train = metrics_train$mse_log,
            mae_mxt_test = metrics_test$mae_mxt,
            mae_mxt_train = metrics_train$mae_mxt,
            mse_mxt_test = metrics_test$mse_mxt,
            mse_mxt_train = metrics_train$mse_mxt,
            stringsAsFactors = FALSE
        )
    }
}

per_age_metrics_table <- do.call(rbind, per_age_rows)
per_age_metrics_table <- per_age_metrics_table %>%
    dplyr::select(age, model, mse_mxt_train, mse_mxt_test)

# Compute baseline per-age for comparison
baseline_per_age <- per_age_metrics_table %>%
    dplyr::filter(model == "BASELINE") %>%
    dplyr::select(age, baseline_mse_mxt_test = mse_mxt_test)

per_age_metrics_table <- per_age_metrics_table %>%
    dplyr::left_join(baseline_per_age, by = "age") %>%
    dplyr::mutate(
        diff_vs_baseline_test = mse_mxt_test - baseline_mse_mxt_test,
        pct_improve_vs_baseline = ifelse(
            is.na(baseline_mse_mxt_test) | baseline_mse_mxt_test == 0,
            NA_real_,
            (baseline_mse_mxt_test - mse_mxt_test) / baseline_mse_mxt_test * 100
        )
    )

per_region_rows <- list()
for (reg in regions) {
    for (mname in names(models_data)) {
        if (is.na(mname) || is.null(mname) || nchar(mname) == 0) {
            next
        }

        region_cols_train <- get_region_metric_cols(
            models_data[[mname]]$train,
            reg
        )
        region_cols_test <- get_region_metric_cols(
            models_data[[mname]]$test,
            reg
        )
        pred_train_df <- models_data[[mname]]$train %>%
            dplyr::select(Week, dplyr::all_of(region_cols_train)) %>%
            dplyr::arrange(Week)
        pred_test_df <- models_data[[mname]]$test %>%
            dplyr::select(Week, dplyr::all_of(region_cols_test)) %>%
            dplyr::arrange(Week)
        obs_train_df <- observed_train %>%
            dplyr::select(
                Week,
                dplyr::all_of(get_region_metric_cols(observed_train, reg))
            ) %>%
            dplyr::arrange(Week)
        obs_test_df <- observed_test %>%
            dplyr::select(
                Week,
                dplyr::all_of(get_region_metric_cols(observed_test, reg))
            ) %>%
            dplyr::arrange(Week)

        metrics_train <- compute_metrics(pred_train_df, obs_train_df)
        metrics_test <- compute_metrics(pred_test_df, obs_test_df)

        rowkey <- paste(reg, mname, sep = "_")
        per_region_rows[[rowkey]] <- data.frame(
            region = reg,
            model = mname,
            mse_log_train = metrics_train$mse_log,
            mse_log_test = metrics_test$mse_log,
            mse_mxt_train = metrics_train$mse_mxt,
            mse_mxt_test = metrics_test$mse_mxt,
            stringsAsFactors = FALSE
        )
    }
}

per_region_metrics_table <- do.call(rbind, per_region_rows)

################################################################################
# IMPROVEMENT CHECKS (MortFCNet, CNN-LSTM, GNN-LSTM vs BASELINE)
################################################################################
focus_models <- c("MortFCNet", "CNN-LSTM", "GNN-LSTM")

age_improvement_summary <- per_age_metrics_table %>%
    dplyr::filter(model %in% focus_models) %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
        n_ages = dplyr::n(),
        n_ages_improved = sum(diff_vs_baseline_test < 0, na.rm = TRUE),
        pct_ages_improved = 100 * n_ages_improved / n_ages,
        mean_pct_improve = mean(pct_improve_vs_baseline, na.rm = TRUE),
        .groups = "drop"
    )

baseline_region_test <- per_region_metrics_table %>%
    dplyr::filter(model == "BASELINE") %>%
    dplyr::select(region, baseline_mse_mxt_test = mse_mxt_test)

region_improvement_detail <- per_region_metrics_table %>%
    dplyr::filter(model %in% focus_models) %>%
    dplyr::left_join(baseline_region_test, by = "region") %>%
    dplyr::mutate(
        diff_vs_baseline_test = mse_mxt_test - baseline_mse_mxt_test,
        pct_improve_vs_baseline = ifelse(
            is.na(baseline_mse_mxt_test) | baseline_mse_mxt_test == 0,
            NA_real_,
            (baseline_mse_mxt_test - mse_mxt_test) /
                baseline_mse_mxt_test *
                100
        )
    )

region_improvement_summary <- region_improvement_detail %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
        n_regions = dplyr::n(),
        n_regions_improved = sum(diff_vs_baseline_test < 0, na.rm = TRUE),
        pct_regions_improved = 100 * n_regions_improved / n_regions,
        mean_pct_improve = mean(pct_improve_vs_baseline, na.rm = TRUE),
        .groups = "drop"
    )

improvement_check <- age_improvement_summary %>%
    dplyr::rename(mean_age_pct_improve = mean_pct_improve) %>%
    dplyr::left_join(
        region_improvement_summary %>%
            dplyr::rename(mean_region_pct_improve = mean_pct_improve),
        by = "model"
    ) %>%
    dplyr::mutate(
        improves_age_and_region = (mean_age_pct_improve > 0) &
            (mean_region_pct_improve > 0)
    )


if (any(!improvement_check$improves_age_and_region, na.rm = TRUE)) {
    warning(
        "Some focus models do not improve over BASELINE on both age and ",
        "region mean pct metrics; inspect improvement_check output."
    )
}

################################################################################
# FINAL RESULTS SUMMARY
################################################################################
desired_order <- c("BASELINE", "MortFCNet", "CNN-LSTM", "GNN-LSTM")
available_models <- unique(c(
    "BASELINE",
    names(all_results),
    as.character(metrics_table$model),
    as.character(per_region_metrics_table$model)
))
present_order <- desired_order[desired_order %in% available_models]

# Reorder overall metrics table
metrics_table <- metrics_table[
    match(present_order, metrics_table$model),
    ,
    drop = FALSE
]
# Reorder per-age table (by age then requested model order)
per_age_metrics_table <- per_age_metrics_table %>%
    dplyr::mutate(model = factor(model, levels = present_order)) %>%
    dplyr::arrange(age, model)
# Ensure per-region table uses factor levels for plot/legend order
per_region_metrics_table <- per_region_metrics_table %>%
    dplyr::mutate(model = factor(model, levels = present_order))

invisible(list(
    metrics_table = metrics_table,
    per_age_metrics_table = per_age_metrics_table,
    per_region_metrics_table = per_region_metrics_table
))
