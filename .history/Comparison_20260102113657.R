rm(list = ls(all.names = TRUE))
gc()
set.seed(42L)
options(dplyr.summarise.inform = FALSE) # Suppress summarise info

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
    'keras'
)
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
    install.packages(new.packages)
}
suppressMessages(sapply(packages, require, character.only = TRUE))


library(keras3)
library(ggplot2)

# Set working directory to the script's directory
setwd("C:\1 My Code\Climate-Mortality-Forecasting")

# Input parameters that need to be changed
ISOstart = 2013
ISOend = 2019
start_date <- as.Date("2013-01-01")
end_date <- as.Date("2019-12-31")

# Weekly Mortality Model
t_min <- as.Date('2013-01-01', format = '%Y-%m-%d') # Start of ISO week 1 in 2013
t_max <- as.Date('2024-06-30', format = '%Y-%m-%d') # End of ISO week 26 in 2024

time_frame <- seq(t_min, t_max, by = "days")
df_itime <- data.frame(
    'Date' = time_frame,
    'ISOWeek' = isoweek(time_frame),
    'ISOYear' = isoyear(time_frame)
)

# Load shape file NUTS 2 European regions
shapefile <- read_sf('../../Data/NUTS/NUTS_RG_20M_2021_3035.shp')

# Extract shapefile for NUTS 2 regions in the countries of interest
nuts.spec <- 2
ctry.spec <- 'FR'
shapef <- shapefile[
    shapefile$CNTR_CODE %in% ctry.spec & shapefile$LEVL_CODE == nuts.spec,
]

# Remove overseas areas
overseas <- c(paste0('FRY', 1:5), 'FRM0')
ind.rm <- unlist(sapply(overseas, function(x) {
    which(grepl(x, shapef$NUTS_ID, fixed = TRUE))
}))
shapef <- shapef[-c(ind.rm), ] %>%
    dplyr::arrange(NUTS_ID)
regions <- sort(unique(shapef$NUTS_ID))

# Load weekly deaths and preprocess
d.xtw.all <- get_eurostat(
    id = 'demo_r_mweek3',
    cache = FALSE,
    compress_file = FALSE,
    time_format = 'raw',
    filters = list(geo = c(shapef$NUTS_ID))
)

colnames(d.xtw.all) <- c(
    'Freq',
    'Unit',
    'Sex',
    'Age',
    'Region',
    'Time',
    'Deaths'
)

# Filter and pre-process
d.xtw.c <- d.xtw.all %>%
    dplyr::filter(Region %in% shapef$NUTS_ID & !Age %in% c('UNK', 'TOTAL')) %>%
    mutate(
        ISOYear = as.integer(substr(Time, 1, 4)),
        ISOWeek = as.integer(substr(Time, 7, 9)),
        Age = factor(
            Age,
            levels = c(
                'Y_LT5',
                paste0('Y', seq(5, 85, 5), '-', seq(9, 89, 5)),
                'Y_GE90'
            )
        ),
        .before = 1
    ) %>%
    dplyr::select(-c('Freq', 'Unit', 'Time')) %>%
    dplyr::arrange(Sex, ISOYear, ISOWeek, Age) %>%
    dplyr::filter(
        !grepl(':', Deaths),
        as.integer(Age) >= 14,
        ISOYear >= year(t_min) &
            ((ISOYear < year(t_max)) | (ISOYear == year(t_max) & ISOWeek <= 26))
    )
# Drop levels
d.xtw.c$Age <- droplevels(d.xtw.c$Age)

# Focus on age groups above 65 / combined data on both sexes
d.xtw.c <- d.xtw.c %>%
    dplyr::filter(
        Age %in%
            c('Y65-69', 'Y70-74', 'Y75-79', 'Y80-84', 'Y85-89', 'Y_GE90') &
            Sex == 'T'
    )

# Add date
wdate <- ISOweek2date(sprintf(
    "%d-W%02d-%d",
    d.xtw.c$ISOYear,
    d.xtw.c$ISOWeek,
    1
))
d.xtw.c <- d.xtw.c %>%
    mutate('Date' = wdate, .before = 1)

# Remove Sex column
df.d <- d.xtw.c %>%
    dplyr::select(-c(Sex))

# Annual exposures from Eurostat
P.xt.all <- get_eurostat(
    id = 'demo_r_d2jan',
    cache = FALSE,
    compress_file = FALSE,
    time_format = 'raw',
    filters = list(geo = c(shapef$NUTS_ID))
)

colnames(P.xt.all) <- c('Freq', 'Unit', 'Sex', 'Age', 'Region', 'Time', 'Pop')

# Annual exposures- countries of interest
P.xt.c <- P.xt.all %>%
    dplyr::filter(Region %in% shapef$NUTS_ID, !Age %in% c('UNK', 'TOTAL')) %>%
    mutate(
        ISOYear = as.integer(as.character(Time)),
        Age = factor(
            Age,
            levels = c('Y_LT1', paste0('Y', seq(1, 99, 1)), 'Y_OPEN')
        ),
        .before = 1
    ) %>%
    dplyr::select(-c('Freq', 'Unit', 'Time')) %>%
    arrange(Sex, ISOYear, Age) %>%
    dplyr::filter(
        !grepl(':', Pop),
        as.numeric(Age) >= 66,
        ISOYear <= year(t_max) + 1,
        ISOYear >= year(t_min)
    )

# Drop levels
P.xt.c$Age <- droplevels(P.xt.c$Age)

# Focus on unisex data
P.xt.c <- P.xt.c %>%
    dplyr::filter(Sex == "T") %>%
    dplyr::select(-c(Sex))

# Extrapolation for up to the year 2025 (January 1)
list.FR <- P.xt.c %>%
    split(~ Region + Age)
newdf <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(newdf) <- c('ISOYear', 'Age', 'Region', 'Pop')

for (s in names(list.FR)) {
    # Population counts for region r
    sub <- list.FR[[s]]

    # LM fit
    lmr <- mgcv::gam(Pop ~ s(ISOYear), data = sub)

    # Missing year(s)
    years <- year(t_min):(year(t_max) + 1)
    missing <- years[!years %in% sub$ISOYear]

    # Prediction on missing year(s)
    pred <- predict(lmr, newdata = data.frame(ISOYear = missing))

    # Add to new data frame
    newdf <- rbind(
        newdf,
        data.frame(
            'ISOYear' = missing,
            'Region' = substr(s, 1, 4),
            'Age' = substr(s, 6, 15),
            'Pop' = pred
        )
    )
}

P.xt.c <- rbind(P.xt.c, newdf)

# Create weekly exposures from population count
E.xt.c <- P.xt.c %>%
    dplyr::group_by(Age, Region) %>%
    arrange(ISOYear) %>%
    reframe(
        ISOYear,
        'Expo' = c(
            (Pop[-length(ISOYear)] + Pop[-1]) / 2,
            Pop[length(ISOYear)]
        ) /
            52.1775
    ) %>%
    dplyr::filter(ISOYear <= year(t_max))

# Group per age group
E.xt.c$Age <- plyr::mapvalues(
    E.xt.c$Age,
    from = levels(E.xt.c$Age),
    to = c(
        rep('Y65-69', 5),
        rep('Y70-74', 5),
        rep('Y75-79', 5),
        rep('Y80-84', 5),
        rep('Y85-89', 5),
        rep('Y_GE90', 11)
    )
)
E.xt.c <- E.xt.c %>%
    group_by(ISOYear, Age, Region) %>%
    reframe(Expo = sum(Expo)) %>%
    ungroup()

# Harmonize types for join keys and Expo
E.xt.c <- E.xt.c %>%
    dplyr::mutate(
        ISOYear = as.integer(ISOYear),
        Age = as.character(Age),
        Region = as.character(Region),
        Expo = as.numeric(Expo)
    )

df.d <- df.d %>%
    dplyr::mutate(
        ISOYear = as.integer(ISOYear),
        Age = as.character(Age),
        Region = as.character(Region)
    )

# Add exposure to death counts data and compute rates
df.d <- df.d %>%
    dplyr::left_join(E.xt.c, by = c('ISOYear', 'Age', 'Region')) %>%
    dplyr::mutate(
        Deaths = as.numeric(Deaths),
        mxt = dplyr::if_else(!is.na(Expo) & Expo > 0, Deaths / Expo, NA_real_),
        logmxt = dplyr::if_else(!is.na(mxt), log(pmax(mxt, 1e-12)), NA_real_)
    )

# Fit Lee-Carter Model
source('Functions/wLC-NB.R')

# Filter data for model fitting
DfM <- df.d %>%
    dplyr::filter(
        !is.na(Expo),
        !is.na(Deaths),
        ISOYear >= ISOstart,
        ISOYear <= ISOend,
        ISOWeek <= 52
    ) %>%
    dplyr::arrange(Region, Age, Date)

# Filter data for observed mortality up to 2019
DfM_full <- df.d %>%
    dplyr::filter(
        !is.na(Expo),
        !is.na(Deaths),
        ISOYear >= ISOstart,
        ISOYear <= ISOend,
        ISOWeek <= 52
    ) %>%
    dplyr::arrange(Region, Age, Date)

# Convert age labels to numeric midpoints for plotting
age_map <- c(
    'Y65-69' = 67,
    'Y70-74' = 72,
    'Y75-79' = 77,
    'Y80-84' = 82,
    'Y85-89' = 87,
    'Y_GE90' = 92
)
DfM$age_numeric <- age_map[DfM$Age]
DfM_full$age_numeric <- age_map[DfM_full$Age]

# Map numeric center ages to human-readable age-group labels and a short note
age_group_map <- c(
    "67" = "65-69",
    "72" = "70-74",
    "77" = "75-79",
    "82" = "80-84",
    "87" = "85-89",
    "92" = "90+"
)
age_map_note <- paste0(
    sapply(names(age_group_map), function(k) {
        paste0(k, "=", age_group_map[[k]])
    }),
    collapse = ", "
)

# Create wide dataframe for observed log mxt from 2013-2019
observed_logmxt_wide <- DfM_full %>%
    dplyr::select(Date, age_numeric, Region, logmxt) %>%
    tidyr::pivot_wider(
        id_cols = Date,
        names_from = c(age_numeric, Region),
        values_from = logmxt,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Date) %>%
    dplyr::rename(Week = Date)

# Create matrices: deaths and exposures by (time, age, region)
dtxr <- stats::xtabs(Deaths ~ Date + age_numeric + Region, data = DfM)
etxr <- stats::xtabs(Expo ~ Date + age_numeric + Region, data = DfM)

# Avoid numerical issues with zero deaths
dtxr[dtxr == 0] <- 0.001

# Extract dimensions
xv <- as.integer(colnames(dtxr)) # Age values
isoyv <- lubridate::isoyear(as.Date(rownames(dtxr))) - min(DfM$ISOYear) + 1
isowv <- lubridate::isoweek(as.Date(rownames(dtxr)))
R <- dim(dtxr)[3] # Number of regions

lca <- fit701M.nb(xv, isoyv, isowv, etxr, dtxr, xv * 0, 'ALL')

# fit701M.nb does not always attach dimnames. Add them so downstream
# forecasting code can map ages/regions/weeks reliably.
regions_fit <- dimnames(dtxr)[[3]]
ages_fit <- as.character(xv)
weeks_fit <- sort(unique(isowv))

if (!is.null(lca$beta1)) {
    if (is.null(rownames(lca$beta1)) && length(ages_fit) == nrow(lca$beta1)) {
        rownames(lca$beta1) <- ages_fit
    }
    if (
        is.null(colnames(lca$beta1)) && length(regions_fit) == ncol(lca$beta1)
    ) {
        colnames(lca$beta1) <- regions_fit
    }
}

if (!is.null(lca$kappa2) && is.null(colnames(lca$kappa2))) {
    if (length(regions_fit) == ncol(lca$kappa2)) {
        colnames(lca$kappa2) <- regions_fit
    }
}

if (!is.null(lca$kappa3)) {
    if (
        is.null(colnames(lca$kappa3)) && length(regions_fit) == ncol(lca$kappa3)
    ) {
        colnames(lca$kappa3) <- regions_fit
    }
    if (
        is.null(rownames(lca$kappa3)) && length(weeks_fit) == nrow(lca$kappa3)
    ) {
        rownames(lca$kappa3) <- as.character(weeks_fit)
    }
}

# Extract parameter estimates with confidence intervals
I <- solve(-lca$H.base)
sd1 <- sqrt(diag(I))

# Forecast Kappa2
kappa2_mat <- as.matrix(lca$kappa2)
region_names <- sort(unique(DfM$Region))
if (is.null(colnames(kappa2_mat))) {
    colnames(kappa2_mat) <- region_names
}

years_kappa2 <- ISOstart - 1 + seq_len(nrow(kappa2_mat))
if (length(years_kappa2) != nrow(kappa2_mat)) {
    stop("Could not map kappa2 rows to calendar years")
}

forecast_kappa2_random_walk_1y <- function(
    kappa2_mat,
    year_last,
    n_sim = 100,
    seed = 42L,
    level = 0.95
) {
    # Minimal validation (keep errors readable)
    if (!is.matrix(kappa2_mat) || nrow(kappa2_mat) < 2) {
        stop("kappa2_mat must be a matrix with at least 2 rows")
    }
    if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
        stop("level must be a single number in (0, 1)")
    }

    drift <- apply(kappa2_mat, 2, function(x) mean(diff(x), na.rm = TRUE))
    noise_sd <- apply(kappa2_mat, 2, function(x) {
        stats::sd(diff(x), na.rm = TRUE)
    })
    noise_sd[is.na(noise_sd) | noise_sd == 0] <- 0

    last_val <- kappa2_mat[nrow(kappa2_mat), ]
    mu <- as.numeric(last_val + drift)

    if (is.null(n_sim) || !is.numeric(n_sim) || length(n_sim) != 1) {
        stop("n_sim must be a single numeric value")
    }

    if (n_sim <= 0) {
        alpha <- (1 - level) / 2
        lower <- stats::qnorm(alpha, mean = mu, sd = noise_sd)
        upper <- stats::qnorm(1 - alpha, mean = mu, sd = noise_sd)
    } else {
        set.seed(seed)
        sims <- sapply(seq_along(mu), function(j) {
            mu[j] + stats::rnorm(n_sim, mean = 0, sd = noise_sd[j])
        })
        alpha <- (1 - level) / 2
        lower <- apply(sims, 2, stats::quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(
            sims,
            2,
            stats::quantile,
            probs = 1 - alpha,
            na.rm = TRUE
        )
    }

    data.frame(
        year = year_last + 1,
        region = colnames(kappa2_mat),
        kappa2_mean = mu,
        kappa2_lower = as.numeric(lower),
        kappa2_upper = as.numeric(upper),
        stringsAsFactors = FALSE
    )
}

kappa2_forecast_1y <- forecast_kappa2_random_walk_1y(
    kappa2_mat = kappa2_mat,
    year_last = max(years_kappa2),
    n_sim = 100,
    seed = 42L,
    level = 0.95
)

kappa2_mean_1y <- stats::setNames(
    kappa2_forecast_1y$kappa2_mean,
    kappa2_forecast_1y$region
)

# Forecast Weekly Mortality Rates
# Model formula in JensCode/wLC-NB.R:
#   mhat[week, age, region] = exp(beta1[age,region] + beta2[age]*kappa2[year,region] + beta3[age]*kappa3[week,region])

forecast_weekly_mortality_rate_1y <- function(
    lca,
    kappa2_next_mean,
    year_forecast,
    iso_weeks = 1:52
) {
    if (!is.list(lca)) {
        stop("lca must be a list")
    }
    if (is.null(names(kappa2_next_mean))) {
        stop(
            "kappa2_next_mean must be a named numeric vector (names = regions)"
        )
    }

    beta1 <- as.matrix(lca$beta1)
    beta2 <- as.numeric(lca$beta2)
    beta3 <- as.numeric(lca$beta3)
    kappa3 <- as.matrix(lca$kappa3)

    regions <- colnames(beta1)
    if (is.null(regions) || length(regions) == 0) {
        # Fall back to kappa2_next_mean names if beta1 has no colnames
        regions <- names(kappa2_next_mean)
    }
    if (is.null(regions) || length(regions) == 0) {
        stop(
            "Could not determine region names (beta1 has no colnames and kappa2_next_mean has no names)"
        )
    }

    if (is.null(rownames(beta1))) {
        rownames(beta1) <- as.character(seq_len(nrow(beta1)))
    }
    ages_numeric <- rownames(beta1)

    if (is.null(colnames(kappa3)) && ncol(kappa3) == length(regions)) {
        colnames(kappa3) <- regions
    }

    if (!all(regions %in% names(kappa2_next_mean))) {
        missing <- setdiff(regions, names(kappa2_next_mean))
        stop(
            "kappa2_next_mean is missing regions: ",
            paste(missing, collapse = ", ")
        )
    }
    kappa2_next <- as.numeric(kappa2_next_mean[regions])

    # If kappa3 has rownames like "1","2",..., match by name; else fall back to positional.
    iso_weeks <- as.integer(iso_weeks)
    if (!is.null(rownames(kappa3))) {
        weeks_avail <- as.integer(rownames(kappa3))
        iso_weeks <- iso_weeks[iso_weeks %in% weeks_avail]
        if (length(iso_weeks) == 0) {
            stop(
                "None of the requested iso_weeks are present in rownames(lca$kappa3)"
            )
        }
    } else {
        iso_weeks <- iso_weeks[iso_weeks >= 1 & iso_weeks <= nrow(kappa3)]
        if (length(iso_weeks) == 0) {
            stop("iso_weeks is empty after filtering to available kappa3 rows")
        }
    }

    # Week start (Monday) dates for the forecast ISO year.
    week_dates <- ISOweek::ISOweek2date(sprintf(
        "%d-W%02d-1",
        year_forecast,
        iso_weeks
    ))

    n_weeks <- length(iso_weeks)
    n_ages <- nrow(beta1)
    n_regions <- length(regions)

    eta_forecast <- array(
        NA_real_,
        dim = c(n_weeks, n_ages, n_regions),
        dimnames = list(
            Week = as.character(week_dates),
            Age = ages_numeric,
            Region = regions
        )
    )

    mhat_forecast <- array(
        NA_real_,
        dim = c(n_weeks, n_ages, n_regions),
        dimnames = list(
            Week = as.character(week_dates),
            Age = ages_numeric,
            Region = regions
        )
    )

    for (w_idx in seq_len(n_weeks)) {
        w <- iso_weeks[w_idx]
        k3_w <- if (!is.null(rownames(kappa3))) {
            as.numeric(kappa3[as.character(w), regions])
        } else {
            as.numeric(kappa3[w, regions])
        }

        for (j in seq_len(n_ages)) {
            eta <- as.numeric(beta1[j, regions]) +
                beta2[j] * kappa2_next +
                beta3[j] * k3_w
            eta_forecast[w_idx, j, ] <- eta
            mhat_forecast[w_idx, j, ] <- exp(eta)
        }
    }

    forecast_long <- expand.grid(
        Week = week_dates,
        AgeNumeric = ages_numeric,
        Region = regions,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    forecast_long$ISOYear <- year_forecast
    forecast_long$ISOWeek <- rep(iso_weeks, times = n_ages * n_regions)
    forecast_long$logmxt_forecast <- as.vector(eta_forecast)
    forecast_long$mxt_forecast <- as.vector(mhat_forecast)

    # Optional: attach age-group labels when the script's age_map is available
    if (exists("age_map", inherits = TRUE)) {
        age_label_by_num <- stats::setNames(
            names(age_map),
            as.character(age_map)
        )
        forecast_long$AgeGroup <- unname(age_label_by_num[as.character(
            forecast_long$AgeNumeric
        )])
    }

    list(
        mxt_array = mhat_forecast,
        logmxt_array = eta_forecast,
        forecast_long = forecast_long
    )
}

combine_logmxt_and_residuals <- function(
    baseline_logmxt_wide,
    residuals_wide,
    week_col = "Week"
) {
    if (
        !is.data.frame(baseline_logmxt_wide) || !is.data.frame(residuals_wide)
    ) {
        stop("baseline_logmxt_wide and residuals_wide must be data.frames")
    }
    if (
        !(week_col %in% names(baseline_logmxt_wide)) ||
            !(week_col %in% names(residuals_wide))
    ) {
        stop("Both inputs must contain week_col = '", week_col, "'")
    }

    b <- baseline_logmxt_wide
    r <- residuals_wide
    b[[week_col]] <- as.Date(b[[week_col]])
    r[[week_col]] <- as.Date(r[[week_col]])

    common_weeks <- sort(intersect(b[[week_col]], r[[week_col]]))
    if (length(common_weeks) > 0) {
        b <- b[b[[week_col]] %in% common_weeks, , drop = FALSE]
        r <- r[r[[week_col]] %in% common_weeks, , drop = FALSE]
        b <- b[order(b[[week_col]]), , drop = FALSE]
        r <- r[order(r[[week_col]]), , drop = FALSE]
        week_out <- b[[week_col]]
    } else {
        # Assume rows correspond in order; truncate to common length.
        n <- min(nrow(b), nrow(r))
        if (n == 0) {
            stop("No rows available to combine")
        }
        b <- b[seq_len(n), , drop = FALSE]
        r <- r[seq_len(n), , drop = FALSE]
        week_out <- r[[week_col]]
    }

    b_cols <- setdiff(names(b), week_col)
    r_cols <- setdiff(names(r), week_col)

    if (identical(b_cols, r_cols)) {
        common_cols <- b_cols
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    } else if (setequal(b_cols, r_cols)) {
        common_cols <- b_cols
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    } else if (length(b_cols) == length(r_cols)) {
        # Assume same column order / dimensions; ignore names.
        common_cols <- b_cols
        b_mat <- as.matrix(b[, b_cols, drop = FALSE])
        r_mat <- as.matrix(r[, r_cols, drop = FALSE])
        colnames(r_mat) <- common_cols
    } else {
        common_cols <- intersect(b_cols, r_cols)
        if (length(common_cols) == 0) {
            stop(
                "No overlapping forecast columns and column counts differ; cannot combine"
            )
        }
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    }

    out_mat <- exp(b_mat + r_mat)

    out <- as.data.frame(out_mat)
    out[[week_col]] <- week_out
    out <- out[, c(week_col, common_cols), drop = FALSE]
    out
}

year_last_fit <- max(years_kappa2)
year_forecast <- year_last_fit + 1

mortality_rate_forecast_1y <- forecast_weekly_mortality_rate_1y(
    lca = lca,
    kappa2_next_mean = kappa2_mean_1y,
    year_forecast = year_forecast,
    iso_weeks = 1:52
)

# Tidy long table: one row per (week, age, region)
mortality_rate_forecast_1y_long <- mortality_rate_forecast_1y$forecast_long

# Convenience wide table for log(mxt) using the same columns as mxt wide.
mortality_log_rate_forecast_1y_wide <- mortality_rate_forecast_1y_long %>%
    dplyr::select(Week, AgeNumeric, Region, logmxt_forecast) %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(AgeNumeric, Region),
        values_from = logmxt_forecast,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

# If AgeGroup exists, also provide AgeGroup_Region columns (matches residuals_matrix naming).
if ("AgeGroup" %in% names(mortality_rate_forecast_1y_long)) {
    mortality_log_rate_forecast_1y_wide_agegroup <- mortality_rate_forecast_1y_long %>%
        dplyr::select(Week, AgeGroup, Region, logmxt_forecast) %>%
        tidyr::pivot_wider(
            id_cols = Week,
            names_from = c(AgeGroup, Region),
            values_from = logmxt_forecast,
            names_sep = "_"
        ) %>%
        dplyr::arrange(Week)
}

# Wide table aligned to residuals_matrix column naming: "<AgeNumeric>_<Region>"
mortality_rate_forecast_1y_wide <- mortality_rate_forecast_1y_long %>%
    dplyr::select(Week, AgeNumeric, Region, mxt_forecast) %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(AgeNumeric, Region),
        values_from = mxt_forecast,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

# Combine baseline forecast with test residuals for adjusted forecast
# Function to combine log mxt and residuals (returns adjusted log mxt)
combine_logmxt_and_residuals_log <- function(
    baseline_logmxt_wide,
    residuals_wide,
    week_col = "Week"
) {
    if (
        !is.data.frame(baseline_logmxt_wide) || !is.data.frame(residuals_wide)
    ) {
        stop("baseline_logmxt_wide and residuals_wide must be data.frames")
    }
    if (
        !(week_col %in% names(baseline_logmxt_wide)) ||
            !(week_col %in% names(residuals_wide))
    ) {
        stop("Both inputs must contain week_col = '", week_col, "'")
    }

    b <- baseline_logmxt_wide
    r <- residuals_wide
    b[[week_col]] <- as.Date(b[[week_col]])
    r[[week_col]] <- as.Date(r[[week_col]])

    common_weeks <- sort(intersect(b[[week_col]], r[[week_col]]))
    if (length(common_weeks) > 0) {
        b <- b[b[[week_col]] %in% common_weeks, , drop = FALSE]
        r <- r[r[[week_col]] %in% common_weeks, , drop = FALSE]
        b <- b[order(b[[week_col]]), , drop = FALSE]
        r <- r[order(r[[week_col]]), , drop = FALSE]
        week_out <- b[[week_col]]
    } else {
        # Assume rows correspond in order; truncate to common length.
        n <- min(nrow(b), nrow(r))
        if (n == 0) {
            stop("No rows available to combine")
        }
        b <- b[seq_len(n), , drop = FALSE]
        r <- r[seq_len(n), , drop = FALSE]
        week_out <- r[[week_col]]
    }

    b_cols <- setdiff(names(b), week_col)
    r_cols <- setdiff(names(r), week_col)

    if (identical(b_cols, r_cols)) {
        common_cols <- b_cols
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    } else if (setequal(b_cols, r_cols)) {
        common_cols <- b_cols
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    } else if (length(b_cols) == length(r_cols)) {
        # Assume same column order / dimensions; ignore names.
        common_cols <- b_cols
        b_mat <- as.matrix(b[, b_cols, drop = FALSE])
        r_mat <- as.matrix(r[, r_cols, drop = FALSE])
        colnames(r_mat) <- common_cols
    } else {
        common_cols <- intersect(b_cols, r_cols)
        if (length(common_cols) == 0) {
            stop(
                "No overlapping forecast columns and column counts differ; cannot combine"
            )
        }
        b_mat <- as.matrix(b[, common_cols, drop = FALSE])
        r_mat <- as.matrix(r[, common_cols, drop = FALSE])
    }

    out_mat <- b_mat + r_mat # Add for adjusted log mxt

    out <- as.data.frame(out_mat)
    out[[week_col]] <- week_out
    out <- out[, c(week_col, common_cols), drop = FALSE]
    out
}

##############################################################
#################### CALCULATE RESIDUALS #####################
##############################################################
log_obs <- log(
    stats::xtabs(Deaths ~ Date + age_numeric + Region, data = DfM) /
        stats::xtabs(Expo ~ Date + age_numeric + Region, data = DfM)
) # vector: <obs, age group, region>
log_est = log(lca$mhat) # vector: <obs, age group, region>
additive_res <- log_obs - log_est

residuals_aligned <- additive_res

# Reshape additive_res
residuals_matrix <- matrix(
    as.vector(additive_res),
    nrow = dim(additive_res)[1],
    ncol = dim(additive_res)[2] * dim(additive_res)[3]
)

# set column names
age_groups <- dimnames(additive_res)[[2]]
region_names <- dimnames(additive_res)[[3]]
col_labels <- as.vector(outer(age_groups, region_names, paste, sep = "_"))
colnames(residuals_matrix) <- col_labels

# Set rownames to weekly date
date_vec <- dimnames(additive_res)[[1]]
region_vec <- region_names
rownames(residuals_matrix) <- dimnames(additive_res)[[1]]

# Create baseline log mxt wide from Lee-Carter fitted values
baseline_long <- expand.grid(
    Week = as.Date(dimnames(additive_res)[[1]]),
    AgeNumeric = dimnames(additive_res)[[2]],
    Region = dimnames(additive_res)[[3]],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
baseline_long$logmxt_baseline <- as.vector(log_est)

baseline_logmxt_wide_2013_2018 <- baseline_long %>%
    tidyr::pivot_wider(
        id_cols = Week,
        names_from = c(AgeNumeric, Region),
        values_from = logmxt_baseline,
        names_sep = "_"
    ) %>%
    dplyr::arrange(Week)

# Combine fitted 2013-2018 with forecasted 2019
baseline_logmxt_wide <- dplyr::bind_rows(
    baseline_logmxt_wide_2013_2018,
    mortality_log_rate_forecast_1y_wide
)

##############################################################
###################### Climate Data ##########################
##############################################################
source('Functions/K.ReadCAMS.R')
source('Functions/K.ReadEOBS.R')

# CAMS Data
agg <- DatasetAggregator$new(dir = "c:/1 My Code/FullDataset/Data/CAMS")
agg$aggregate_files()
CAMS_avg <- agg$select_avg_cols()

# EOBS Data
EOBS <- read_eobs_files(combine = TRUE)

# Merging the Datasets
# Ensure Date columns are Date class
if ("Date" %in% names(CAMS_avg) && !inherits(CAMS_avg$Date, "Date")) {
    CAMS_avg$Date <- as.Date(CAMS_avg$Date)
}
if ("Date" %in% names(EOBS) && !inherits(EOBS$Date, "Date")) {
    EOBS$Date <- as.Date(EOBS$Date)
}

# Find a common region column name present in both datasets (adjust candidates if needed)
get_region_col <- function(
    a,
    b,
    candidates = c(
        "region",
        "Region",
        "NUTS_ID",
        "nuts_id",
        "id",
        "ID",
        "cell",
        "CellID"
    )
) {
    ca <- intersect(candidates, names(a))
    cb <- intersect(candidates, names(b))
    common <- intersect(ca, cb)
    if (length(common) == 0) {
        stop("No common region column found. Set region_col manually.")
    }
    common[1]
}
region_col <- get_region_col(CAMS_avg, EOBS)

# Compute intersection of dates.
common_dates <- intersect(unique(CAMS_avg$Date), unique(EOBS$Date))
common_dates <- common_dates[
    common_dates >= start_date & common_dates <= end_date
]
if (length(common_dates) == 0) {
    stop("No common dates found")
}

CAMS_f <- CAMS_avg[CAMS_avg$Date %in% common_dates, , drop = FALSE]
EOBS_f <- EOBS[EOBS$Date %in% common_dates, , drop = FALSE]

# Merge / join by Date and region (inner join => intersection)
climate_df <- merge(CAMS_f, EOBS_f, by = c("Date", region_col), all = FALSE)

# rename pollutant columns to short names
nm <- colnames(climate_df)
nm[grepl("no2", nm, ignore.case = TRUE)] <- "no2"
nm[grepl("\\bo3\\b|(^|_)o3(_|\\.|$)|o3\\.", nm, ignore.case = TRUE)] <- "o3"
nm[grepl("pm10", nm, ignore.case = TRUE)] <- "pm10"
nm[grepl("pm2p5|pm2\\.5|pm25", nm, ignore.case = TRUE)] <- "pm2p5"
colnames(climate_df) <- make.unique(nm)

# Drop FRM0 region (Corsica)
climate_df <- climate_df[climate_df[[region_col]] != "FRM0", ]

# Standardize region column name for downstream use.
names(climate_df)[names(climate_df) == region_col] <- "Region"

# Aggregate climate to weekly per Region (NUTS_ID) and add Country_Code.

climate_df <- climate_df %>%
    dplyr::mutate(
        Week = floor_date(Date, unit = "week", week_start = 1),
        Country_Code = substr(Region, 1, 2)
    )

measurement_cols <- setdiff(
    names(climate_df),
    c("Date", "Week", "Region", "Country_Code")
)

climate_weekly_region <- climate_df %>%
    dplyr::group_by(Week, Region, Country_Code) %>%
    dplyr::summarise(
        dplyr::across(tidyselect::all_of(measurement_cols), mean, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::arrange(Region, Week) %>%
    dplyr::filter(Week >= as.Date("2013-01-01"))

# Standardize climate variables (keep one-hot vars unscaled).
climate_weekly_region_z <- climate_weekly_region
climate_weekly_region_z[, measurement_cols] <- scale(
    climate_weekly_region[, measurement_cols, drop = FALSE]
)

##############################################################
################### CNN-GRU/LSTM Forecast Residuals ###################
##############################################################
cnn_mort_fc_net <- function(
    seq_length,
    input_size,
    cnn_filters = c(128, 64, 32),
    hidden_size_1 = 128,
    hidden_size_2 = 64,
    hidden_size_3 = 32,
    output_size = 1,
    cell_type = c("lstm", "gru"),
    use_cnn = TRUE,
    use_dense = TRUE
) {
    cell_type <- match.arg(cell_type)
    inputs <- layer_input(shape = c(seq_length, input_size))

    x <- inputs
    if (use_cnn) {
        # CNN layers for dimension reduction
        for (f in cnn_filters) {
            x <- x %>%
                layer_conv_1d(
                    filters = f,
                    kernel_size = 3,
                    padding = 'same',
                    activation = 'relu'
                ) %>%
                layer_batch_normalization() %>%
                layer_dropout(rate = 0.1)
        }
        lstm_units <- if (use_dense) tail(cnn_filters, 1) else output_size
    } else {
        # Skip CNN layers
        lstm_units <- if (use_dense) input_size else output_size
    }

    # Recurrent layer: choose LSTM or GRU by `cell_type`
    if (cell_type == "gru") {
        x <- x %>%
            layer_gru(units = lstm_units, return_sequences = FALSE)
    } else {
        x <- x %>%
            layer_lstm(units = lstm_units, return_sequences = FALSE)
    }

    if (use_dense) {
        # FC1 -> LayerNorm -> LeakyReLU -> Dropout
        x <- x %>%
            layer_dense(units = hidden_size_1, use_bias = TRUE) %>%
            layer_layer_normalization() %>%
            layer_activation_leaky_relu(alpha = 0.3) %>%
            layer_dropout(rate = 0.1)

        # FC2 -> LayerNorm -> LeakyReLU -> Dropout
        x <- x %>%
            layer_dense(units = hidden_size_2, use_bias = TRUE) %>%
            layer_layer_normalization() %>%
            layer_activation_leaky_relu(alpha = 0.3) %>%
            layer_dropout(rate = 0.1)

        # FC3 -> LayerNorm -> LeakyReLU -> Dropout
        x <- x %>%
            layer_dense(units = hidden_size_3, use_bias = TRUE) %>%
            layer_layer_normalization() %>%
            layer_activation_leaky_relu(alpha = 0.3) %>%
            layer_dropout(rate = 0.1)

        # FC4 (output)
        outputs <- x %>%
            layer_dense(units = output_size, use_bias = TRUE)
    } else {
        # No dense layers, output is the recurrent output
        outputs <- x
    }

    model <- keras_model(inputs = inputs, outputs = outputs)
    return(model)
}

make_sequences <- function(x_mat, y_mat, seq_length) {
    n <- nrow(x_mat)
    n_samples <- n - seq_length
    n_features <- ncol(x_mat)
    n_targets <- ncol(y_mat)

    x_arr <- array(0, dim = c(n_samples, seq_length, n_features))
    y_out <- matrix(0, nrow = n_samples, ncol = n_targets)
    target_row_index <- integer(n_samples)

    for (i in seq_len(n_samples)) {
        start <- i
        end <- i + seq_length - 1
        target <- end + 1
        x_arr[i, , ] <- x_mat[start:end, , drop = FALSE]
        y_out[i, ] <- y_mat[target, , drop = TRUE]
        target_row_index[i] <- target
    }

    list(X = x_arr, y = y_out, target_row_index = target_row_index)
}

make_sequences_by_region <- function(df_key, X_mat, y_mat, seq_length) {
    regions <- sort(unique(df_key$Region))
    n_features <- ncol(X_mat)
    n_targets <- ncol(y_mat)

    n_samples_by_region <- sapply(regions, function(r) {
        n_r <- sum(df_key$Region == r)
        max(n_r - seq_length, 0)
    })
    n_total <- sum(n_samples_by_region)

    X_all <- array(0, dim = c(n_total, seq_length, n_features))
    y_all <- matrix(0, nrow = n_total, ncol = n_targets)
    target_week <- as.Date(rep(NA, n_total))
    target_region <- rep(NA_character_, n_total)

    pos <- 1
    for (r in regions) {
        idx <- which(df_key$Region == r)
        idx <- idx[order(df_key$Week[idx])]

        x_sub <- X_mat[idx, , drop = FALSE]
        y_sub <- y_mat[idx, , drop = FALSE]
        wk_sub <- df_key$Week[idx]

        seqs <- make_sequences(x_sub, y_sub, seq_length = seq_length)
        n_r <- dim(seqs$X)[1]

        if (n_r > 0) {
            rng <- pos:(pos + n_r - 1)
            X_all[rng, , ] <- seqs$X
            y_all[rng, ] <- seqs$y
            target_week[rng] <- wk_sub[seqs$target_row_index]
            target_region[rng] <- r
            pos <- pos + n_r
        }
    }

    list(
        X = X_all,
        y = y_all,
        target_week = target_week,
        Region = target_region
    )
}

# Prepare targets (y): residuals per (Week, Region), columns = age groups.
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

df_model <- climate_weekly_region_z %>%
    dplyr::inner_join(residual_wide_region, by = c("Week", "Region")) %>%
    dplyr::arrange(Region, Week)

target_cols <- setdiff(names(residual_wide_region), c("Week", "Region"))

# Prepare features (X): climate variables + one-hot NUTS_ID and Country_Code.
X_num <- as.matrix(df_model[, measurement_cols, drop = FALSE])

one_hot_safe <- function(x, prefix) {
    x <- as.factor(x)
    lv <- levels(x)

    if (length(lv) <= 1) {
        m <- matrix(1, nrow = length(x), ncol = 1)
        colnames(m) <- paste0(prefix, lv[1])
        return(m)
    }

    m <- stats::model.matrix(~ 0 + x)
    colnames(m) <- sub("^x", prefix, colnames(m))
    m
}

X_region <- one_hot_safe(df_model$Region, prefix = "Region_")
X_country <- one_hot_safe(df_model$Country_Code, prefix = "Country_Code_")
X_mat <- cbind(X_num, X_region, X_country)
y_mat <- as.matrix(df_model[, target_cols, drop = FALSE])

# IMPORTANT PARAMETERS
forecast_horizon <- 52
seq_length <- 12

# Standardize targets for training stability; keep params to invert later.
y_center <- colMeans(y_mat, na.rm = TRUE)
y_scale <- apply(y_mat, 2, stats::sd, na.rm = TRUE)
y_scale[y_scale == 0 | is.na(y_scale)] <- 1
y_mat_z <- scale(y_mat, center = y_center, scale = y_scale)

unscale_targets <- function(scaled, center, scale) {
    sweep(scaled, 2, scale, "*") + center
}

seqs <- make_sequences_by_region(
    df_key = df_model[, c("Week", "Region"), drop = FALSE],
    X_mat = X_mat,
    y_mat = y_mat_z,
    seq_length = seq_length
)

# Hold out last 52 target weeks.
weeks_all <- sort(unique(seqs$target_week))
train_end_week <- length(weeks_all) - forecast_horizon
cutoff_week <- weeks_all[train_end_week]

train_mask <- seqs$target_week <= cutoff_week
test_mask <- seqs$target_week > cutoff_week

X_train <- seqs$X[train_mask, , , drop = FALSE]
y_train <- seqs$y[train_mask, , drop = FALSE]
X_test <- seqs$X[test_mask, , , drop = FALSE]
y_test <- seqs$y[test_mask, , drop = FALSE]

n_features <- dim(X_train)[3]
n_targets <- ncol(y_train)

# Train and evaluate both CNN-GRU and CNN-LSTM variants, returning adjusted forecasts.
run_variant <- function(
    cell_type,
    epochs = 80,
    use_cnn = TRUE,
    use_dense = TRUE
) {
    mdl <- cnn_mort_fc_net(
        seq_length = seq_length,
        input_size = n_features,
        output_size = n_targets,
        cell_type = cell_type,
        use_cnn = use_cnn,
        use_dense = use_dense
    )

    base_lr <- 3e-5
    lr_sched <- callback_learning_rate_scheduler(function(epoch, lr) {
        base_lr * (0.5^(epoch %/% 20))
    })

    mdl$compile(
        optimizer = keras$optimizers$Adam(learning_rate = base_lr),
        loss = "mse",
        metrics = list("mse")
    )

    keras::fit(
        mdl,
        x = X_train,
        y = y_train,
        epochs = epochs,
        batch_size = 8,
        callbacks = list(lr_sched),
        verbose = 2
    )

    # Test predictions
    pred_test_z <- predict(mdl, X_test)
    pred_test <- unscale_targets(
        pred_test_z,
        center = y_center,
        scale = y_scale
    )

    test_pred_df <- data.frame(
        Week = seqs$target_week[test_mask],
        Region = seqs$Region[test_mask],
        pred_test,
        check.names = FALSE
    )
    colnames(test_pred_df)[3:ncol(test_pred_df)] <- paste0(
        "AgeGroup_",
        colnames(test_pred_df)[3:ncol(test_pred_df)]
    )

    test_residuals_wide <- test_pred_df %>%
        tidyr::pivot_longer(
            cols = tidyselect::all_of(colnames(test_pred_df)[
                3:ncol(test_pred_df)
            ]),
            names_to = "AgeNumeric",
            values_to = "residual_pred"
        ) %>%
        dplyr::mutate(AgeNumeric = sub("AgeGroup_", "", AgeNumeric)) %>%
        tidyr::pivot_wider(
            id_cols = Week,
            names_from = c(AgeNumeric, Region),
            values_from = residual_pred,
            names_sep = "_"
        ) %>%
        dplyr::arrange(Week)

    # Forecast residuals for 2019 using climate forecast sequences
    pred_2019_z <- predict(mdl, seqs_forecast$X)
    pred_2019 <- unscale_targets(
        pred_2019_z,
        center = y_center,
        scale = y_scale
    )

    forecast_pred_df <- data.frame(
        Week = seqs_forecast$target_week,
        Region = seqs_forecast$Region,
        pred_2019,
        check.names = FALSE
    )
    colnames(forecast_pred_df)[3:ncol(forecast_pred_df)] <- paste0(
        "AgeGroup_",
        colnames(forecast_pred_df)[3:ncol(forecast_pred_df)]
    )

    forecast_residuals_wide <- forecast_pred_df %>%
        tidyr::pivot_longer(
            cols = tidyselect::all_of(colnames(forecast_pred_df)[
                3:ncol(forecast_pred_df)
            ]),
            names_to = "AgeNumeric",
            values_to = "residual_pred"
        ) %>%
        dplyr::mutate(AgeNumeric = sub("AgeGroup_", "", AgeNumeric)) %>%
        tidyr::pivot_wider(
            id_cols = Week,
            names_from = c(AgeNumeric, Region),
            values_from = residual_pred,
            names_sep = "_"
        ) %>%
        dplyr::arrange(Week)

    adjusted_test_wide <- combine_logmxt_and_residuals_log(
        baseline_logmxt_wide = mortality_log_rate_forecast_1y_wide,
        residuals_wide = forecast_residuals_wide,
        week_col = "Week"
    )

    # Train predictions -> adjusted train
    pred_train_z <- predict(mdl, X_train)
    pred_train <- unscale_targets(
        pred_train_z,
        center = y_center,
        scale = y_scale
    )
    train_pred_df <- data.frame(
        Week = seqs$target_week[train_mask],
        Region = seqs$Region[train_mask],
        pred_train,
        check.names = FALSE
    )
    colnames(train_pred_df)[3:ncol(train_pred_df)] <- paste0(
        "AgeGroup_",
        colnames(train_pred_df)[3:ncol(train_pred_df)]
    )
    train_residuals_wide <- train_pred_df %>%
        tidyr::pivot_longer(
            cols = tidyselect::all_of(colnames(train_pred_df)[
                3:ncol(train_pred_df)
            ]),
            names_to = "AgeNumeric",
            values_to = "residual_pred"
        ) %>%
        dplyr::mutate(AgeNumeric = sub("AgeGroup_", "", AgeNumeric)) %>%
        tidyr::pivot_wider(
            id_cols = Week,
            names_from = c(AgeNumeric, Region),
            values_from = residual_pred,
            names_sep = "_"
        ) %>%
        dplyr::arrange(Week)

    adjusted_train_wide <- combine_logmxt_and_residuals_log(
        baseline_logmxt_wide = baseline_logmxt_wide_2013_2018,
        residuals_wide = train_residuals_wide,
        week_col = "Week"
    )

    list(
        cell_type = cell_type,
        use_cnn = use_cnn,
        use_dense = use_dense,
        adjusted_test_wide = adjusted_test_wide,
        adjusted_train_wide = adjusted_train_wide
    )
}

# Precompute forecast sequences used by variants (re-using earlier logic)
forecast_weeks <- ISOweek::ISOweek2date(sprintf("2019-W%02d-1", 1:52))
start_forecast <- min(forecast_weeks) - seq_length * 7
climate_forecast <- climate_weekly_region_z %>%
    dplyr::filter(Week >= start_forecast & Week <= max(forecast_weeks)) %>%
    dplyr::arrange(Region, Week)
X_num_forecast <- as.matrix(climate_forecast[, measurement_cols, drop = FALSE])
X_region_forecast <- one_hot_safe(climate_forecast$Region, prefix = "Region_")
X_country_forecast <- one_hot_safe(
    climate_forecast$Country_Code,
    prefix = "Country_Code_"
)
X_mat_forecast <- cbind(X_num_forecast, X_region_forecast, X_country_forecast)
dummy_y_forecast <- matrix(0, nrow = nrow(climate_forecast), ncol = n_targets)
seqs_forecast <- make_sequences_by_region(
    df_key = climate_forecast[, c("Week", "Region"), drop = FALSE],
    X_mat = X_mat_forecast,
    y_mat = dummy_y_forecast,
    seq_length = seq_length
)

# Run all variants (CNN-GRU, CNN-LSTM, LSTM-only, GRU-only)
# use fewer epochs for quick comparison; adjust `epochs_compare` as needed
variant_configs <- expand.grid(
    cell_type = c("gru", "lstm"),
    use_cnn = c(TRUE, FALSE),
    use_dense = c(TRUE, FALSE),
    stringsAsFactors = FALSE
)
epochs_compare <- 2
results <- lapply(1:nrow(variant_configs), function(i) {
    config <- variant_configs[i, ]
    run_variant(
        cell_type = config$cell_type,
        epochs = epochs_compare,
        use_cnn = config$use_cnn,
        use_dense = config$use_dense
    )
})

# Function to compute mean logmxt per week (inserted before variant summarization)
compute_mean_logmxt <- function(df, type_name) {
    df %>%
        dplyr::mutate(
            mean_logmxt = rowMeans(dplyr::select(., -Week), na.rm = TRUE)
        ) %>%
        dplyr::select(Week, mean_logmxt) %>%
        dplyr::mutate(type = type_name)
}

# Precompute observed and baseline mean series used in the combined plot
observed_mean <- compute_mean_logmxt(
    observed_logmxt_wide,
    "Observed (2013-2019)"
)
baseline_2019_mean <- compute_mean_logmxt(
    baseline_logmxt_wide,
    "Baseline Forecast (Lee-Carter, 2013-2019)"
)

# Compute mean logmxt per variant for plotting (use display names)
display_cell_name <- function(cell) {
    toupper(cell)
}

variant_means <- lapply(results, function(res) {
    model_type <- if (res$use_cnn) {
        "CNN-"
    } else if (res$use_dense) {
        "Encoder-"
    } else {
        ""
    }
    train_mean <- compute_mean_logmxt(
        res$adjusted_train_wide,
        paste0(
            "Adjusted Train (",
            model_type,
            display_cell_name(res$cell_type),
            ")"
        )
    )
    test_mean <- compute_mean_logmxt(
        res$adjusted_test_wide,
        paste0(
            "Adjusted Test (",
            model_type,
            display_cell_name(res$cell_type),
            ")"
        )
    )
    list(train = train_mean, test = test_mean)
})

# Bind all plot data together
plot_data <- dplyr::bind_rows(
    observed_mean,
    baseline_2019_mean,
    lapply(variant_means, function(vm) vm$train),
    lapply(variant_means, function(vm) vm$test)
) %>%
    dplyr::distinct(Week, type, .keep_all = TRUE) %>%
    dplyr::mutate(
        linetype = dplyr::case_when(
            grepl("Adjusted Test", type) ~ "dashed",
            TRUE ~ "solid"
        )
    )

## Build combined plot for observed, baseline, CNN-GRU, CNN-LSTM, GRU-only and LSTM-only variants
windows()
plot_data_combined <- plot_data %>%
    dplyr::mutate(type = as.character(type))

ggplot(
    plot_data_combined,
    aes(x = Week, y = mean_logmxt, color = type, linetype = linetype)
) +
    geom_line(aes(alpha = type), size = 1.1, na.rm = TRUE) +
    scale_color_manual(
        values = c(
            "Observed (2013-2019)" = "black",
            "Baseline Forecast (Lee-Carter, 2013-2019)" = "blue",
            "Adjusted Train (CNN-GRU)" = "red",
            "Adjusted Test (CNN-GRU)" = "red",
            "Adjusted Train (CNN-LSTM)" = "forestgreen",
            "Adjusted Test (CNN-LSTM)" = "forestgreen",
            "Adjusted Train (Encoder-GRU)" = "orange",
            "Adjusted Test (Encoder-GRU)" = "orange",
            "Adjusted Train (Encoder-LSTM)" = "purple",
            "Adjusted Test (Encoder-LSTM)" = "purple"
        )
    ) +
    scale_alpha_manual(
        values = c(
            "Observed (2013-2019)" = 1,
            "Baseline Forecast (Lee-Carter, 2013-2019)" = 0.6,
            "Adjusted Train (CNN-GRU)" = 0.6,
            "Adjusted Test (CNN-GRU)" = 0.6,
            "Adjusted Train (CNN-LSTM)" = 0.6,
            "Adjusted Test (CNN-LSTM)" = 0.6,
            "Adjusted Train (Encoder-GRU)" = 0.6,
            "Adjusted Test (Encoder-GRU)" = 0.6,
            "Adjusted Train (Encoder-LSTM)" = 0.6,
            "Adjusted Test (Encoder-LSTM)" = 0.6
        )
    ) +
    scale_linetype_identity() +
    scale_x_date(limits = as.Date(c("2013-01-01", "2019-12-31"))) +
    labs(
        title = "Mean Log Mortality: Observed vs Baseline vs CNN-GRU/CNN-LSTM/Encoder-GRU/Encoder-LSTM/GRU/LSTM Adjusted",
        subtitle = age_map_note,
        x = "Week",
        y = "Mean Log Mxt",
        color = "Series"
    ) +
    guides(alpha = "none", linetype = "none") +
    theme_minimal()

# ggsave("combined_logmortality_gru_lstm.png", width = 12, height = 8, dpi = 300)

# Plot for a specific age group and region (random selection)
available_ages <- c("67", "72", "77", "82", "87", "92")
available_regions <- regions
specific_age <- sample(available_ages, 1)
specific_region <- sample(available_regions, 1)
specific_col <- paste0(specific_age, "_", specific_region)

# Helper to safely extract a single age-region column and label it
compute_specific_logmxt <- function(df, type_name, col_name) {
    if (!("Week" %in% names(df))) {
        stop("Input df must contain 'Week' column")
    }
    if (!(col_name %in% names(df))) {
        # return NA series if the column does not exist
        tibble::tibble(Week = df$Week, logmxt = NA_real_, type = type_name)
    } else {
        tibble::tibble(
            Week = df$Week,
            logmxt = df[[col_name]],
            type = type_name
        )
    }
}

# Observed and baseline
observed_specific <- compute_specific_logmxt(
    observed_logmxt_wide,
    "Observed (2013-2019)",
    specific_col
)
baseline_specific <- compute_specific_logmxt(
    baseline_logmxt_wide,
    "Baseline Forecast (Lee-Carter, 2013-2019)",
    specific_col
)

# Find results for GRU and LSTM (results list contains cell_type)
get_result_by_type <- function(
    results_list,
    cell,
    use_cnn = TRUE,
    use_dense = TRUE
) {
    for (r in results_list) {
        if (
            !is.null(r$cell_type) &&
                tolower(r$cell_type) == tolower(cell) &&
                r$use_cnn == use_cnn &&
                r$use_dense == use_dense
        ) {
            return(r)
        }
    }
    NULL
}

res_gru <- get_result_by_type(results, "gru", TRUE, TRUE)
res_lstm <- get_result_by_type(results, "lstm", TRUE, TRUE)
res_gru_encoder <- get_result_by_type(results, "gru", FALSE, TRUE)
res_lstm_encoder <- get_result_by_type(results, "lstm", FALSE, TRUE)

# Extract adjusted train/test for each variant (safe-if-missing)
get_specific_from_adjusted <- function(adj_df, label, col_name) {
    if (is.null(adj_df)) {
        return(tibble::tibble(
            Week = as.Date(NA),
            logmxt = NA_real_,
            type = label
        ))
    }
    compute_specific_logmxt(adj_df, label, col_name)
}

gru_train_specific <- if (!is.null(res_gru)) {
    get_specific_from_adjusted(
        res_gru$adjusted_train_wide,
        paste0(
            "Adjusted Train (",
            if (res_gru$use_cnn) "CNN-" else "",
            toupper(res_gru$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}
gru_test_specific <- if (!is.null(res_gru)) {
    get_specific_from_adjusted(
        res_gru$adjusted_test_wide,
        paste0(
            "Adjusted Test (",
            if (res_gru$use_cnn) "CNN-" else "",
            toupper(res_gru$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}

lstm_train_specific <- if (!is.null(res_lstm)) {
    get_specific_from_adjusted(
        res_lstm$adjusted_train_wide,
        paste0(
            "Adjusted Train (",
            if (res_lstm$use_cnn) "CNN-" else "",
            toupper(res_lstm$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}
lstm_test_specific <- if (!is.null(res_lstm)) {
    get_specific_from_adjusted(
        res_lstm$adjusted_test_wide,
        paste0(
            "Adjusted Test (",
            if (res_lstm$use_cnn) "CNN-" else "",
            toupper(res_lstm$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}

gru_encoder_train_specific <- if (!is.null(res_gru_encoder)) {
    get_specific_from_adjusted(
        res_gru_encoder$adjusted_train_wide,
        paste0(
            "Adjusted Train (Encoder-",
            toupper(res_gru_encoder$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}
gru_encoder_test_specific <- if (!is.null(res_gru_encoder)) {
    get_specific_from_adjusted(
        res_gru_encoder$adjusted_test_wide,
        paste0(
            "Adjusted Test (Encoder-",
            toupper(res_gru_encoder$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}

lstm_encoder_train_specific <- if (!is.null(res_lstm_encoder)) {
    get_specific_from_adjusted(
        res_lstm_encoder$adjusted_train_wide,
        paste0(
            "Adjusted Train (Encoder-",
            toupper(res_lstm_encoder$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}
lstm_encoder_test_specific <- if (!is.null(res_lstm_encoder)) {
    get_specific_from_adjusted(
        res_lstm_encoder$adjusted_test_wide,
        paste0(
            "Adjusted Test (Encoder-",
            toupper(res_lstm_encoder$cell_type),
            ")"
        ),
        specific_col
    )
} else {
    NULL
}

# Combine into plot data and clean
plot_data_specific <- dplyr::bind_rows(
    observed_specific,
    baseline_specific,
    gru_train_specific,
    gru_test_specific,
    lstm_train_specific,
    lstm_test_specific,
    gru_encoder_train_specific,
    gru_encoder_test_specific,
    lstm_encoder_train_specific,
    lstm_encoder_test_specific
) %>%
    dplyr::distinct(Week, type, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(Week)) %>%
    dplyr::filter(!is.na(logmxt) & is.finite(logmxt)) %>%
    dplyr::mutate(
        linetype = dplyr::case_when(
            grepl("Adjusted Test", type) ~ "dashed",
            TRUE ~ "solid"
        )
    )

# Create the plot
plot_data_specific <- plot_data_specific %>%
    dplyr::mutate(type = as.character(type))
ggplot(
    plot_data_specific,
    aes(x = Week, y = logmxt, color = type, linetype = linetype)
) +
    geom_line(aes(alpha = type), size = 1.2, na.rm = TRUE) +
    scale_color_manual(
        values = c(
            "Observed (2013-2019)" = "black",
            "Baseline Forecast (Lee-Carter, 2013-2019)" = "blue",
            "Adjusted Train (CNN-GRU)" = "red",
            "Adjusted Test (CNN-GRU)" = "red",
            "Adjusted Train (CNN-LSTM)" = "forestgreen",
            "Adjusted Test (CNN-LSTM)" = "forestgreen",
            "Adjusted Train (GRU)" = "orange",
            "Adjusted Test (GRU)" = "orange",
            "Adjusted Train (LSTM)" = "purple",
            "Adjusted Test (LSTM)" = "purple"
        )
    ) +
    scale_alpha_manual(
        values = c(
            "Observed (2013-2019)" = 1,
            "Baseline Forecast (Lee-Carter, 2013-2019)" = 0.6,
            "Adjusted Train (CNN-GRU)" = 0.6,
            "Adjusted Test (CNN-GRU)" = 0.6,
            "Adjusted Train (CNN-LSTM)" = 0.6,
            "Adjusted Test (CNN-LSTM)" = 0.6,
            "Adjusted Train (GRU)" = 0.6,
            "Adjusted Test (GRU)" = 0.6,
            "Adjusted Train (LSTM)" = 0.6,
            "Adjusted Test (LSTM)" = 0.6
        )
    ) +
    scale_linetype_identity() +
    scale_x_date(limits = as.Date(c("2013-01-01", "2019-12-31"))) +
    labs(
        title = paste(
            "Log Mortality Rates for Age",
            specific_age,
            "in",
            specific_region
        ),
        subtitle = age_map_note,
        x = "Week",
        y = "Log Mxt",
        color = "Series"
    ) +
    guides(alpha = "none", linetype = "none") +
    theme_minimal()

# Save the specific plot as PNG
# ggsave("specific_logmortality_gru_lstm.png", width = 12, height = 8, dpi = 300)

# Is this train and test?
# Compute MAE between predicted and observed log mortality
compute_metrics <- function(pred_wide, obs_wide) {
    if (is.null(pred_wide) || is.null(obs_wide)) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            rmse_log = NA_real_,
            rmse_mxt = NA_real_,
            mape_log = NA_real_,
            mape_mxt = NA_real_
        ))
    }
    pred_wide$Week <- as.Date(pred_wide$Week)
    obs_wide$Week <- as.Date(obs_wide$Week)

    common_weeks <- sort(intersect(pred_wide$Week, obs_wide$Week))
    if (length(common_weeks) == 0) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            rmse_log = NA_real_,
            rmse_mxt = NA_real_,
            mape_log = NA_real_,
            mape_mxt = NA_real_
        ))
    }

    obs_sub <- obs_wide[obs_wide$Week %in% common_weeks, , drop = FALSE]
    pred_sub <- pred_wide[pred_wide$Week %in% common_weeks, , drop = FALSE]

    obs_sub <- obs_sub[order(obs_sub$Week), , drop = FALSE]
    pred_sub <- pred_sub[order(pred_sub$Week), , drop = FALSE]

    cols_obs <- setdiff(names(obs_sub), "Week")
    cols_pred <- setdiff(names(pred_sub), "Week")
    common_cols <- intersect(cols_obs, cols_pred)
    if (length(common_cols) == 0) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            rmse_log = NA_real_,
            rmse_mxt = NA_real_,
            mape_log = NA_real_,
            mape_mxt = NA_real_
        ))
    }

    obs_mat <- as.matrix(obs_sub[, common_cols, drop = FALSE])
    pred_mat <- as.matrix(pred_sub[, common_cols, drop = FALSE])

    # Ensure same dims
    if (!all(dim(obs_mat) == dim(pred_mat))) {
        nrow_min <- min(nrow(obs_mat), nrow(pred_mat))
        obs_mat <- obs_mat[seq_len(nrow_min), , drop = FALSE]
        pred_mat <- pred_mat[seq_len(nrow_min), , drop = FALSE]
    }

    # Compute errors on log scale
    diff_log <- pred_mat - obs_mat
    abs_diff_log <- abs(diff_log)
    mae_log <- mean(abs_diff_log[is.finite(abs_diff_log)], na.rm = TRUE)
    rmse_log <- sqrt(mean(diff_log^2[is.finite(diff_log)], na.rm = TRUE))
    mape_log <- mean(
        abs_diff_log /
            abs(obs_mat)[
                is.finite(abs_diff_log) & is.finite(obs_mat) & abs(obs_mat) > 0
            ],
        na.rm = TRUE
    ) *
        100

    # Compute errors on original scale (mortality rates)
    obs_mxt <- exp(obs_mat)
    pred_mxt <- exp(pred_mat)
    diff_mxt <- pred_mxt - obs_mxt
    abs_diff_mxt <- abs(diff_mxt)
    mae_mxt <- mean(abs_diff_mxt[is.finite(abs_diff_mxt)], na.rm = TRUE)
    rmse_mxt <- sqrt(mean(diff_mxt^2[is.finite(diff_mxt)], na.rm = TRUE))
    mape_mxt <- mean(
        abs_diff_mxt /
            abs(obs_mxt)[
                is.finite(abs_diff_mxt) & is.finite(obs_mxt) & abs(obs_mxt) > 0
            ],
        na.rm = TRUE
    ) *
        100

    list(
        mae_log = mae_log,
        mae_mxt = mae_mxt,
        rmse_log = rmse_log,
        rmse_mxt = rmse_mxt,
        mape_log = mape_log,
        mape_mxt = mape_mxt
    )
}

# Calculate for both variants and save results
metrics_rows <- list()
for (res in results) {
    model_name <- if (res$use_cnn) {
        paste0("CNN-", toupper(res$cell_type))
    } else if (res$use_dense) {
        paste0("Encoder-", toupper(res$cell_type))
    } else {
        toupper(res$cell_type)
    }
    metrics_test <- compute_metrics(
        res$adjusted_test_wide,
        observed_logmxt_wide
    )
    metrics_train <- compute_metrics(
        res$adjusted_train_wide,
        observed_logmxt_wide
    )
    msg <- sprintf(
        "%s - MAE (log): test=%.6g, train=%.6g; MAE (mxt): test=%.6g, train=%.6g; RMSE (log): test=%.6g, train=%.6g; RMSE (mxt): test=%.6g, train=%.6g; MAPE (log): test=%.4g%%, train=%.4g%%; MAPE (mxt): test=%.4g%%, train=%.4g%%",
        model_name,
        metrics_test$mae_log,
        metrics_train$mae_log,
        metrics_test$mae_mxt,
        metrics_train$mae_mxt,
        metrics_test$rmse_log,
        metrics_train$rmse_log,
        metrics_test$rmse_mxt,
        metrics_train$rmse_mxt,
        metrics_test$mape_log,
        metrics_train$mape_log,
        metrics_test$mape_mxt,
        metrics_train$mape_mxt
    )
    cat(msg, "\n")
    metrics_rows[[model_name]] <- data.frame(
        model = model_name,
        mae_log_test = metrics_test$mae_log,
        mae_log_train = metrics_train$mae_log,
        mae_mxt_test = metrics_test$mae_mxt,
        mae_mxt_train = metrics_train$mae_mxt,
        rmse_log_test = metrics_test$rmse_log,
        rmse_log_train = metrics_train$rmse_log,
        rmse_mxt_test = metrics_test$rmse_mxt,
        rmse_mxt_train = metrics_train$rmse_mxt,
        mape_log_test = metrics_test$mape_log,
        mape_log_train = metrics_train$mape_log,
        mape_mxt_test = metrics_test$mape_mxt,
        mape_mxt_train = metrics_train$mape_mxt
    )
}

# Calculate metrics for baseline (Lee-Carter)
metrics_baseline_test <- compute_metrics(
    baseline_logmxt_wide,
    observed_logmxt_wide
)
metrics_baseline_train <- compute_metrics(
    baseline_logmxt_wide_2013_2018,
    observed_logmxt_wide
)
msg_baseline <- sprintf(
    "BASELINE - MAE (log): test=%.6g, train=%.6g; MAE (mxt): test=%.6g, train=%.6g; RMSE (log): test=%.6g, train=%.6g; RMSE (mxt): test=%.6g, train=%.6g; MAPE (log): test=%.4g%%, train=%.4g%%; MAPE (mxt): test=%.4g%%, train=%.4g%%",
    metrics_baseline_test$mae_log,
    metrics_baseline_train$mae_log,
    metrics_baseline_test$mae_mxt,
    metrics_baseline_train$mae_mxt,
    metrics_baseline_test$rmse_log,
    metrics_baseline_train$rmse_log,
    metrics_baseline_test$rmse_mxt,
    metrics_baseline_train$rmse_mxt,
    metrics_baseline_test$mape_log,
    metrics_baseline_train$mape_log,
    metrics_baseline_test$mape_mxt,
    metrics_baseline_train$mape_mxt
)
cat(msg_baseline, "\n")
metrics_rows[["BASELINE"]] <- data.frame(
    model = "BASELINE",
    mae_log_test = metrics_baseline_test$mae_log,
    mae_log_train = metrics_baseline_train$mae_log,
    mae_mxt_test = metrics_baseline_test$mae_mxt,
    mae_mxt_train = metrics_baseline_train$mae_mxt,
    rmse_log_test = metrics_baseline_test$rmse_log,
    rmse_log_train = metrics_baseline_train$rmse_log,
    rmse_mxt_test = metrics_baseline_test$rmse_mxt,
    rmse_mxt_train = metrics_baseline_train$rmse_mxt,
    mape_log_test = metrics_baseline_test$mape_log,
    mape_log_train = metrics_baseline_train$mape_log,
    mape_mxt_test = metrics_baseline_test$mape_mxt,
    mape_mxt_train = metrics_baseline_train$mape_mxt
)

metrics_table <- do.call(rbind, metrics_rows)
# Round numeric columns to appropriate precision
numeric_cols <- sapply(metrics_table, is.numeric)
metrics_table[, numeric_cols] <- round(metrics_table[, numeric_cols], 4)
# Round MAPE columns to 2 decimal places for percentage display
mape_cols <- grep("mape", names(metrics_table))
if (length(mape_cols) > 0) {
    metrics_table[, mape_cols] <- round(metrics_table[, mape_cols], 2)
}
metrics_table

# Compute MAE per age group
ages <- c("67", "72", "77", "82", "87", "92")

compute_metrics_per_age <- function(pred_wide, obs_wide, age) {
    pred_cols <- grep(paste0("^", age, "_"), names(pred_wide), value = TRUE)
    obs_cols <- grep(paste0("^", age, "_"), names(obs_wide), value = TRUE)
    common_cols <- intersect(pred_cols, obs_cols)

    if (length(common_cols) == 0) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            rmse_log = NA_real_,
            rmse_mxt = NA_real_,
            mape_log = NA_real_,
            mape_mxt = NA_real_
        ))
    }

    pred_sub <- pred_wide[, c("Week", common_cols), drop = FALSE]
    obs_sub <- obs_wide[, c("Week", common_cols), drop = FALSE]

    compute_metrics(pred_sub, obs_sub)
}

# Metrics per age for each model
metrics_per_age_rows <- list()
for (age in c("67", "72", "77", "82", "87", "92")) {
    for (res in results) {
        model_name <- if (res$use_cnn) {
            paste0("CNN-", toupper(res$cell_type))
        } else if (res$use_dense) {
            paste0("Encoder-", toupper(res$cell_type))
        } else {
            toupper(res$cell_type)
        }
        metrics_test <- compute_metrics_per_age(
            res$adjusted_test_wide,
            observed_logmxt_wide,
            age
        )
        metrics_train <- compute_metrics_per_age(
            res$adjusted_train_wide,
            observed_logmxt_wide,
            age
        )
        metrics_per_age_rows[[paste(age, model_name, sep = "_")]] <- data.frame(
            age = age,
            model = model_name,
            mae_log_test = metrics_test$mae_log,
            mae_log_train = metrics_train$mae_log,
            mae_mxt_test = metrics_test$mae_mxt,
            mae_mxt_train = metrics_train$mae_mxt,
            rmse_log_test = metrics_test$rmse_log,
            rmse_log_train = metrics_train$rmse_log,
            rmse_mxt_test = metrics_test$rmse_mxt,
            rmse_mxt_train = metrics_train$rmse_mxt,
            mape_log_test = metrics_test$mape_log,
            mape_log_train = metrics_train$mape_log,
            mape_mxt_test = metrics_test$mape_mxt,
            mape_mxt_train = metrics_train$mape_mxt
        )
    }

    # Baseline
    metrics_baseline_test <- compute_metrics_per_age(
        baseline_logmxt_wide,
        observed_logmxt_wide,
        age
    )
    metrics_baseline_train <- compute_metrics_per_age(
        baseline_logmxt_wide_2013_2018,
        observed_logmxt_wide,
        age
    )
    metrics_per_age_rows[[paste(age, "BASELINE", sep = "_")]] <- data.frame(
        age = age,
        model = "BASELINE",
        mae_log_test = metrics_baseline_test$mae_log,
        mae_log_train = metrics_baseline_train$mae_log,
        mae_mxt_test = metrics_baseline_test$mae_mxt,
        mae_mxt_train = metrics_baseline_train$mae_mxt,
        rmse_log_test = metrics_baseline_test$rmse_log,
        rmse_log_train = metrics_baseline_train$rmse_log,
        rmse_mxt_test = metrics_baseline_test$rmse_mxt,
        rmse_mxt_train = metrics_baseline_train$rmse_mxt,
        mape_log_test = metrics_baseline_test$mape_log,
        mape_log_train = metrics_baseline_train$mape_log,
        mape_mxt_test = metrics_baseline_test$mape_mxt,
        mape_mxt_train = metrics_baseline_train$mape_mxt
    )
}

metrics_per_age_table <- do.call(rbind, metrics_per_age_rows)
# Round numeric columns
numeric_cols <- sapply(metrics_per_age_table, is.numeric)
metrics_per_age_table[, numeric_cols] <- round(
    metrics_per_age_table[, numeric_cols],
    4
)
# Round MAPE columns to 2 decimal places
mape_cols <- grep("mape", names(metrics_per_age_table))
if (length(mape_cols) > 0) {
    metrics_per_age_table[, mape_cols] <- round(
        metrics_per_age_table[, mape_cols],
        2
    )
}
metrics_per_age_table

# Plot mean log mortality per age group
compute_mean_logmxt_per_age <- function(df, type_name) {
    df_long <- df %>%
        tidyr::pivot_longer(
            cols = -Week,
            names_to = "age_region",
            values_to = "logmxt"
        ) %>%
        dplyr::mutate(
            age = sub("_.*", "", age_region)
        ) %>%
        dplyr::group_by(Week, age) %>%
        dplyr::summarise(
            mean_logmxt = mean(logmxt, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(type = type_name)
    return(df_long)
}

# Compute per age for observed, baseline, GRU, LSTM
observed_per_age <- compute_mean_logmxt_per_age(
    observed_logmxt_wide,
    "Observed (2013-2019)"
)
baseline_per_age <- compute_mean_logmxt_per_age(
    baseline_logmxt_wide,
    "Baseline Forecast (Lee-Carter, 2013-2019)"
)

variant_per_age <- lapply(results, function(res) {
    model_type <- if (res$use_cnn) {
        "CNN-"
    } else if (res$use_dense) {
        "Encoder-"
    } else {
        ""
    }
    train_per_age <- compute_mean_logmxt_per_age(
        res$adjusted_train_wide,
        paste0("Adjusted Train (", model_type, toupper(res$cell_type), ")")
    )
    test_per_age <- compute_mean_logmxt_per_age(
        res$adjusted_test_wide,
        paste0("Adjusted Test (", model_type, toupper(res$cell_type), ")")
    )
    list(train = train_per_age, test = test_per_age)
})

# Bind all per age data
plot_data_per_age <- dplyr::bind_rows(
    observed_per_age,
    baseline_per_age,
    lapply(variant_per_age, function(vm) vm$train),
    lapply(variant_per_age, function(vm) vm$test)
) %>%
    dplyr::distinct(Week, age, type, .keep_all = TRUE) %>%
    dplyr::mutate(
        linetype = dplyr::case_when(
            grepl("Adjusted Test", type) ~ "dashed",
            TRUE ~ "solid"
        )
    ) %>%
    dplyr::mutate(
        age = ifelse(!is.na(age_group_map[age]), age_group_map[age], age)
    )

# Plot faceted by age
ggplot(
    plot_data_per_age,
    aes(x = Week, y = mean_logmxt, color = type, linetype = linetype)
) +
    geom_line(aes(alpha = type), size = 0.8, na.rm = TRUE) +
    facet_wrap(~age, scales = "free_y", ncol = 2) +
    scale_color_manual(
        values = c(
            "Observed (2013-2019)" = "black",
            "Baseline Forecast (Lee-Carter, 2013-2019)" = "blue",
            "Adjusted Train (CNN-GRU)" = "red",
            "Adjusted Test (CNN-GRU)" = "red",
            "Adjusted Train (CNN-LSTM)" = "forestgreen",
            "Adjusted Test (CNN-LSTM)" = "forestgreen",
            "Adjusted Train (GRU)" = "orange",
            "Adjusted Test (GRU)" = "orange",
            "Adjusted Train (LSTM)" = "purple",
            "Adjusted Test (LSTM)" = "purple"
        )
    ) +
    scale_alpha_manual(
        values = c(
            "Observed (2013-2019)" = 1,
            "Baseline Forecast (Lee-Carter, 2013-2019)" = 0.6,
            "Adjusted Train (CNN-GRU)" = 0.6,
            "Adjusted Test (CNN-GRU)" = 0.6,
            "Adjusted Train (CNN-LSTM)" = 0.6,
            "Adjusted Test (CNN-LSTM)" = 0.6,
            "Adjusted Train (GRU)" = 0.6,
            "Adjusted Test (GRU)" = 0.6,
            "Adjusted Train (LSTM)" = 0.6,
            "Adjusted Test (LSTM)" = 0.6
        )
    ) +
    scale_linetype_identity() +
    scale_x_date(limits = as.Date(c("2013-01-01", "2019-12-31"))) +
    labs(
        title = "Mean Log Mortality per Age Group: Observed vs Baseline vs CNN-GRU/CNN-LSTM/Encoder-GRU/Encoder-LSTM/GRU/LSTM Adjusted",
        subtitle = age_map_note,
        x = "Week",
        y = "Mean Log Mxt",
        color = "Series"
    ) +
    guides(alpha = "none", linetype = "none") +
    theme_minimal()

# ggsave("logmortality_per_age_gru_lstm.png", width = 16, height = 12, dpi = 300)
