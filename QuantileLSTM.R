# ==============================================================================
# PHASE 8: UNCERTAINTY DECOMPOSITION PIPELINE
# ==============================================================================
# This section keeps only the uncertainty components requested by the user:
# 1) LC parameter bootstrap from simulated deaths
# 2) kappa_2 forecast bootstraps from RW + AR(1)
# 3) quantile-LSTM split-normal residual sampling using 0.05/0.50/0.95
# 4) empirical 5% / 95% intervals for log-mu on the combined sample paths

# Override with environment variable `BOOT_B`, e.g. Sys.setenv(BOOT_B = "1000").
B_boot <- as.integer(Sys.getenv("BOOT_B", "10000"))

tau_levels <- c(0.05, 0.50, 0.95)
n_quantiles <- length(tau_levels)
pi_lo_idx <- 1L
pi_med_idx <- 2L
pi_hi_idx <- 3L

epochs_q_mort <- 70L
epochs_q_cnn <- 70L
epochs_q_gnn <- 200L
patience_q <- 15L

make_pinball_loss <- function(tau_vec) {
    force(tau_vec)
    tau_arr <- array(tau_vec, dim = c(1L, length(tau_vec)))

    function(y_true, y_pred) {
        tau_t <- keras3::op_cast(
            keras3::op_convert_to_tensor(tau_arr),
            dtype = "float32"
        )
        y_t <- keras3::op_reshape(
            keras3::op_cast(y_true, dtype = "float32"),
            newshape = c(-1L, 1L)
        )
        y_t <- keras3::op_broadcast_to(y_t, keras3::op_shape(y_pred))
        err <- y_t - y_pred
        keras3::op_mean(keras3::op_maximum(tau_t * err, (tau_t - 1.0) * err))
    }
}

pinball_loss <- make_pinball_loss(tau_levels)

sort_quantiles <- function(mat_n_by_q) {
    t(apply(mat_n_by_q, 1, sort))
}

unscale_quantile_matrix <- function(
    q_mat_z,
    col_names,
    y_scale,
    y_center,
    intercept_fix
) {
    scale_vec <- y_scale[col_names]
    center_vec <- y_center[col_names]
    intercept_vec <- intercept_fix[col_names]
    sweep(q_mat_z, 1, scale_vec, "*") +
        matrix(
            center_vec + intercept_vec,
            nrow = nrow(q_mat_z),
            ncol = ncol(q_mat_z)
        )
}

build_uq_long <- function(
    q_mat,
    target_weeks,
    regions_vec,
    ages_vec,
    baseline_logmxt_wide,
    tau_vec,
    model_label
) {
    colnames(q_mat) <- paste0("q", tau_vec * 100)

    base_df <- data.frame(
        Week = target_weeks,
        Region = regions_vec,
        Age = ages_vec,
        stringsAsFactors = FALSE
    )
    base_df <- cbind(base_df, as.data.frame(q_mat))

    long_df <- tidyr::pivot_longer(
        base_df,
        cols = starts_with("q"),
        names_to = "quantile_label",
        values_to = "residual_pred"
    ) %>%
        dplyr::mutate(
            tau = as.numeric(sub("^q", "", quantile_label)) / 100,
            AgeNumeric = safe_map_age(Age),
            model = model_label
        )

    baseline_long <- baseline_logmxt_wide %>%
        tidyr::pivot_longer(
            cols = -Week,
            names_to = c("age_str", "Region"),
            names_sep = "_(?=[^_]+$)",
            values_to = "logmxt_baseline"
        ) %>%
        dplyr::mutate(
            Week = as.Date(Week),
            AgeNumeric = age_str
        ) %>%
        dplyr::select(Week, Region, AgeNumeric, logmxt_baseline)

    long_df %>%
        dplyr::mutate(Week = as.Date(Week)) %>%
        dplyr::left_join(
            baseline_long,
            by = c("Week", "Region", "AgeNumeric")
        ) %>%
        dplyr::mutate(
            mxt_baseline = exp(logmxt_baseline),
            mxt_adjusted_q = mxt_baseline + residual_pred,
            logmxt_adjusted_q = log(pmax(mxt_adjusted_q, 1e-12))
        )
}

compute_test_logmu_array <- function(
    lca_obj,
    kappa2_full_b,
    test_dates,
    forecast_years_b,
    n_ages,
    n_regions
) {
    n_t <- length(test_dates)
    out <- array(NA_real_, dim = c(n_t, n_ages, n_regions))
    for (t in seq_len(n_t)) {
        year_t <- lubridate::isoyear(test_dates[t])
        week_t <- lubridate::isoweek(test_dates[t])
        yi <- resolve_year_index(year_t, forecast_years_b, "boot test LC year")
        for (a in seq_len(n_ages)) {
            out[t, a, ] <- lca_obj$beta1[a, ] +
                lca_obj$beta2[a] * kappa2_full_b[yi, ] +
                lca_obj$beta3[a] * lca_obj$kappa3[week_t, ]
        }
    }
    out
}

start_q_cnn <- Sys.time()

mdl_q_cnn <- cnn_lstm_age_aware_quantile(
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
    use_age_embedding = use_age_embedding_cnn,
    use_region_embedding = use_region_embedding_cnn,
    use_age_numeric = TRUE,
    n_quantiles = n_quantiles
)

mdl_q_cnn$compile(
    optimizer = keras::optimizer_adam(
        learning_rate = best_cnn$learning_rate,
        clipnorm = 0.5
    ),
    loss = pinball_loss
)

x_fit_q_cnn <- select_model_inputs(
    mdl_q_cnn,
    list(
        climate_input = X_fit_aa,
        region_input = region_fit_aa,
        age_input = age_fit_aa,
        age_numeric_input = age_mid_fit_aa
    )
)

history_q_cnn <- keras::fit(
    mdl_q_cnn,
    x = x_fit_q_cnn,
    y = y_fit_aa,
    sample_weight = sample_weights_fit,
    shuffle = TRUE,
    epochs = epochs_q_cnn,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_q,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

time_q_cnn <- difftime(Sys.time(), start_q_cnn, units = "mins")

start_q_gnn <- Sys.time()

mdl_q_gnn <- gnn_lstm_scalar_output_quantile(
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
    use_year_numeric = TRUE,
    n_quantiles = n_quantiles
)

mdl_q_gnn$compile(
    optimizer = keras::optimizer_adam(
        learning_rate = best_gnn$learning_rate,
        clipnorm = 0.5
    ),
    loss = pinball_loss
)

x_fit_q_gnn <- select_model_inputs(
    mdl_q_gnn,
    list(
        climate_input = X_fit_gnn,
        age_input = age_idx_fit_gnn,
        region_input = region_idx_fit_gnn,
        age_numeric_input = age_mid_fit_gnn,
        year_numeric_input = year_numeric_fit_gnn
    )
)

history_q_gnn <- keras::fit(
    mdl_q_gnn,
    x = x_fit_q_gnn,
    y = y_fit_gnn,
    sample_weight = sample_weights_fit_gnn,
    shuffle = TRUE,
    epochs = epochs_q_gnn,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_q,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

mdl_q_gnn$compile(
    optimizer = keras::optimizer_adam(
        learning_rate = best_gnn$learning_rate,
        clipnorm = 0.5
    ),
    loss = pinball_loss
)

time_q_gnn <- difftime(Sys.time(), start_q_gnn, units = "mins")

# -- 8.4A  Quantile MortFCNet -----------------------------------------------
start_q_mort <- Sys.time()

mdl_q_mort <- cnn_mort_fc_net(
    seq_length = seq_length_mort,
    input_size = n_climate_features_mort,
    lstm_units = best_mort$lstm_units,
    hidden_size_1 = best_mort$hidden_size_1,
    hidden_size_2 = best_mort$hidden_size_2,
    hidden_size_3 = best_mort$hidden_size_3,
    output_size = n_quantiles,
    cell_type = "gru",
    use_cnn = FALSE,
    use_dense = TRUE,
    use_age_embedding = FALSE,
    use_age_embedding_pre_lstm = FALSE,
    use_region_embedding = FALSE,
    dropout_rate = best_mort$dropout_rate
)

mdl_q_mort$compile(
    optimizer = keras::optimizer_adam(
        learning_rate = best_mort$learning_rate,
        clipnorm = 0.5
    ),
    loss = pinball_loss
)

history_q_mort <- keras::fit(
    mdl_q_mort,
    x = X_fit_mort,
    y = y_fit_mort,
    shuffle = TRUE,
    epochs = epochs_q_mort,
    batch_size = batch_size,
    callbacks = list(
        keras::callback_early_stopping(
            monitor = "loss",
            patience = patience_q,
            restore_best_weights = TRUE,
            verbose = 1
        )
    ),
    verbose = 2
)

time_q_mort <- difftime(Sys.time(), start_q_mort, units = "mins")

# -- 8.4B  Lee-Carter only baseline (no ML residual, just LC) ---------------
# LC only has point estimates. Uncertainty comes from bootstrap components only.
# We represent LC output as a matrix with the median quantile = LC value,
# and low/high quantiles = LC value (to be combined with bootstrap spreads).
logmu_lc_only_test <- array(
    NA_real_,
    dim = c(n_weeks_test, n_ages, n_regions),
    dimnames = list(as.character(test_dates), age_names, region_names)
)
for (t in seq_len(n_weeks_test)) {
    year_t <- lubridate::isoyear(test_dates[t])
    week_t <- lubridate::isoweek(test_dates[t])
    year_idx <- resolve_year_index(year_t, forecast_years, "LC only test year")
    for (a in seq_len(n_ages)) {
        logmu_lc_only_test[t, a, ] <- lca$beta1[a, ] +
            lca$beta2[a] * kappa2_full[year_idx, ] +
            lca$beta3[a] * lca$kappa3[week_t, ]
    }
}

raw_q_cnn <- sort_quantiles(predict(
    mdl_q_cnn,
    select_model_inputs(
        mdl_q_cnn,
        list(
            climate_input = X_test_aa,
            region_input = region_test_aa,
            age_input = age_test_aa,
            age_numeric_input = age_mid_test_aa
        )
    )
))
q_cnn_test <- unscale_quantile_matrix(
    raw_q_cnn,
    col_names_test_aa,
    y_scale,
    y_center,
    intercept_fix_per_age
)

raw_q_gnn <- sort_quantiles(predict(
    mdl_q_gnn,
    select_model_inputs(
        mdl_q_gnn,
        list(
            climate_input = X_test_gnn,
            age_input = age_idx_test_gnn,
            region_input = region_idx_test_gnn,
            age_numeric_input = age_mid_test_gnn,
            year_numeric_input = year_numeric_test_gnn
        )
    )
))
q_gnn_test <- unscale_quantile_matrix(
    raw_q_gnn,
    col_names_test_gnn,
    y_scale,
    y_center,
    intercept_fix_gnn
)

# -- MortFCNet quantile test predictions
raw_q_mort <- sort_quantiles(predict(mdl_q_mort, X_test_mort))
# MortFCNet outputs region-wide residual quantiles; scale each quantile column
q_mort_test <- sweep(raw_q_mort, 2, y_scale_mort_common, "*") +
    y_center_mort_common

uq_cnn_long <- build_uq_long(
    q_mat = q_cnn_test,
    target_weeks = seqs_age_aware$target_week[test_mask_aa],
    regions_vec = seqs_age_aware$Region[test_mask_aa],
    ages_vec = seqs_age_aware$Age[test_mask_aa],
    baseline_logmxt_wide = baseline_logmxt_wide,
    tau_vec = tau_levels,
    model_label = "CNN-LSTM-Q"
)

uq_gnn_long <- build_uq_long(
    q_mat = q_gnn_test,
    target_weeks = target_week_test_gnn,
    regions_vec = region_name_test_gnn,
    ages_vec = age_name_test_gnn,
    baseline_logmxt_wide = baseline_logmxt_wide,
    tau_vec = tau_levels,
    model_label = "GNN-LSTM-Q"
)

# -- MortFCNet UQ long format (broadcast region residuals to age groups)
test_pred_df_mort_q <- data.frame(
    Week = seqs_age_aware_mort$target_week[test_mask_mort],
    Region = seqs_age_aware_mort$Region[test_mask_mort],
    residual_q05 = q_mort_test[, 1],
    residual_q50 = q_mort_test[, 2],
    residual_q95 = q_mort_test[, 3],
    stringsAsFactors = FALSE
)
# Broadcast to all age groups for this region-week
uq_mort_long_list <- list()
for (age_val in xv) {
    age_str <- as.character(age_val)
    uq_mort_long_list[[age_str]] <- test_pred_df_mort_q %>%
        dplyr::mutate(Age = age_str) %>%
        tidyr::pivot_longer(
            cols = starts_with("residual_"),
            names_to = "quantile_label",
            values_to = "residual_pred"
        ) %>%
        dplyr::mutate(
            tau = as.numeric(sub("^residual_q", "", quantile_label)) / 100,
            AgeNumeric = safe_map_age(Age),
            model = "MortFCNet-Q"
        ) %>%
        dplyr::select(Week, Region, Age, AgeNumeric, tau, residual_pred, model)
}
mort_long_combined <- dplyr::bind_rows(uq_mort_long_list)

baseline_long_mort <- baseline_logmxt_wide %>%
    tidyr::pivot_longer(
        cols = -Week,
        names_to = c("age_str", "Region"),
        names_sep = "_(?=[^_]+$)",
        values_to = "logmxt_baseline"
    ) %>%
    dplyr::mutate(
        Week = as.Date(Week),
        AgeNumeric = age_str
    ) %>%
    dplyr::select(Week, Region, AgeNumeric, logmxt_baseline)

uq_mort_long <- mort_long_combined %>%
    dplyr::mutate(Week = as.Date(Week)) %>%
    dplyr::left_join(
        baseline_long_mort,
        by = c("Week", "Region", "AgeNumeric")
    ) %>%
    dplyr::mutate(
        mxt_baseline = exp(logmxt_baseline),
        mxt_adjusted_q = mxt_baseline + residual_pred,
        logmxt_adjusted_q = log(pmax(mxt_adjusted_q, 1e-12))
    ) %>%
    dplyr::select(
        Week,
        Region,
        Age,
        AgeNumeric,
        tau,
        residual_pred,
        logmxt_adjusted_q,
        model
    )

# -- Lee-Carter only UQ (point logmxt from LC, no residuals; uncertainty only from bootstrap)
uq_lc_long_list <- list()
for (age_val in xv) {
    age_str <- as.character(age_val)
    base_df <- expand.grid(
        Week = as.Date(test_dates),
        Region = region_names,
        Age = age_str,
        tau = tau_levels,
        stringsAsFactors = FALSE
    )
    week_idx <- match(base_df$Week, as.Date(test_dates))
    age_idx <- match(base_df$Age, age_names)
    region_idx <- match(base_df$Region, region_names)
    base_df$logmxt_adjusted_q <- logmu_lc_only_test[cbind(
        week_idx,
        age_idx,
        region_idx
    )]
    base_df$model <- "LC-only"
    base_df$AgeNumeric <- safe_map_age(base_df$Age)
    uq_lc_long_list[[age_str]] <- base_df
}
uq_lc_long <- dplyr::bind_rows(uq_lc_long_list) %>%
    dplyr::select(Week, Region, Age, AgeNumeric, tau, logmxt_adjusted_q, model)

uq_results_long <- dplyr::bind_rows(
    uq_cnn_long,
    uq_gnn_long,
    uq_mort_long,
    uq_lc_long
)

compute_lca_dispersion_array <- function(lca_obj, n_weeks) {
    phi_mat <- exp(outer(lca_obj$phix, lca_obj$phir, "+"))
    array(phi_mat, dim = c(n_weeks, nrow(phi_mat), ncol(phi_mat)))
}

compute_lca_mean_deaths <- lca$mhat * etxr
compute_lca_mean_deaths[!is.finite(compute_lca_mean_deaths)] <- 1e-6
disp_array <- compute_lca_dispersion_array(lca, dim(compute_lca_mean_deaths)[1])

boot_logmu_test_lc <- array(
    NA_real_,
    dim = c(B_boot, n_weeks_test, n_ages, n_regions),
    dimnames = list(NULL, as.character(test_dates), age_names, region_names)
)

ncores_boot <- parallel::detectCores() - 1L
if (!is.finite(ncores_boot) || ncores_boot < 1L) {
    ncores_boot <- 1L
}

cl <- snow::makeCluster(ncores_boot)
doSNOW::registerDoSNOW(cl)
snow::clusterEvalQ(cl, {
    library(dplyr)
    library(MASS)
    library(lubridate)
})
snow::clusterExport(
    cl,
    varlist = c(
        "compute_lca_mean_deaths",
        "disp_array",
        "dtxr",
        "etxr",
        "xv",
        "isoyv",
        "isowv",
        "lca",
        "kappa2_full",
        "test_dates",
        "forecast_years",
        "n_ages",
        "n_regions",
        "ind.rm",
        "fit701M.nb",
        "compute_test_logmu_array",
        "resolve_year_index"
    ),
    envir = environment()
)

pb <- progress::progress_bar$new(
    format = "Bootstrapping = :letter [:bar] :elapsed | eta: :eta",
    total = B_boot,
    width = 60
)
progress_letter <- seq_len(B_boot)
progress <- function(n) {
    pb$tick(tokens = list(letter = progress_letter[n]))
}
opts <- list(progress = progress)

boot_logmu_test_lc_list <- tryCatch(
    {
        foreach::foreach(
            b = seq_len(B_boot),
            .options.snow = opts,
            .combine = "c"
        ) %dopar%
            {
                set.seed(b)
                d_boot <- array(
                    rnbinom(
                        n = prod(dim(compute_lca_mean_deaths)),
                        mu = as.vector(compute_lca_mean_deaths),
                        size = as.vector(disp_array)
                    ),
                    dim = dim(compute_lca_mean_deaths),
                    dimnames = dimnames(compute_lca_mean_deaths)
                )
                d_boot[ind.rm, , ] <- dtxr[ind.rm, , ]
                d_boot[d_boot <= 0] <- 0.001

                lca_b <- tryCatch(
                    fit701M.nb(xv, isoyv, isowv, etxr, d_boot, xv * 0, "ALL"),
                    error = function(e) NULL
                )
                if (is.null(lca_b)) {
                    list(NULL)
                } else {
                    list(compute_test_logmu_array(
                        lca_obj = lca_b,
                        kappa2_full_b = kappa2_full,
                        test_dates = test_dates,
                        forecast_years_b = forecast_years,
                        n_ages = n_ages,
                        n_regions = n_regions
                    ))
                }
            }
    },
    finally = {
        snow::stopCluster(cl)
    }
)

lc_boot_success <- 0L
for (b in seq_len(B_boot)) {
    if (is.null(boot_logmu_test_lc_list[[b]])) {
        next
    }
    boot_logmu_test_lc[b, , , ] <- boot_logmu_test_lc_list[[b]]
    lc_boot_success <- lc_boot_success + 1L
}

kappa_fr_fitted <- coherent_res$kappa_fr
n_train_years <- length(kappa_fr_fitted)
u_fitted_mat <- coherent_res$u_matrix

rw_innovations <- diff(kappa_fr_fitted) - coherent_res$drift
rw_sigma <- stats::sd(rw_innovations, na.rm = TRUE)
if (!is.finite(rw_sigma) || rw_sigma == 0) {
    rw_sigma <- 1e-4
}

ar1_coef <- setNames(numeric(n_regions), region_names)
ar1_sigma <- setNames(numeric(n_regions), region_names)
for (r in seq_len(n_regions)) {
    u_r <- u_fitted_mat[, r]
    ar1_fit <- tryCatch(
        stats::arima(
            u_r,
            order = c(1, 0, 0),
            include.mean = FALSE,
            method = "ML"
        ),
        error = function(e) NULL
    )
    if (!is.null(ar1_fit)) {
        ar1_coef[r] <- as.numeric(ar1_fit$coef["ar1"])
        ar1_sigma[r] <- sqrt(ar1_fit$sigma2)
    } else {
        ar1_coef[r] <- 0
        ar1_sigma[r] <- stats::sd(u_r, na.rm = TRUE)
    }
    if (!is.finite(ar1_sigma[r]) || ar1_sigma[r] == 0) {
        ar1_sigma[r] <- 1e-4
    }
}

boot_logmu_test_kappa <- array(
    NA_real_,
    dim = c(B_boot, n_weeks_test, n_ages, n_regions),
    dimnames = list(NULL, as.character(test_dates), age_names, region_names)
)

kappa_fr_last <- tail(kappa_fr_fitted, 1L)
u_last <- u_fitted_mat[n_train_years, ]

start_kappa_boot <- Sys.time()
for (b in seq_len(B_boot)) {
    if (b %% 100L == 0L) {
        elapsed <- as.numeric(difftime(
            Sys.time(),
            start_kappa_boot,
            units = "secs"
        ))
        rate_per_sec <- b / elapsed
        est_total_secs <- B_boot / rate_per_sec
        est_remaining_secs <- est_total_secs - elapsed
        elapsed
    }

    eps_rw <- rnorm(n_years_forecast, 0, rw_sigma)
    kappa_fr_b <- numeric(n_years_forecast)
    kappa_fr_b[1] <- kappa_fr_last + coherent_res$drift + eps_rw[1]
    if (n_years_forecast > 1L) {
        for (s in 2:n_years_forecast) {
            kappa_fr_b[s] <- kappa_fr_b[s - 1L] + coherent_res$drift + eps_rw[s]
        }
    }

    u_b <- matrix(NA_real_, nrow = n_years_forecast, ncol = n_regions)
    for (s in seq_len(n_years_forecast)) {
        prev_u <- if (s == 1L) u_last else u_b[s - 1L, ]
        eps_ar <- rnorm(n_regions, 0, ar1_sigma)
        u_b[s, ] <- ar1_coef * prev_u + eps_ar
    }

    kappa2_b_forecast <- u_b + kappa_fr_b
    colnames(kappa2_b_forecast) <- region_names
    kappa2_b_full <- rbind(kappa2_fitted, kappa2_b_forecast)
    rownames(kappa2_b_full) <- as.character(
        train_start:(train_end + n_years_forecast)
    )
    forecast_years_b <- as.integer(rownames(kappa2_b_full))

    boot_logmu_test_kappa[b, , , ] <- compute_test_logmu_array(
        lca_obj = lca,
        kappa2_full_b = kappa2_b_full,
        test_dates = test_dates,
        forecast_years_b = forecast_years_b,
        n_ages = n_ages,
        n_regions = n_regions
    )
}

z_lo_sn <- qnorm(0.05)
z_hi_sn <- qnorm(0.95)

extract_tau <- function(uq_long_df, tau_val) {
    uq_long_df %>%
        dplyr::filter(tau == tau_val) %>%
        dplyr::mutate(
            Week = as.Date(Week),
            AgeNumeric = as.character(AgeNumeric)
        ) %>%
        dplyr::select(Week, Region, AgeNumeric, residual_pred)
}

build_splitnorm_params <- function(uq_long_df, model_label) {
    # If the input doesn't contain residual predictions (e.g. LC-only),
    # return an empty params data frame so downstream code treats it as
    # having no split-normal residuals.
    if (!"residual_pred" %in% names(uq_long_df)) {
        return(dplyr::tibble(
            Week = as.Date(character()),
            Region = character(),
            AgeNumeric = character(),
            q05 = double(),
            q50 = double(),
            q95 = double(),
            mu = double(),
            sigma_L = double(),
            sigma_R = double(),
            t_idx = integer(),
            a_idx = integer(),
            r_idx = integer(),
            model = character()
        ))
    }

    q05 <- extract_tau(uq_long_df, 0.05)
    q50 <- extract_tau(uq_long_df, 0.50)
    q95 <- extract_tau(uq_long_df, 0.95)

    q50 %>%
        dplyr::rename(q50 = residual_pred) %>%
        dplyr::inner_join(
            q05 %>% dplyr::rename(q05 = residual_pred),
            by = c("Week", "Region", "AgeNumeric")
        ) %>%
        dplyr::inner_join(
            q95 %>% dplyr::rename(q95 = residual_pred),
            by = c("Week", "Region", "AgeNumeric")
        ) %>%
        dplyr::mutate(
            mu = q50,
            sigma_L = pmax((mu - q05) / (-z_lo_sn), 1e-10),
            sigma_R = pmax((q95 - mu) / z_hi_sn, 1e-10),
            t_idx = match(Week, test_dates),
            a_idx = match(AgeNumeric, age_names),
            r_idx = match(Region, region_names),
            model = model_label
        ) %>%
        dplyr::filter(!is.na(t_idx), !is.na(a_idx), !is.na(r_idx))
}

sample_splitnorm_matrix <- function(params_df, B = B_boot) {
    n_cells <- nrow(params_df)
    smp <- matrix(NA_real_, nrow = n_cells, ncol = B)
    for (i in seq_len(n_cells)) {
        mu <- params_df$mu[i]
        sigma_L <- params_df$sigma_L[i]
        sigma_R <- params_df$sigma_R[i]
        left_half <- mu - abs(stats::rnorm(B, 0, sigma_L))
        right_half <- mu + abs(stats::rnorm(B, 0, sigma_R))
        mask <- stats::runif(B) < 0.5
        smp[i, ] <- ifelse(mask, left_half, right_half)
    }
    smp
}

sn_params_cnn <- build_splitnorm_params(uq_cnn_long, "CNN-LSTM")
sn_params_gnn <- build_splitnorm_params(uq_gnn_long, "GNN-LSTM")

sn_samples_cnn <- sample_splitnorm_matrix(sn_params_cnn)

sn_samples_gnn <- sample_splitnorm_matrix(sn_params_gnn)

sn_params_mort <- build_splitnorm_params(uq_mort_long, "MortFCNet")
sn_params_lc <- build_splitnorm_params(uq_lc_long, "LC-only")

sn_samples_mort <- sample_splitnorm_matrix(sn_params_mort)

sn_samples_lc <- sample_splitnorm_matrix(sn_params_lc)

probs_out <- c(0.05, 0.50, 0.95)

compute_combined_empirical_quantiles <- function(
    boot_lc,
    boot_kappa,
    logmu_point,
    sn_params,
    sn_samples,
    probs = probs_out
) {
    B <- dim(boot_lc)[1]
    Tt <- dim(boot_lc)[2]
    A <- dim(boot_lc)[3]
    R <- dim(boot_lc)[4]

    q_out <- lapply(setNames(probs, paste0("q", probs * 100)), function(p) {
        array(NA_real_, dim = c(Tt, A, R))
    })

    has_sn <- !is.null(sn_params) && nrow(sn_params) > 0L
    cell_lookup <- array(NA_integer_, dim = c(Tt, A, R))
    if (has_sn) {
        for (i in seq_len(nrow(sn_params))) {
            cell_lookup[
                sn_params$t_idx[i],
                sn_params$a_idx[i],
                sn_params$r_idx[i]
            ] <- i
        }
    }

    for (a in seq_len(A)) {
        for (r in seq_len(R)) {
            lc_point_ar <- logmu_point[, a, r]
            d_param <- boot_lc[,, a, r] -
                matrix(lc_point_ar, nrow = B, ncol = Tt, byrow = TRUE)
            d_kappa <- boot_kappa[,, a, r] -
                matrix(lc_point_ar, nrow = B, ncol = Tt, byrow = TRUE)
            logmu_lc_b <- matrix(
                lc_point_ar,
                nrow = B,
                ncol = Tt,
                byrow = TRUE
            ) +
                d_param +
                d_kappa

            for (t in seq_len(Tt)) {
                row_i <- if (has_sn) cell_lookup[t, a, r] else NA_integer_
                R_b <- if (has_sn && !is.na(row_i)) {
                    sn_samples[row_i, ]
                } else {
                    rep(0, B)
                }
                mxt_total_b <- exp(logmu_lc_b[, t]) + R_b
                logmu_total_b <- log(pmax(mxt_total_b, 1e-12))
                for (k in seq_along(probs)) {
                    q_out[[k]][t, a, r] <- stats::quantile(
                        logmu_total_b,
                        probs[k],
                        na.rm = TRUE
                    )
                }
            }
        }
    }

    q_out
}

logmu_point_test <- logmxt_test

q_combined_cnn <- compute_combined_empirical_quantiles(
    boot_lc = boot_logmu_test_lc,
    boot_kappa = boot_logmu_test_kappa,
    logmu_point = logmu_point_test,
    sn_params = sn_params_cnn,
    sn_samples = sn_samples_cnn
)

q_combined_gnn <- compute_combined_empirical_quantiles(
    boot_lc = boot_logmu_test_lc,
    boot_kappa = boot_logmu_test_kappa,
    logmu_point = logmu_point_test,
    sn_params = sn_params_gnn,
    sn_samples = sn_samples_gnn
)

q_combined_mort <- compute_combined_empirical_quantiles(
    boot_lc = boot_logmu_test_lc,
    boot_kappa = boot_logmu_test_kappa,
    logmu_point = logmu_point_test,
    sn_params = sn_params_mort,
    sn_samples = sn_samples_mort
)

q_combined_lc <- compute_combined_empirical_quantiles(
    boot_lc = boot_logmu_test_lc,
    boot_kappa = boot_logmu_test_kappa,
    logmu_point = logmu_point_test,
    sn_params = sn_params_lc,
    sn_samples = sn_samples_lc
)

quantile_arrays_to_long <- function(
    q_list,
    test_dates,
    age_names,
    region_names,
    model_label
) {
    pnames <- names(q_list)
    base <- expand.grid(
        t_idx = seq_len(dim(q_list[[1]])[1]),
        a_idx = seq_len(dim(q_list[[1]])[2]),
        r_idx = seq_len(dim(q_list[[1]])[3]),
        KEEP.OUT.ATTRS = FALSE
    ) %>%
        dplyr::mutate(
            Week = as.Date(test_dates[t_idx]),
            age = as.integer(age_names[a_idx]),
            region = region_names[r_idx],
            model = model_label
        )
    for (k in seq_along(q_list)) {
        base[[pnames[k]]] <- as.vector(q_list[[k]])
    }
    base %>%
        dplyr::select(-t_idx, -a_idx, -r_idx) %>%
        dplyr::mutate(
            mxt_q5 = exp(q5),
            mxt_q50 = exp(q50),
            mxt_q95 = exp(q95),
            interval_width_log = q95 - q5
        )
}

uq_combined_cnn_long <- quantile_arrays_to_long(
    q_combined_cnn,
    test_dates,
    age_names,
    region_names,
    "CNN-LSTM"
)
uq_combined_gnn_long <- quantile_arrays_to_long(
    q_combined_gnn,
    test_dates,
    age_names,
    region_names,
    "GNN-LSTM"
)
uq_combined_mort_long <- quantile_arrays_to_long(
    q_combined_mort,
    test_dates,
    age_names,
    region_names,
    "MortFCNet"
)
uq_combined_lc_long <- quantile_arrays_to_long(
    q_combined_lc,
    test_dates,
    age_names,
    region_names,
    "LC-only"
)
uq_combined_long <- dplyr::bind_rows(
    uq_combined_cnn_long,
    uq_combined_gnn_long,
    uq_combined_mort_long,
    uq_combined_lc_long
) %>%
    dplyr::mutate(Week = as.Date(Week))

obs_for_coverage <- observed_logmxt_wide %>%
    tidyr::pivot_longer(
        cols = -Week,
        names_to = c("age_str", "region"),
        names_sep = "_(?=[^_]+$)",
        values_to = "logmxt_obs"
    ) %>%
    dplyr::mutate(
        Week = as.Date(Week),
        age = as.integer(age_str)
    ) %>%
    dplyr::filter(
        lubridate::isoyear(Week) %in% test_years,
        is.finite(logmxt_obs)
    ) %>%
    dplyr::select(Week, age, region, logmxt_obs)

coverage_combined <- uq_combined_long %>%
    dplyr::left_join(obs_for_coverage, by = c("Week", "age", "region")) %>%
    dplyr::filter(is.finite(logmxt_obs)) %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
        n_obs = dplyr::n(),
        empirical_coverage = mean(
            logmxt_obs >= q5 & logmxt_obs <= q95,
            na.rm = TRUE
        ),
        mean_width_log = mean(interval_width_log, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        nominal_coverage = 0.90,
        coverage_gap = empirical_coverage - nominal_coverage
    )

# Store combined coverage metrics for later analysis
uq_results <- list(
    "CNN-LSTM-Q" = list(
        model_name = "CNN-LSTM-Q",
        tau_levels = tau_levels,
        predictions_long = uq_cnn_long,
        combined_quantiles = uq_combined_cnn_long,
        time_minutes = as.numeric(time_q_cnn)
    ),
    "GNN-LSTM-Q" = list(
        model_name = "GNN-LSTM-Q",
        tau_levels = tau_levels,
        predictions_long = uq_gnn_long,
        combined_quantiles = uq_combined_gnn_long,
        time_minutes = as.numeric(time_q_gnn)
    ),
    "MortFCNet-Q" = list(
        model_name = "MortFCNet-Q",
        tau_levels = tau_levels,
        predictions_long = uq_mort_long,
        combined_quantiles = uq_combined_mort_long,
        time_minutes = as.numeric(time_q_mort)
    ),
    "LC-only" = list(
        model_name = "LC-only",
        predictions_long = uq_lc_long,
        combined_quantiles = uq_combined_lc_long
    )
)

################################################################################
# PREDICTION INTERVAL PERFORMANCE: PICP AND MPIW
################################################################################
interval_plot_dir <- file.path("Results", "plots", "prediction_intervals")
interval_table_dir <- file.path("Results", "tables")
dir.create(interval_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interval_table_dir, recursive = TRUE, showWarnings = FALSE)

prediction_interval_table <- uq_combined_long %>%
    dplyr::left_join(obs_for_coverage, by = c("Week", "age", "region")) %>%
    dplyr::filter(is.finite(logmxt_obs)) %>%
    dplyr::mutate(mxt_obs = exp(logmxt_obs)) %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
        n_obs = dplyr::n(),
        PICP = mean(mxt_obs >= mxt_q5 & mxt_obs <= mxt_q95, na.rm = TRUE),
        MPIW = mean(mxt_q95 - mxt_q5, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        nominal_coverage = 0.90,
        coverage_gap = PICP - nominal_coverage
    ) %>%
    dplyr::arrange(match(
        model,
        c("CNN-LSTM", "GNN-LSTM", "MortFCNet", "LC-only")
    ))

# Save prediction interval evaluation results to disk
write.csv(
    prediction_interval_table,
    file = file.path(interval_table_dir, "prediction_interval_performance.csv"),
    row.names = FALSE
)

picp_plot <- ggplot2::ggplot(
    prediction_interval_table,
    ggplot2::aes(x = model, y = PICP, fill = model)
) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9) +
    ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed", color = "red") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
        title = "Prediction Interval Coverage Probability (PICP)",
        x = NULL,
        y = "PICP"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")

mpiw_plot <- ggplot2::ggplot(
    prediction_interval_table,
    ggplot2::aes(x = model, y = MPIW, fill = model)
) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9) +
    ggplot2::labs(
        title = "Mean Prediction Interval Width (MPIW)",
        x = NULL,
        y = "MPIW"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")

################################################################################
# MODEL-SPECIFIC UNCERTAINTY PLOTS FOR BEST AGE GROUP / BEST REGION
################################################################################
uncertainty_plot_dir <- file.path("Results", "plots", "uncertainty_profiles")
dir.create(uncertainty_plot_dir, recursive = TRUE, showWarnings = FALSE)

observed_long_all <- observed_logmxt_wide %>%
    tidyr::pivot_longer(
        cols = -Week,
        names_to = c("age_str", "region"),
        names_sep = "_(?=[^_]+$)",
        values_to = "logmxt_obs"
    ) %>%
    dplyr::mutate(
        Week = as.Date(Week),
        age = as.integer(age_str),
        region = as.character(region)
    ) %>%
    dplyr::select(Week, age, region, logmxt_obs)

get_best_selection <- function(model_name) {
    # Fixed selection: 90+ age group and region FRG0
    age_mse <- per_age_metrics_table %>%
        dplyr::filter(model == model_name, age == 90) %>%
        dplyr::pull(mse_mxt_test)
    region_mse <- per_region_metrics_table %>%
        dplyr::filter(model == model_name, region == "FRG0") %>%
        dplyr::pull(mse_mxt_test)

    list(
        best_age = 90L,
        best_region = "FRG0",
        best_age_mse = if (length(age_mse) > 0) age_mse[1] else NA_real_,
        best_region_mse = if (length(region_mse) > 0) {
            region_mse[1]
        } else {
            NA_real_
        }
    )
}

build_uncertainty_profile <- function(
    model_name,
    selection_type,
    selection_value
) {
    if (selection_type == "age") {
        profile_long <- uq_combined_long %>%
            dplyr::filter(model == model_name, age == selection_value)
        panel_label <- sprintf("Best age group: %s", selection_value)
        panel_note <- "Averaged across regions"
    } else if (selection_type == "region") {
        profile_long <- uq_combined_long %>%
            dplyr::filter(model == model_name, region == selection_value)
        panel_label <- sprintf("Best region: %s", selection_value)
        panel_note <- "Averaged across age groups"
    } else {
        stop("selection_type must be 'age' or 'region'.")
    }
    profile_long %>%
        dplyr::left_join(observed_long_all, by = c("Week", "age", "region")) %>%
        dplyr::group_by(Week) %>%
        dplyr::summarise(
            q5 = mean(q5, na.rm = TRUE),
            q50 = mean(q50, na.rm = TRUE),
            q95 = mean(q95, na.rm = TRUE),
            logmxt_obs = mean(logmxt_obs, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::arrange(Week) %>%
        dplyr::mutate(
            model = model_name,
            panel = panel_label,
            panel_note = panel_note
        )
}

plot_model_uncertainty <- function(model_name) {
    selection <- get_best_selection(model_name)
    age_profile <- build_uncertainty_profile(
        model_name,
        "age",
        selection$best_age
    )
    region_profile <- build_uncertainty_profile(
        model_name,
        "region",
        selection$best_region
    )
    plot_df <- dplyr::bind_rows(age_profile, region_profile)

    title_text <- sprintf(
        "%s combined uncertainty vs observed mortality",
        model_name
    )
    subtitle_text <- sprintf(
        "Best age group = %s (test MSE %.4f); best region = %s (test MSE %.4f). Shaded band combines quantile LSTM, bootstrap, and kappa simulations.",
        selection$best_age,
        selection$best_age_mse,
        selection$best_region,
        selection$best_region_mse
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Week)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = q5, ymax = q95),
            fill = "#2A6F97",
            alpha = 0.20
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = q50),
            color = "#2A6F97",
            linewidth = 0.8
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = logmxt_obs),
            color = "#111111",
            linewidth = 0.7
        ) +
        ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
        ggplot2::labs(
            title = title_text,
            subtitle = subtitle_text,
            x = "Week",
            y = "Log mortality"
        ) +
        ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            strip.text = ggplot2::element_text(face = "bold"),
            legend.position = "none"
        )

    out_file <- file.path(
        uncertainty_plot_dir,
        paste0(
            gsub("[^A-Za-z0-9]+", "_", model_name),
            "_best_age_region_uncertainty.png"
        )
    )
    invisible(list(selection = selection, plot = p, file = out_file))
}

uncertainty_plot_summaries <- list(
    "CNN-LSTM" = plot_model_uncertainty("CNN-LSTM"),
    "GNN-LSTM" = plot_model_uncertainty("GNN-LSTM"),
    "MortFCNet" = plot_model_uncertainty("MortFCNet")
)

uncertainty_selection_table <- dplyr::bind_rows(
    lapply(names(uncertainty_plot_summaries), function(model_name) {
        sel <- uncertainty_plot_summaries[[model_name]]$selection
        data.frame(
            model = model_name,
            best_age = sel$best_age,
            best_age_mse = sel$best_age_mse,
            best_region = sel$best_region,
            best_region_mse = sel$best_region_mse,
            stringsAsFactors = FALSE
        )
    })
)

# Example-style plots will be generated below for the selected best age and region
################################################################################
# EXAMPLE-STYLE QUANTILE FORECAST PLOTS AT BEST AGE / REGION
################################################################################

save_best_combined_uq_plot <- function(
    model_name,
    uq_long,
    best_age,
    best_region
) {
    plot_df <- uq_long %>%
        dplyr::filter(
            model == model_name,
            age == best_age,
            region == best_region
        ) %>%
        dplyr::left_join(observed_long_all, by = c("Week", "age", "region")) %>%
        dplyr::arrange(Week)

    if (nrow(plot_df) == 0L) {
        warning(sprintf(
            "No rows for %s | age %s | region %s.",
            model_name,
            best_age,
            best_region
        ))
        return(invisible(NULL))
    }

    # ── Per-plot tight y-limits (the actual fix) ──────────────────────────────
    local_limits <- plot_df %>%
        dplyr::summarise(
            y_min = min(c(q5, logmxt_obs), na.rm = TRUE),
            y_max = max(c(q95, logmxt_obs), na.rm = TRUE)
        ) %>%
        dplyr::mutate(
            padding = 0.08 * (y_max - y_min),
            y_min = y_min - padding,
            y_max = y_max + padding
        )

    model_color <- dplyr::case_when(
        grepl("GNN", model_name, ignore.case = TRUE) ~ "#7B1FA2",
        grepl("CNN", model_name, ignore.case = TRUE) ~ "#d95f02",
        TRUE ~ "#1b9e77"
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Week)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = q5, ymax = q95, fill = "90% interval"),
            alpha = 0.20
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = q50, color = "Forecast median"),
            linewidth = 0.85
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = logmxt_obs, color = "Observed"),
            linewidth = 0.65,
            alpha = 0.60,
            na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(
            values = c("Forecast median" = model_color, "Observed" = "#333333"),
            guide = ggplot2::guide_legend(order = 1)
        ) +
        ggplot2::scale_fill_manual(
            values = c("90% interval" = model_color),
            guide = ggplot2::guide_legend(order = 2)
        ) +
        ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
        ggplot2::scale_y_continuous(
            limits = c(local_limits$y_min, local_limits$y_max),
            expand = ggplot2::expansion(0)
        ) +
        ggplot2::labs(
            title = sprintf("%s \u2014 Forecast Uncertainty", model_name),
            subtitle = sprintf(
                "Region: %s  |  Age group: %s",
                best_region,
                best_age
            ),
            x = NULL,
            y = "Log Mortality",
            color = NULL,
            fill = NULL
        ) +
        ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = 13,
                hjust = 0.5
            ),
            plot.title.position = "plot",
            plot.subtitle = ggplot2::element_text(
                size = 10,
                colour = "grey40",
                hjust = 0.5
            ),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.width = ggplot2::unit(1.4, "lines"),
            legend.text = ggplot2::element_text(size = 9),
            legend.spacing.x = ggplot2::unit(0.5, "lines"),
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1,
                size = 9
            ),
            axis.text.y = ggplot2::element_text(size = 9),
            axis.title.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 8)
            ),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                colour = "grey88",
                linewidth = 0.4
            ),
            plot.margin = ggplot2::margin(12, 16, 10, 12)
        )

    out_file <- file.path(
        uncertainty_plot_dir,
        paste0(
            gsub("[^A-Za-z0-9]+", "_", model_name),
            "_best_age_region_forecast.png"
        )
    )
    ggplot2::ggsave(out_file, plot = p, width = 10, height = 5.5, dpi = 300)
    invisible(out_file)
}

purrr::pwalk(
    list(
        model_name = uncertainty_selection_table$model,
        uq_long = list(
            uq_combined_cnn_long,
            uq_combined_gnn_long,
            uq_combined_mort_long
        ),
        best_age = uncertainty_selection_table$best_age,
        best_region = uncertainty_selection_table$best_region
    ),
    save_best_combined_uq_plot
)

################################################################################
# SIMPLE FRG0 YEAR-BY-YEAR IN-SAMPLE FIT PLOTS
################################################################################
frg0_plot_dir <- file.path("Results", "plots", "frg0_yearly_in_sample")
frg0_series_table <- file.path(
    "Results",
    "tables",
    "frg0_yearly_in_sample_series.csv"
)
frg0_index_table <- file.path(
    "Results",
    "tables",
    "frg0_yearly_in_sample_plots_index.csv"
)
dir.create(frg0_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(frg0_series_table), recursive = TRUE, showWarnings = FALSE)

frg0_year_start <- 1990L
frg0_year_end <- 2014L
frg0_years <- frg0_year_start:frg0_year_end
frg0_region <- "FRG0"

aggregate_region_weekly_logmxt <- function(
    logmxt_wide,
    exposure_wide,
    region_value,
    years_keep
) {
    region_cols <- get_region_metric_cols(logmxt_wide, region_value)
    if (length(region_cols) == 0L) {
        stop(
            sprintf("No columns found for region %s.", region_value),
            call. = FALSE
        )
    }

    log_long <- logmxt_wide %>%
        dplyr::select(Week, dplyr::all_of(region_cols)) %>%
        tidyr::pivot_longer(
            cols = -Week,
            names_to = "series_key",
            values_to = "logmxt"
        )

    expo_long <- exposure_wide %>%
        dplyr::select(Week, dplyr::all_of(region_cols)) %>%
        tidyr::pivot_longer(
            cols = -Week,
            names_to = "series_key",
            values_to = "expo"
        )

    dplyr::inner_join(log_long, expo_long, by = c("Week", "series_key")) %>%
        dplyr::mutate(
            Week = as.Date(Week),
            ISOYear = lubridate::isoyear(Week),
            ISOWeek = lubridate::isoweek(Week)
        ) %>%
        dplyr::filter(ISOYear %in% years_keep) %>%
        dplyr::group_by(Week, ISOYear, ISOWeek) %>%
        dplyr::summarise(
            total_expo = sum(expo[is.finite(expo)], na.rm = TRUE),
            total_mxt = sum(exp(logmxt) * expo, na.rm = TRUE),
            log_mortality = ifelse(
                total_expo > 0,
                log(total_mxt / total_expo),
                NA_real_
            ),
            .groups = "drop"
        ) %>%
        dplyr::arrange(Week)
}

frg0_obs_train <- aggregate_region_weekly_logmxt(
    observed_train,
    exposure_train,
    frg0_region,
    frg0_years
)
frg0_cnn_train <- aggregate_region_weekly_logmxt(
    all_results[["CNN-LSTM"]]$adjusted_train_wide,
    exposure_train,
    frg0_region,
    frg0_years
)
frg0_gnn_train <- aggregate_region_weekly_logmxt(
    all_results[["GNN-LSTM"]]$adjusted_train_wide,
    exposure_train,
    frg0_region,
    frg0_years
)

frg0_yearly_series <- dplyr::bind_rows(
    frg0_obs_train %>% dplyr::mutate(model = "Observed"),
    frg0_cnn_train %>% dplyr::mutate(model = "CNN-LSTM"),
    frg0_gnn_train %>% dplyr::mutate(model = "GNN-LSTM")
) %>%
    dplyr::select(Week, ISOYear, ISOWeek, model, log_mortality) %>%
    dplyr::arrange(ISOYear, ISOWeek, model)

write.csv(frg0_yearly_series, frg0_series_table, row.names = FALSE)

save_frg0_year_plot <- function(year_value) {
    year_df <- frg0_yearly_series %>%
        dplyr::filter(ISOYear == year_value) %>%
        dplyr::mutate(
            model = factor(
                model,
                levels = c("Observed", "CNN-LSTM", "GNN-LSTM")
            )
        )

    if (nrow(year_df) == 0L) {
        return(NA_character_)
    }

    finite_vals <- year_df$log_mortality[is.finite(year_df$log_mortality)]
    if (length(finite_vals) == 0L) {
        return(NA_character_)
    }

    y_min <- min(finite_vals)
    y_max <- max(finite_vals)
    y_pad <- 0.08 * max(y_max - y_min, 1e-6)

    p <- ggplot2::ggplot(
        year_df,
        ggplot2::aes(
            x = ISOWeek,
            y = log_mortality,
            colour = model,
            group = model
        )
    ) +
        ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
        ggplot2::scale_colour_manual(
            values = c(
                "Observed" = "#111111",
                "CNN-LSTM" = "#D95F02",
                "GNN-LSTM" = "#1B9E77"
            )
        ) +
        ggplot2::scale_x_continuous(
            breaks = seq(1, 52, by = 4),
            limits = c(1, 52)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(y_min - y_pad, y_max + y_pad),
            expand = ggplot2::expansion(mult = 0)
        ) +
        ggplot2::labs(
            title = sprintf(
                "FRG0 in-sample weekly fit by year: %d",
                year_value
            ),
            subtitle = "Observed vs CNN-LSTM vs GNN-LSTM",
            x = "ISO week",
            y = "Log mortality",
            colour = NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "top",
            plot.title = ggplot2::element_text(face = "bold"),
            axis.text.x = ggplot2::element_text(size = 9),
            axis.text.y = ggplot2::element_text(size = 9)
        )

    out_file <- file.path(
        frg0_plot_dir,
        sprintf("frg0_in_sample_weekly_%d.png", year_value)
    )
    ggplot2::ggsave(out_file, plot = p, width = 10, height = 5.5, dpi = 300)
    out_file
}

frg0_year_plot_files <- vapply(
    frg0_years,
    save_frg0_year_plot,
    character(1)
)

frg0_year_plot_index <- data.frame(
    year = frg0_years,
    plot_file = frg0_year_plot_files,
    stringsAsFactors = FALSE
)

write.csv(
    frg0_year_plot_index,
    frg0_index_table,
    row.names = FALSE
)
