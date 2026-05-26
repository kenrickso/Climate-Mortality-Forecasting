# Compute MAE and MSE between predicted and observed mortality.
compute_metrics <- function(pred_wide, obs_wide) {
    if (is.null(pred_wide) || is.null(obs_wide)) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            mse_log = NA_real_,
            mse_mxt = NA_real_
        ))
    }
    pred_wide$Week <- as.Date(pred_wide$Week)
    obs_wide$Week <- as.Date(obs_wide$Week)

    common_weeks <- sort(intersect(pred_wide$Week, obs_wide$Week))
    if (length(common_weeks) == 0) {
        return(list(
            mae_log = NA_real_,
            mae_mxt = NA_real_,
            mse_log = NA_real_,
            mse_mxt = NA_real_
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
            mse_log = NA_real_,
            mse_mxt = NA_real_
        ))
    }

    obs_mat <- as.matrix(obs_sub[, common_cols, drop = FALSE])
    pred_mat <- as.matrix(pred_sub[, common_cols, drop = FALSE])

    if (!all(dim(obs_mat) == dim(pred_mat))) {
        nrow_min <- min(nrow(obs_mat), nrow(pred_mat))
        obs_mat <- obs_mat[seq_len(nrow_min), , drop = FALSE]
        pred_mat <- pred_mat[seq_len(nrow_min), , drop = FALSE]
    }

    diff_log <- pred_mat - obs_mat
    abs_diff_log <- abs(diff_log)
    mae_log <- mean(abs_diff_log[is.finite(abs_diff_log)], na.rm = TRUE)

    mse_log <- mean(diff_log[is.finite(diff_log)]^2, na.rm = TRUE)

    obs_mxt <- exp(obs_mat)
    pred_mxt <- exp(pred_mat)
    diff_mxt <- pred_mxt - obs_mxt
    abs_diff_mxt <- abs(diff_mxt)
    mae_mxt <- mean(abs_diff_mxt[is.finite(abs_diff_mxt)], na.rm = TRUE)

    mse_mxt <- mean(diff_mxt[is.finite(diff_mxt)]^2, na.rm = TRUE)

    list(
        mae_log = mae_log,
        mae_mxt = mae_mxt,
        mse_log = mse_log,
        mse_mxt = mse_mxt
    )
}

# Function: predict_logmxt
# Kept as used by Final.R via mapply. It expects required objects to be provided
# in the calling environment (or available globally) as in the original pipeline.
predict_logmxt <- function(
    year,
    week,
    age_idx,
    region_idx,
    beta1,
    beta2,
    beta3,
    kappa2_full,
    kappa3
) {
    year_idx <- year - train_start + 1
    if (year_idx < 1 || year_idx > nrow(kappa2_full)) {
        return(NA_real_)
    }
    beta1[age_idx, region_idx] +
        beta2[age_idx] * kappa2_full[year_idx, region_idx] +
        beta3[age_idx] * kappa3[week, region_idx]
}

# Helper: map age labels to numeric ages. Used in Final.R for Age -> numeric mapping.
safe_map_age <- function(a) {
    a_chr <- as.character(a)
    res <- vapply(
        a_chr,
        function(x) {
            if (is.na(x) || nzchar(trimws(x)) == FALSE) {
                return(NA_character_)
            }
            xt <- trimws(x)
            n <- suppressWarnings(as.numeric(xt))
            if (!is.na(n)) {
                n_int <- as.integer(n)
                if (n_int >= 90L) {
                    n_int <- 90L
                }
                return(as.character(n_int))
            }
            if (grepl("^\\d+\\s*\\+$", xt)) {
                n2 <- as.numeric(gsub("\\s*\\+$", "", xt))
                if (!is.na(n2)) {
                    n2_int <- as.integer(n2)
                    if (n2_int >= 90L) {
                        n2_int <- 90L
                    }
                    return(as.character(n2_int))
                }
            }
            if (grepl("^\\s*\\d+\\s*-\\s*\\d+", xt)) {
                m <- regmatches(xt, regexpr("\\d+", xt))
                if (length(m) && nzchar(m)) {
                    return(as.character(as.integer(as.numeric(m))))
                }
            }
            m <- regmatches(xt, regexpr("\\d+", xt))
            if (length(m) && nzchar(m)) {
                return(as.character(as.integer(as.numeric(m))))
            }
            xt
        },
        FUN.VALUE = character(1)
    )
    unname(res)
}

# ------------------------------------------------------------------------------
# Reading the EOBS Files

# Utility to list and read E-OBS region files like: FR_NUTS2_xx_daily_2013-2024.rds

# List available variables (the "xx") for files matching COUNTRY_NUTS_XX_daily_2013-2024.rds

# ------------------------------------------------------------------------------
# Function: list_eobs_files
# Description: Lists all raw E-OBS parameter data files for a specified country and NUTS level.
# Inputs:
#   - path: Base directory.
#   - country: The character code indicating the country.
#   - nuts: The NUTS level.
# Output:
#   - A named list containing paths to the matching RDS files.
# ------------------------------------------------------------------------------
list_eobs_files <- function(
    path = "Data/E-OBS",
    country = "FR",
    nuts = "NUTS2"
) {
    # Recursively list .rds files under `path` and group by variable code
    all_files <- list.files(
        path = path,
        pattern = "\\.rds$",
        full.names = TRUE,
        recursive = TRUE
    )
    if (length(all_files) == 0) {
        return(list())
    }

    # Keep only files that start with the expected country and nuts prefix
    keep <- grepl(paste0("^", country, "_", nuts), basename(all_files))
    all_files <- all_files[keep]
    if (length(all_files) == 0) {
        return(list())
    }

    # Extract variable code using the convention: ..._<var>_daily...rds
    vars <- sub(".*_([A-Za-z0-9]+)_daily.*\\.rds$", "\\1", basename(all_files))
    # Only keep those where the substitution changed the string (i.e., matched)
    matched <- vars != basename(all_files)
    all_files <- all_files[matched]
    vars <- vars[matched]
    if (length(all_files) == 0) {
        return(list())
    }

    # Group into a named list: name = variable code, value = character vector of file paths
    res <- split(all_files, vars)
    res
}

# Read one or more E-OBS region files.
# - vars: character vector of variable codes (e.g. "hu","tg"); if NULL reads all available.
# - combine: if TRUE returns a single data.frame with an added "variable" column; else returns a named list of data.frames.

# ------------------------------------------------------------------------------
# Function: read_eobs_files
# Description: Reads E-OBS data files for specified variables, optionally combining them into a single data frame.
# Inputs:
#   - path: Base directory.
#   - country: The character code indicating the country.
#   - nuts: The NUTS level string.
#   - vars: Optional character vector of E-OBS variable codes to read.
#   - combine: Logical flag indicating whether to bind all variables into one data block.
# Output:
#   - Depending on combine, either a named list of data frames or a single combined data frame.
# ------------------------------------------------------------------------------
read_eobs_files <- function(
    path = "Data/E-OBS",
    country = "FR",
    nuts = "NUTS2",
    vars = NULL,
    combine = FALSE
) {
    files_named <- list_eobs_files(path = path, country = country, nuts = nuts)
    if (length(files_named) == 0) {
        warning("No matching E-OBS files found in ", path)
        return(if (combine) data.frame() else list())
    }
    if (!is.null(vars)) {
        # keep only requested vars that exist
        missing_vars <- setdiff(vars, names(files_named))
        if (length(missing_vars)) {
            warning(
                "Requested vars not found and will be ignored: ",
                paste(missing_vars, collapse = ", ")
            )
        }
        files_named <- files_named[intersect(vars, names(files_named))]
    }

    # read files and normalize each data.frame: single climate column -> "value", ensure Date is Date, add variable column
    read_list <- lapply(names(files_named), function(v) {
        paths <- files_named[[v]]
        # read and bind all period-split files for this variable
        dfs_var <- lapply(paths, function(p) {
            dat_raw <- readRDS(p)
            dat <- as.data.frame(dat_raw)

            # normalize Date and Region column names (case-insensitive common variants)
            nm <- names(dat)
            date_idx <- which(tolower(nm) %in% c("date", "day"))
            region_idx <- which(
                tolower(nm) %in% c("region", "nuts_id", "nuts", "region_id")
            )
            if (length(date_idx) == 1) {
                names(dat)[date_idx] <- "Date"
            }
            if (length(region_idx) == 1) {
                names(dat)[region_idx] <- "Region"
            }
            if (!"Date" %in% names(dat)) {
                # try to infer Date from rownames if possible (leave as is otherwise)
                # no strong fallback; keep as-is
            } else {
                dat$Date <- as.Date(dat$Date)
            }

            # attach a standard column name for the variable code
            dat$variable <- v

            # identify candidate climate columns (everything except Date, Region, variable)
            climate_cols <- setdiff(names(dat), c("Date", "Region", "variable"))

            # If the file uses region codes as column names (e.g. ES11, ES12...) and there
            # is no explicit `Region` column, pivot those columns into long format so that
            # we produce `Region` + `value` rows instead of columns full of NAs.
            if (
                length(climate_cols) > 1 &&
                    !("Region" %in% names(dat)) &&
                    all(grepl("^[A-Z]{2}[0-9]{1,3}$", climate_cols))
            ) {
                # stack the climate columns into a long form
                n <- nrow(dat)
                vals <- unlist(dat[climate_cols], use.names = FALSE)
                # replicate non-climate columns for each climate column
                base_cols <- setdiff(names(dat), climate_cols)
                base_rep <- dat[
                    rep(seq_len(n), times = length(climate_cols)),
                    base_cols,
                    drop = FALSE
                ]
                base_rep$Region <- rep(climate_cols, each = n)
                base_rep$value <- vals
                # ensure Date stays Date
                if ("Date" %in% names(base_rep)) {
                    base_rep$Date <- as.Date(base_rep$Date)
                }
                # keep variable column
                base_rep$variable <- v
                dat <- base_rep
                climate_cols <- "value"
            }

            # choose the climate column robustly:
            chosen <- NULL
            if (length(climate_cols) == 0) {
                dat$value <- NA
            } else if (length(climate_cols) == 1) {
                chosen <- climate_cols[1]
                names(dat)[names(dat) == chosen] <- "value"
            } else {
                # prefer a column named exactly as the variable code
                if (v %in% climate_cols) {
                    chosen <- v
                } else {
                    # prefer numeric columns and pick the one with the fewest NAs (most data)
                    numeric_cols <- climate_cols[sapply(
                        dat[climate_cols],
                        is.numeric
                    )]
                    if (length(numeric_cols) > 0) {
                        na_counts <- sapply(dat[numeric_cols], function(col) {
                            sum(is.na(col))
                        })
                        chosen <- numeric_cols[which.min(na_counts)]
                    } else {
                        # fallback: choose the first climate col
                        chosen <- climate_cols[1]
                    }
                }
                dat$value <- dat[[chosen]]
                # drop the original climate columns (except the new 'value')
                to_keep <- c(
                    "Date",
                    "Region",
                    "variable",
                    "value",
                    setdiff(names(dat), c(climate_cols))
                )
                dat <- dat[to_keep]
            }

            # ensure mandatory columns exist
            if (!"Region" %in% names(dat)) {
                dat$Region <- NA
            }
            if (!"Date" %in% names(dat)) {
                dat$Date <- NA
            }
            dat
        })
        df_combined <- do.call(rbind, dfs_var)
        rownames(df_combined) <- NULL
        df_combined
    })
    names(read_list) <- names(files_named)

    if (combine) {
        # keep only data.frames (drop unexpected objects)
        dfs <- Filter(function(x) is.data.frame(x), read_list)
        if (length(dfs) == 0) {
            return(data.frame())
        }

        # ensure each df has the core columns and align order
        dfs_core <- lapply(dfs, function(df) {
            # ensure Date, Region, value, variable exist
            if (!"Date" %in% names(df)) {
                df$Date <- NA
            }
            if (!"Region" %in% names(df)) {
                df$Region <- NA
            }
            if (!"value" %in% names(df)) {
                # try to detect a single remaining climate-like column
                other_cols <- setdiff(
                    names(df),
                    c("Date", "Region", "variable")
                )
                if (length(other_cols) >= 1) {
                    df$value <- df[[other_cols[1]]]
                } else {
                    df$value <- NA
                }
            }
            if (!"variable" %in% names(df)) {
                df$variable <- NA
            }
            # keep only the core long-format columns
            df[, intersect(c("Date", "Region", "variable", "value"), names(df))]
        })

        # bind into a single long table
        long <- do.call(rbind, dfs_core)
        rownames(long) <- NULL
        # ensure Date is Date
        if ("Date" %in% names(long)) {
            long$Date <- as.Date(long$Date)
        }

        # pivot to wide: variables become column names
        if (requireNamespace("tidyr", quietly = TRUE)) {
            wide <- tidyr::pivot_wider(
                long,
                names_from = "variable",
                values_from = "value"
            )
        } else {
            # base R fallback using reshape
            wide <- reshape(
                long[, c("Date", "Region", "variable", "value")],
                idvar = c("Date", "Region"),
                timevar = "variable",
                direction = "wide"
            )
            # reshape creates column names like value.<var>; clean them
            names(wide) <- sub("^value\\.", "", names(wide))
        }
        # return wide data.frame (Date, Region, <vars...>)
        return(wide)
    }
    read_list
}

# For breakpoint calculations
calc_piecewise_rss <- function(series_vec, break_idx) {
    n <- length(series_vec)
    if (n == 0L) {
        return(NA_real_)
    }
    break_idx <- as.integer(break_idx)
    break_idx <- break_idx[is.finite(break_idx)]
    break_idx <- break_idx[break_idx >= 1L & break_idx < n]
    break_idx <- sort(unique(break_idx))

    seg_starts <- c(1L, break_idx + 1L)
    seg_ends <- c(break_idx, n)

    rss <- 0
    for (i in seq_along(seg_starts)) {
        seg <- series_vec[seg_starts[i]:seg_ends[i]]
        mu <- mean(seg, na.rm = TRUE)
        rss <- rss + sum((seg - mu)^2, na.rm = TRUE)
    }
    rss
}

find_optimal_break_idx <- function(series_vec, m, min_seg) {
    n <- length(series_vec)
    if (n == 0L || m < 0L) {
        return(list(break_idx = integer(0), rss = NA_real_))
    }
    if (m == 0L) {
        return(list(
            break_idx = integer(0),
            rss = calc_piecewise_rss(series_vec, integer(0))
        ))
    }

    candidate_idx <- seq.int(1L, n - 1L)
    if (length(candidate_idx) < m) {
        return(list(break_idx = integer(0), rss = NA_real_))
    }

    n_combos <- choose(length(candidate_idx), m)
    if (!is.finite(n_combos) || n_combos > 2e5) {
        return(list(break_idx = integer(0), rss = NA_real_))
    }

    combos <- utils::combn(candidate_idx, m)
    best_idx <- integer(0)
    best_rss <- Inf

    for (j in seq_len(ncol(combos))) {
        idx <- sort(as.integer(combos[, j]))
        seg_starts <- c(1L, idx + 1L)
        seg_ends <- c(idx, n)
        seg_lengths <- seg_ends - seg_starts + 1L
        if (any(seg_lengths < min_seg)) {
            next
        }
        rss <- calc_piecewise_rss(series_vec, idx)
        if (is.finite(rss) && rss < best_rss) {
            best_rss <- rss
            best_idx <- idx
        }
    }

    if (!is.finite(best_rss)) {
        return(list(break_idx = integer(0), rss = NA_real_))
    }
    list(break_idx = best_idx, rss = best_rss)
}


# Keep the raw age-specific residual pattern instead of
# subtracting the age mean before training. The model still gets the same
# climate and time features, but the target flow is rebuilt without the age
# bias correction used by CNN-LSTM and GNN-LSTM.
restore_age_bias_to_targets <- function(
    target_matrix,
    target_columns,
    age_bias_values
) {
    age_labels <- as.character(as.integer(sub("_.*", "", target_columns)))
    sweep(target_matrix, 2, age_bias_values[age_labels], "+")
}

make_recency_sample_weights <- function(week_vec, strength = 0.04) {
    year_offset <- lubridate::isoyear(week_vec) -
        min(lubridate::isoyear(week_vec), na.rm = TRUE)
    weights <- 1 + strength * year_offset
    weights / mean(weights)
}

get_age_metric_cols <- function(df, age_value) {
    names(df)[startsWith(names(df), paste0(age_value, "_"))]
}

get_region_metric_cols <- function(df, region_value) {
    names(df)[endsWith(names(df), paste0("_", region_value))]
}

resolve_year_index <- function(year_values, available_years, context) {
    year_values <- as.integer(year_values)
    available_years <- as.integer(available_years)
    year_idx <- match(year_values, available_years)
    if (anyNA(year_idx)) {
        missing_years <- unique(year_values[is.na(year_idx)])
        stop(
            sprintf(
                "Unable to resolve %s for year(s): %s",
                context,
                paste(missing_years, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    as.integer(year_idx)
}

# ------------------------------------------------------------------------------
# Function: combine_logmxt_and_residuals_mxt
# Description: Combines baseline mortality rates (derived from baseline log mxt)
# with predicted residuals in mxt space, then returns adjusted log mxt wide table.
# Inputs:
#   - baseline_logmxt_wide: baseline log mortality wide table
#   - residuals_wide: predicted residuals in mortality-rate space
#   - week_col: week/date column
#   - eps: lower bound to keep adjusted mxt strictly positive
# Output:
#   - adjusted wide table on log-mortality scale
# ------------------------------------------------------------------------------
combine_logmxt_and_residuals_mxt <- function(
    baseline_logmxt_wide,
    residuals_wide,
    week_col = "Week",
    eps = 1e-12
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

    baseline_mxt <- exp(b_mat)
    out_mxt <- pmax(baseline_mxt + r_mat, eps)
    out_log <- log(out_mxt)

    out <- as.data.frame(out_log)
    out[[week_col]] <- week_out
    out <- out[, c(week_col, common_cols), drop = FALSE]
    out
}

# Keep keras input ordering aligned with the model definition even when
# optional branches are disabled.
select_model_inputs <- function(model, candidate_inputs) {
    if (length(candidate_inputs) == 0) {
        return(list())
    }

    model_input_names <- vapply(model$inputs, function(i) i$name, character(1))
    missing_inputs <- setdiff(model_input_names, names(candidate_inputs))
    if (length(missing_inputs) > 0) {
        stop(
            "select_model_inputs: model expects inputs not found in candidates: ",
            paste(missing_inputs, collapse = ", "),
            "\nAvailable candidate keys: ",
            paste(names(candidate_inputs), collapse = ", ")
        )
    }

    selected_inputs <- candidate_inputs[model_input_names]
    null_inputs <- model_input_names[
        vapply(selected_inputs, is.null, logical(1))
    ]
    if (length(null_inputs) > 0) {
        stop(
            "select_model_inputs: candidate inputs were NULL for: ",
            paste(null_inputs, collapse = ", "),
            "\nAvailable candidate keys: ",
            paste(names(candidate_inputs), collapse = ", ")
        )
    }

    unname(selected_inputs)
}

make_peak_aware_sample_weights <- function(
    target_weeks,
    observed_residuals,
    recency_strength = 0.04,
    peak_strength = 2.0,
    peak_quantile = 0.90
) {
    w <- make_recency_sample_weights(target_weeks, strength = recency_strength)

    abs_res <- abs(observed_residuals)
    abs_res[!is.finite(abs_res)] <- 0
    threshold <- quantile(abs_res, peak_quantile, na.rm = TRUE)

    excess <- pmax(abs_res - threshold, 0) /
        (max(abs_res, na.rm = TRUE) - threshold + 1e-8)
    peak_multiplier <- 1.0 + (peak_strength - 1.0) * excess

    w <- w * peak_multiplier
    w / mean(w)
}


broadcast_region_residuals_to_baseline <- function(
    region_residuals_wide,
    baseline_wide,
    week_col = "Week"
) {
    if (
        !is.data.frame(region_residuals_wide) || !is.data.frame(baseline_wide)
    ) {
        stop("MortFCNet residual broadcasting requires data.frames")
    }

    region_residuals_wide[[week_col]] <- as.Date(region_residuals_wide[[
        week_col
    ]])
    baseline_wide[[week_col]] <- as.Date(baseline_wide[[week_col]])

    common_weeks <- sort(intersect(
        baseline_wide[[week_col]],
        region_residuals_wide[[week_col]]
    ))
    baseline_wide <- baseline_wide[
        baseline_wide[[week_col]] %in% common_weeks,
        ,
        drop = FALSE
    ]
    region_residuals_wide <- region_residuals_wide[
        region_residuals_wide[[week_col]] %in% common_weeks,
        ,
        drop = FALSE
    ]

    baseline_wide <- baseline_wide[
        order(baseline_wide[[week_col]]),
        ,
        drop = FALSE
    ]
    region_residuals_wide <- region_residuals_wide[
        order(region_residuals_wide[[week_col]]),
        ,
        drop = FALSE
    ]

    out <- data.frame(Week = baseline_wide[[week_col]], check.names = FALSE)
    baseline_cols <- setdiff(names(baseline_wide), week_col)

    for (col_name in baseline_cols) {
        region_name <- sub("^.*_", "", col_name)
        if (!region_name %in% names(region_residuals_wide)) {
            stop(
                "MortFCNet residual table is missing region column: ",
                region_name
            )
        }
        out[[col_name]] <- region_residuals_wide[[region_name]]
    }

    out
}
