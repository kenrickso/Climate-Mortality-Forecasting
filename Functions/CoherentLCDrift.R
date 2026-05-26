################################################################################
# COHERENT MULTI-POPULATION LEE-CARTER (Li & Lee–type)
#
# Provides a single exported function `forecast_kappa2_coherent()` used by
# `Final.R` to produce coherent forecasts of regional kappa2 (Li & Lee style).
################################################################################

#' Forecast kappa2 using the coherent Li & Lee decomposition.
#'
#' @param kappa2_matrix  Matrix (n_years x n_regions) of fitted kappa2 values.
#'                        Columns = regions, rows = years. Must have colnames.
#' @param n_ahead         Number of years to forecast.
#' @param weights         Optional vector of region weights (e.g. total
#'                        exposure). If NULL, uses equal weights (rowMeans).
#' @param verbose         Print diagnostics?
#' @return A list with:
#'   - kappa2_forecast: matrix (n_ahead x n_regions) of forecasted kappa2
#'   - kappa_fr:        numeric, common national factor (fitted)
#'   - kappa_fr_forecast: numeric, forecasted national factor
#'   - u_matrix:        matrix (n_years x n_regions) of fitted deviations
#'   - u_forecast:      matrix (n_ahead x n_regions) of forecasted deviations
#'   - phi:             named numeric vector of AR(1) coefficients per region
#'   - drift:           drift of the common national factor
#'   - arima_fr:        fitted ARIMA object for the common factor
#'   - ar_fits:         list of AR(1) fits per region
forecast_kappa2_coherent <- function(
    kappa2_matrix,
    n_ahead,
    weights = NULL,
    verbose = TRUE
) {
    stopifnot(is.matrix(kappa2_matrix), ncol(kappa2_matrix) >= 2)
    R <- ncol(kappa2_matrix)
    region_names <- colnames(kappa2_matrix)
    if (is.null(region_names)) {
        region_names <- paste0("R", seq_len(R))
        colnames(kappa2_matrix) <- region_names
    }

    # === Step 1: Estimate common national κ^(F) ===
    if (!is.null(weights)) {
        stopifnot(length(weights) == R)
        w <- weights / sum(weights)
        kappa_fr <- as.numeric(kappa2_matrix %*% w)
    } else {
        kappa_fr <- rowMeans(kappa2_matrix)
    }

    if (verbose) {
        cat(
            "  Coherent LC: national kappa range =",
            sprintf("[%.3f, %.3f]", min(kappa_fr), max(kappa_fr)),
            "\n"
        )
    }

    # === Step 2: Fit Random Walk with Drift to κ^(F) ===
    arima_fr <- tryCatch(
        {
            forecast::auto.arima(
                kappa_fr,
                seasonal = FALSE,
                max.p = 2,
                max.q = 2
            )
        },
        error = function(e) {
            # Fallback to simple RW with drift
            tryCatch(
                stats::arima(kappa_fr, order = c(0, 1, 0), include.mean = TRUE),
                error = function(e2) NULL
            )
        }
    )

    if (!is.null(arima_fr)) {
        fc_fr <- forecast::forecast(arima_fr, h = n_ahead)
        kappa_fr_forecast <- as.numeric(fc_fr$mean)
        drift <- ifelse(
            "drift" %in% names(coef(arima_fr)),
            coef(arima_fr)["drift"],
            mean(diff(kappa_fr), na.rm = TRUE)
        )
    } else {
        # Ultimate fallback: linear drift
        drift <- mean(diff(kappa_fr), na.rm = TRUE)
        last_val <- tail(kappa_fr, 1)
        kappa_fr_forecast <- last_val + (1:n_ahead) * drift
        warning(
            "Coherent LC: ARIMA failed for national factor; using linear drift."
        )
    }

    if (verbose) {
        cat("  Coherent LC: national drift =", sprintf("%.5f", drift), "\n")
    }

    # === Step 3: Compute regional residuals u_{t,r} = κ_{t,r} - κ_t^(F) ===
    u_matrix <- sweep(kappa2_matrix, 1, kappa_fr, "-")

    # === Step 4: Fit AR(1) per region ===
    phi <- numeric(R)
    names(phi) <- region_names
    ar_fits <- vector("list", R)
    names(ar_fits) <- region_names

    for (r in seq_len(R)) {
        u_r <- u_matrix[, r]
        u_r_finite <- u_r[is.finite(u_r)]

        if (length(u_r_finite) < 4) {
            # Too short for AR(1): assume quick mean reversion
            phi[r] <- 0
            ar_fits[[r]] <- NULL
            next
        }

        ar_fit <- tryCatch(
            {
                stats::arima(u_r, order = c(1, 0, 0), include.mean = TRUE)
            },
            error = function(e) NULL
        )

        if (!is.null(ar_fit)) {
            phi[r] <- coef(ar_fit)["ar1"]
            ar_fits[[r]] <- ar_fit
        } else {
            phi[r] <- 0
            ar_fits[[r]] <- NULL
        }
    }

    # Clamp |phi| < 1 for stationarity
    phi <- pmin(pmax(phi, -0.99), 0.99)

    if (verbose) {
        cat(
            "  Coherent LC: AR(1) phi range =",
            sprintf("[%.3f, %.3f]", min(phi), max(phi)),
            "\n"
        )
    }

    # === Step 5: Forecast regional deviations ===
    u_forecast <- matrix(NA, nrow = n_ahead, ncol = R)
    colnames(u_forecast) <- region_names

    for (r in seq_len(R)) {
        if (!is.null(ar_fits[[r]])) {
            u_forecast[, r] <- as.numeric(
                predict(ar_fits[[r]], n.ahead = n_ahead)$pred
            )
        } else {
            # No AR fit: deviations revert to zero immediately.
            u_forecast[, r] <- rep(0, n_ahead)
        }
    }

    # === Step 6: Reconstruct regional κ forecast ===
    kappa2_forecast <- sweep(u_forecast, 1, kappa_fr_forecast, "+")
    colnames(kappa2_forecast) <- region_names

    if (!is.null(rownames(kappa2_matrix))) {
        last_year <- max(as.numeric(rownames(kappa2_matrix)))
        rownames(kappa2_forecast) <- as.character(
            (last_year + 1):(last_year + n_ahead)
        )
    }

    if (verbose) {
        cat(
            "  Coherent LC: forecast complete.",
            R,
            "regions x",
            n_ahead,
            "years ahead\n"
        )
    }

    list(
        kappa2_forecast = kappa2_forecast,
        kappa_fr = kappa_fr,
        kappa_fr_forecast = kappa_fr_forecast,
        u_matrix = u_matrix,
        u_forecast = u_forecast,
        phi = phi,
        drift = drift,
        arima_fr = arima_fr,
        ar_fits = ar_fits
    )
}
