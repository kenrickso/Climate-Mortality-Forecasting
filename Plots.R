###############################################################################
# Lee Carter
###############################################################################
I <- solve(-lca$H.base)
sd1 <- sqrt(diag(I))
region_legend_guide <- guide_legend(
    title = "Region",
    nrow = 3,
    ncol = 7,
    byrow = TRUE
)

# beta1
se.beta1 <- sd1[grepl('Alpha', names(sd1), fixed = TRUE)]
df.beta1 <- data.frame(
    'age' = as.integer(sub(".*Age([0-9]+).*", "\\1", names(se.beta1))),
    'region' = sub(".*Region", "", names(se.beta1)),
    'value' = as.vector(lca$beta1),
    'lower' = as.vector(lca$beta1) - 1.96 * unname(se.beta1),
    'upper' = as.vector(lca$beta1) + 1.96 * unname(se.beta1)
)
p1 <- ggplot(df.beta1) +
    theme_bw(base_size = 15) +
    geom_line(aes(x = age, y = value, group = region, col = region)) +
    geom_ribbon(
        aes(x = age, ymin = lower, ymax = upper, group = region, fill = region),
        alpha = 0.1
    ) +
    xlab('Age') +
    ylab(bquote(alpha['x,r'])) +
    labs(color = "Region", fill = "Region") +
    theme(
        axis.title.y = element_text(size = 20),
        legend.position = 'bottom',
        legend.direction = "horizontal"
    ) +
    guides(color = region_legend_guide, fill = region_legend_guide)

# beta 2
se.beta2 <- c('Beta2_Age50' = 0, sd1[grepl('Beta2', names(sd1), fixed = TRUE)])
df.beta2 <- data.frame(
    'age' = as.integer(sub(".*Age([0-9]+).*", "\\1", names(se.beta2))),
    'value' = as.vector(lca$beta2),
    'lower' = as.vector(lca$beta2) - 1.96 * as.vector(se.beta2),
    'upper' = as.vector(lca$beta2) + 1.96 * as.vector(se.beta2)
)
p2 <- ggplot(df.beta2) +
    theme_bw(base_size = 15) +
    geom_line(aes(x = age, y = value), col = 'black') +
    geom_ribbon(
        aes(x = age, ymin = lower, ymax = upper),
        fill = 'black',
        alpha = 0.1
    ) +
    xlab('Age') +
    ylab(bquote(beta['x,r'])) +
    theme(axis.title.y = element_text(size = 20))

# kappa2
vec <- rep(0, R)
names(vec) <- paste0('Kappa2_Year1_Region', dimnames(dtxr)[[3]])
vec <- c(vec, sd1[grepl('Kappa2', names(sd1), fixed = TRUE)])
se.kappa2 <- unlist(lapply(dimnames(dtxr)[[3]], function(r) {
    vec[grepl(paste0('Region', r), names(vec), fixed = TRUE)]
}))
df.kappa2 <- data.frame(
    'year' = as.integer(sub(".*Year([0-9]+).*", "\\1", names(se.kappa2))) +
        min(DfM_train$ISOYear) -
        1,
    'region' = sub(".*Region", "", names(se.kappa2)),
    'value' = as.vector(lca$kappa2),
    'lower' = as.vector(lca$kappa2) - 1.96 * as.vector(se.kappa2),
    'upper' = as.vector(lca$kappa2) + 1.96 * as.vector(se.kappa2)
)
p3 <- ggplot(df.kappa2) +
    theme_bw(base_size = 15) +
    geom_line(aes(x = year, y = value, group = region, col = region)) +
    geom_ribbon(
        aes(
            x = year,
            ymin = lower,
            ymax = upper,
            group = region,
            fill = region
        ),
        alpha = 0.1
    ) +
    ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.01, 0.08))
    ) +
    xlab('Year') +
    ylab(bquote(kappa['t,r'])) +
    labs(color = "Region", fill = "Region") +
    theme(
        axis.title.y = element_text(size = 20),
        plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5),
        legend.direction = "horizontal"
    ) +
    guides(color = region_legend_guide, fill = region_legend_guide)

# beta 3
se.beta3 <- c('Beta3_Age50' = 0, sd1[grepl('Beta3', names(sd1), fixed = TRUE)])
df.beta3 <- data.frame(
    'age' = as.integer(sub(".*Age([0-9]+).*", "\\1", names(se.beta3))),
    'value' = as.vector(lca$beta3),
    'lower' = as.vector(lca$beta3) - 1.96 * as.vector(se.beta3),
    'upper' = as.vector(lca$beta3) + 1.96 * as.vector(se.beta3)
)
p4 <- ggplot(df.beta3) +
    theme_bw(base_size = 15) +
    geom_line(aes(x = age, y = value), col = 'black') +
    geom_ribbon(
        aes(x = age, ymin = lower, ymax = upper),
        fill = 'black',
        alpha = 0.1
    ) +
    xlab('Age') +
    ylab(bquote(gamma['x'])) +
    theme(axis.title.y = element_text(size = 20))

# kappa3
vec <- rep(0, R)
names(vec) <- paste0('Kappa3_Week1_Region', dimnames(dtxr)[[3]])
vec <- c(vec, sd1[grepl('Kappa3', names(sd1), fixed = TRUE)])
se.kappa3 <- unlist(lapply(dimnames(dtxr)[[3]], function(r) {
    vec[grepl(paste0('Region', r), names(vec), fixed = TRUE)]
}))
df.kappa3 <- data.frame(
    'week' = as.integer(sub(".*Week([0-9]+).*", "\\1", names(se.kappa3))),
    'region' = sub(".*Region", "", names(se.kappa3)),
    'value' = as.vector(lca$kappa3),
    'lower' = as.vector(lca$kappa3) - 1.96 * as.vector(se.kappa3),
    'upper' = as.vector(lca$kappa3) + 1.96 * as.vector(se.kappa3)
)

p5 <- ggplot(df.kappa3) +
    theme_bw(base_size = 15) +
    geom_line(aes(x = week, y = value, group = region, col = region)) +
    geom_ribbon(
        aes(
            x = week,
            ymin = lower,
            ymax = upper,
            group = region,
            fill = region
        ),
        alpha = 0.1
    ) +
    xlab('Week') +
    ylab(bquote(lambda['w,r'])) +
    labs(color = "Region", fill = "Region") +
    theme(
        axis.title.y = element_text(size = 20),
        legend.direction = "horizontal"
    ) +
    guides(color = region_legend_guide, fill = region_legend_guide)

# Baseline parameters
p12345 <- ggpubr::ggarrange(
    p1,
    ggpubr::ggarrange(p2, p3, ncol = 2, labels = c("B", "C"), legend = 'none'),
    ggpubr::ggarrange(p4, p5, ncol = 2, labels = c("D", "E"), legend = 'none'),
    nrow = 3,
    labels = "A",
    common.legend = TRUE,
    legend = 'bottom'
) +
    theme(legend.direction = "horizontal") +
    guides(color = region_legend_guide, fill = region_legend_guide)

ggsave(
    filename = file.path("Results", "baseline_p12345.png"),
    plot = p12345,
    width = 10,
    height = 10,
    units = "in",
    dpi = 500
)

################################################################################
# Plot: per-region test MSE comparison across all models
################################################################################
lee_carter_blue <- "#1f4e8c"
mortfcnet_green <- "#1b9e77"
cnn_lstm_orange <- "#d95f02"
gnn_lstm_violet <- "#7570b3"
observed_black <- "#000000"

models_for_line <- c("BASELINE", "MortFCNet", "CNN-LSTM", "GNN-LSTM")

model_plot_colors <- c(
    "BASELINE" = lee_carter_blue,
    "MortFCNet" = mortfcnet_green,
    "CNN-LSTM" = cnn_lstm_orange,
    "GNN-LSTM" = gnn_lstm_violet
)

winner_map_colors <- c(
    "BASELINE" = lee_carter_blue,
    "MortFCNet" = mortfcnet_green,
    "CNN-LSTM" = cnn_lstm_orange,
    "GNN-LSTM" = gnn_lstm_violet,
    "TIE" = "#000000"
)

top_region_plot_colors <- c(
    "Observed" = observed_black,
    "LC" = lee_carter_blue,
    "MortFCNet" = mortfcnet_green,
    "CNN-LSTM" = cnn_lstm_orange,
    "GNN-LSTM" = gnn_lstm_violet
)

observed_linewidth <- 1.1
model_linewidth <- 1.0

base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
        plot.subtitle = ggplot2::element_text(
            hjust = 0.5,
            size = 11,
            color = "grey40"
        ),
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 11),
        axis.title = ggplot2::element_text(size = 13)
    )

scientific_x10_labels <- function(digits = 2) {
    function(x) {
        labs <- scales::label_scientific(digits = digits)(x)
        labs <- gsub("e\\+?", " %*% 10^", labs)
        labs <- gsub("^0 %\\*% 10\\^[-0-9]+$", "0", labs)
        labs <- gsub("10\\^0$", "1", labs)
        parse(text = labs)
    }
}

# Keep consistent region ordering (by baseline test MSE for readability)
region_order <- per_region_metrics_table %>%
    dplyr::filter(model == "BASELINE") %>%
    dplyr::arrange(mse_mxt_test) %>%
    dplyr::pull(region)

line_df <- per_region_metrics_table %>%
    dplyr::filter(as.character(model) %in% models_for_line) %>%
    dplyr::mutate(
        model = factor(as.character(model), levels = models_for_line),
        region = factor(region, levels = region_order)
    )

baseline_line_df <- line_df %>%
    dplyr::filter(model == "BASELINE")
other_models_point_df <- line_df %>%
    dplyr::filter(model != "BASELINE")

p_region_mse_compare <- ggplot2::ggplot(
    line_df,
    ggplot2::aes(
        x = region,
        y = mse_mxt_test,
        color = model,
        shape = model,
        group = model
    )
) +
    ggplot2::geom_line(
        data = baseline_line_df,
        linewidth = 1.0,
        lineend = "round"
    ) +
    ggplot2::geom_point(
        data = other_models_point_df,
        size = 2.5,
        alpha = 0.95,
        position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_color_manual(name = "Model", values = model_plot_colors) +
    ggplot2::scale_shape_manual(
        name = "Model",
        values = c(
            "BASELINE" = 1,
            "MortFCNet" = 15,
            "CNN-LSTM" = 16,
            "GNN-LSTM" = 17
        )
    ) +
    ggplot2::scale_y_continuous(
        labels = scientific_x10_labels(digits = 2)
    ) +
    ggplot2::labs(
        title = "Per-region Test MSE",
        x = "Region",
        y = "MSE",
        color = "Model"
    ) +
    base_theme +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(
            size = 11,
            angle = 90,
            vjust = 0.5,
            hjust = 1
        ),
        axis.text.y = ggplot2::element_text(size = 11),
        axis.title.x = ggplot2::element_text(size = 15),
        axis.title.y = ggplot2::element_text(size = 15),
        plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
        legend.title = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 12),
        legend.position = "top"
    )
ggplot2::guides(
    color = ggplot2::guide_legend(
        override.aes = list(
            shape = c(1, 15, 16, 17),
            linetype = c(1, 0, 0, 0)
        )
    )
)

plot_dir <- "Results/plots"
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(
    filename = file.path(plot_dir, "per_region_mse_test_line.png"),
    plot = p_region_mse_compare,
    width = 12,
    height = 6,
    dpi = 300
)


# Format and print per-age metrics (mortality rate) with 6 decimal places for numeric columns (show all rows)
# Determine per-region winners by aggregated MSE (mse_mxt_test)
agg_winners_mxt <- per_region_metrics_table %>%
    dplyr::filter(!is.na(model), !is.na(mse_mxt_test)) %>%
    dplyr::group_by(region) %>%
    dplyr::filter(mse_mxt_test == min(mse_mxt_test, na.rm = TRUE)) %>%
    dplyr::summarise(
        winners = paste(sort(as.character(model)), collapse = "/"),
        n_winners = dplyr::n(),
        .groups = "drop"
    )

# Helper: build summary table of region wins and fractional wins
make_summary <- function(agg) {
    sole_counts <- agg %>%
        dplyr::filter(n_winners == 1) %>%
        dplyr::count(winners) %>%
        dplyr::rename(model = winners, sole_regions = n)

    fractional <- agg %>%
        dplyr::mutate(frac = 1 / n_winners) %>%
        tidyr::separate_rows(winners, sep = "/") %>%
        dplyr::group_by(winners) %>%
        dplyr::summarise(frac_regions = sum(frac), .groups = "drop") %>%
        dplyr::rename(model = winners)

    n_regions_total <- length(unique(per_region_metrics_table$region))

    dplyr::full_join(sole_counts, fractional, by = "model") %>%
        tidyr::replace_na(list(sole_regions = 0, frac_regions = 0)) %>%
        dplyr::arrange(dplyr::desc(frac_regions)) %>%
        dplyr::mutate(
            prop_sole = paste0(sole_regions, "/", n_regions_total),
            prop_frac = paste0(round(frac_regions, 2), "/", n_regions_total)
        )
}

# Produce and print region-level summary using RAW mortality MSE (mxt)
summary_regions_agg_mxt <- make_summary(agg_winners_mxt)

# Create and save a map showing which model is best per NUTS2 region (by mse_mxt_test)
region_winners_map <- agg_winners_mxt %>%
    dplyr::mutate(
        winner_label = ifelse(n_winners == 1, winners, "TIE"),
        winner = winners
    ) %>%
    dplyr::select(region, winner, winner_label)

# Join winners to shapefile and ensure factor ordering
shapef_map <- shapef %>%
    dplyr::left_join(region_winners_map, by = c("NUTS_ID" = "region"))

model_levels <- c("BASELINE", "MortFCNet", "CNN-LSTM", "GNN-LSTM", "TIE")
shapef_map$winner_label <- factor(
    shapef_map$winner_label,
    levels = model_levels
)

# Plot: per-region best model map
color_map <- winner_map_colors

p_region_winner_map <- ggplot2::ggplot(shapef_map) +
    ggplot2::geom_sf(
        ggplot2::aes(fill = winner_label),
        color = "grey40",
        size = 0.1
    ) +
    ggplot2::scale_fill_manual(
        values = color_map,
        na.value = "white",
        name = "Best model\n(by mse_mxt_test)"
    ) +
    ggplot2::labs(title = "Per-region Best Model (mse_mxt_test)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5)
    )

# Create heatmaps (choropleths) of mse_mxt_test per region for each model (NUTS2)
per_region_melt <- per_region_metrics_table %>%
    dplyr::filter(model %in% models_for_line) %>%
    dplyr::select(region, model, mse_mxt_test)

# Add a small epsilon and compute log10 to make color scale interpretable
eps <- 1e-12
per_region_melt <- per_region_melt %>%
    dplyr::mutate(
        mse_clamped = dplyr::if_else(
            is.na(mse_mxt_test),
            NA_real_,
            pmax(mse_mxt_test, eps)
        ),
        mse_log10 = log10(mse_clamped)
    )

# Compute robust (2% - 98%) bounds on the log10 scale and report them
lwr_log <- as.numeric(stats::quantile(
    per_region_melt$mse_log10,
    probs = 0.02,
    na.rm = TRUE
))
upr_log <- as.numeric(stats::quantile(
    per_region_melt$mse_log10,
    probs = 0.98,
    na.rm = TRUE
))

# Join to shapefile (will create repeated geometries for facetting)
shapef_heat <- shapef %>%
    dplyr::left_join(per_region_melt, by = c("NUTS_ID" = "region"))

# Ensure model is factor ordered
shapef_heat$model <- factor(shapef_heat$model, levels = models_for_line)

# Define colorbar breaks in log10 space (use rounded values) and show original MSE as labels
breaks_log <- seq(floor(lwr_log), ceiling(upr_log), by = 0.5)
labels_breaks <- scales::label_number(accuracy = 1e-6)(10^breaks_log)

# Compute percent improvement per region and model relative to BASELINE
baseline_reg <- per_region_metrics_table %>%
    dplyr::filter(model == "BASELINE") %>%
    dplyr::select(region, baseline_mse = mse_mxt_test)

per_region_pct <- per_region_metrics_table %>%
    dplyr::left_join(baseline_reg, by = "region") %>%
    dplyr::mutate(
        pct_improve = ifelse(
            is.na(baseline_mse) | baseline_mse == 0,
            NA_real_,
            (baseline_mse - mse_mxt_test) / baseline_mse * 100
        )
    )

# Join to shapefile for plotting
shp_pct <- shapef %>%
    dplyr::left_join(per_region_pct, by = c("NUTS_ID" = "region"))
shp_pct$model <- factor(shp_pct$model, levels = models_for_line)
shp_pct_faceted <- shp_pct %>%
    dplyr::filter(as.character(model) != "BASELINE")

pct_vals <- shp_pct_faceted$pct_improve[!is.na(shp_pct_faceted$pct_improve)]
pct_min <- floor(min(pct_vals))
pct_max <- ceiling(max(pct_vals))
pct_mid <- round((pct_min + pct_max) / 2, 1)

p_pct_faceted <- ggplot2::ggplot(shp_pct_faceted) +
    ggplot2::geom_sf(
        ggplot2::aes(fill = pct_improve),
        color = "grey40",
        size = 0.1
    ) +
    ggplot2::facet_wrap(~model, ncol = 3) +
    ggplot2::scale_fill_gradient2(
        low = "#d73027",
        mid = "white",
        high = "#1a9850",
        midpoint = 0,
        limits = c(pct_min, pct_max),
        breaks = c(pct_min, pct_mid, pct_max),
        labels = paste0(c(pct_min, pct_mid, pct_max), " %"),
        name = "Improvement",
        guide = ggplot2::guide_colorbar(barheight = grid::unit(6, "cm"))
    ) +
    ggplot2::labs(
        title = "Percent Improvement Against Baseline"
    ) +
    base_theme +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right",
        strip.text = ggplot2::element_text(size = 13, face = "bold"),
        panel.grid = ggplot2::element_blank()
    )

ggplot2::ggsave(
    filename = file.path(plot_dir, "facet_improvement.png"),
    plot = p_pct_faceted,
    width = 14,
    height = 8,
    dpi = 300
)

# Plot: embedding cluster diagnostics
################################################################################
build_region_embedding_pca_df <- function(
    model_obj,
    layer_name,
    region_names,
    model_label,
    k_clusters = 4L
) {
    region_names <- unique(as.character(region_names))
    n_regions <- length(region_names)
    if (n_regions == 0) {
        return(NULL)
    }

    emb_layer <- model_obj$get_layer(name = layer_name)
    emb_weights <- emb_layer$get_weights()[[1]]

    emb_region <- emb_weights[2:(n_regions + 1L), , drop = FALSE]
    rownames(emb_region) <- region_names

    pca_obj <- stats::prcomp(emb_region, center = TRUE, scale. = TRUE)
    pcs <- pca_obj$x[, 1:2, drop = FALSE]

    set.seed(2026)
    km <- stats::kmeans(pcs, centers = as.integer(k_clusters), nstart = 50)

    out <- data.frame(
        region = rownames(pcs),
        comp1 = pcs[, 1],
        comp2 = pcs[, 2],
        cluster = factor(km$cluster),
        model = model_label,
        stringsAsFactors = FALSE
    )

    centroids <- as.data.frame(km$centers)
    names(centroids) <- c("comp1", "comp2")
    centroids$cluster <- factor(seq_len(nrow(centroids)))

    out <- out %>%
        dplyr::left_join(
            centroids,
            by = "cluster",
            suffix = c("", "_centroid")
        )
    out
}

embedding_plot_parts <- list()

cnn_embedding_df <- build_region_embedding_pca_df(
    model_obj = mdl_cnn,
    layer_name = "region_embedding",
    region_names = region_levels,
    model_label = "CNN-LSTM",
    k_clusters = 4L
)
if (!is.null(cnn_embedding_df) && nrow(cnn_embedding_df) > 0) {
    embedding_plot_parts[[
        length(embedding_plot_parts) + 1L
    ]] <- cnn_embedding_df
}

gnn_embedding_df <- build_region_embedding_pca_df(
    model_obj = mdl_gnn,
    layer_name = "gnn_scalar_region_embedding",
    region_names = unique(seqs_gnn$Region),
    model_label = "GNN-LSTM",
    k_clusters = 4L
)
if (!is.null(gnn_embedding_df) && nrow(gnn_embedding_df) > 0) {
    embedding_plot_parts[[
        length(embedding_plot_parts) + 1L
    ]] <- gnn_embedding_df
}

embedding_plot_df <- dplyr::bind_rows(embedding_plot_parts)

embedding_results_plots_dir <- if (exists("results_plots_dir")) {
    results_plots_dir
} else {
    file.path("Results", "plots")
}
if (!dir.exists(embedding_results_plots_dir)) {
    dir.create(
        embedding_results_plots_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
}

# Additional: age 90+ plots for FRK2 and FRF3 at the end of the script
age_focus_plot_dir <- file.path(results_plots_dir, "age_focus")
dir.create(age_focus_plot_dir, recursive = TRUE, showWarnings = FALSE)

if (exists("observed_test") && exists("baseline_logmxt_wide")) {
    age_val <- 90L
    regions_to_plot <- c("FRK2", "FRF3", "FR10", "FRE1", "FRG0", "FRL0")

    get_series_local <- function(wide_df, age_value, region_value) {
        if (
            is.null(wide_df) ||
                !is.data.frame(wide_df) ||
                !("Week" %in% names(wide_df))
        ) {
            return(NULL)
        }
        candidates <- c(
            paste0(age_value, "_", region_value),
            paste0(region_value, "_", age_value)
        )
        col_found <- intersect(candidates, names(wide_df))[1]
        if (is.na(col_found) || !(col_found %in% names(wide_df))) {
            return(NULL)
        }
        s <- wide_df[, c("Week", col_found), drop = FALSE]
        names(s)[2] <- "logmxt"
        s$Week <- as.Date(s$Week)
        s$logmxt <- as.numeric(s$logmxt)
        s <- s[is.finite(s$logmxt), , drop = FALSE]
        if (nrow(s) == 0) {
            return(NULL)
        }
        s[order(s$Week), , drop = FALSE]
    }

    mort_wide <- cnn_wide <- gnn_wide <- NULL
    if (exists("resolve_model_wides", mode = "function")) {
        mort_wide <- tryCatch(
            resolve_model_wides("MortFCNet")$test,
            error = function(e) NULL
        )
        cnn_wide <- tryCatch(
            resolve_model_wides("CNN-LSTM")$test,
            error = function(e) NULL
        )
        gnn_wide <- tryCatch(
            resolve_model_wides("GNN-LSTM")$test,
            error = function(e) NULL
        )
    }

    all_plot_dfs <- list()

    for (region_value in regions_to_plot) {
        parts <- list()
        obs <- get_series_local(observed_test, age_val, region_value)
        lc <- get_series_local(baseline_logmxt_wide, age_val, region_value)
        mort <- if (!is.null(mort_wide)) {
            get_series_local(mort_wide, age_val, region_value)
        } else {
            NULL
        }
        cnn <- if (!is.null(cnn_wide)) {
            get_series_local(cnn_wide, age_val, region_value)
        } else {
            NULL
        }
        gnn <- if (!is.null(gnn_wide)) {
            get_series_local(gnn_wide, age_val, region_value)
        } else {
            NULL
        }

        if (!is.null(obs)) {
            parts[[length(parts) + 1]] <- dplyr::mutate(obs, model = "Observed")
        }
        if (!is.null(lc)) {
            parts[[length(parts) + 1]] <- dplyr::mutate(lc, model = "LC")
        }
        if (!is.null(mort)) {
            parts[[length(parts) + 1]] <- dplyr::mutate(
                mort,
                model = "MortFCNet"
            )
        }
        if (!is.null(cnn)) {
            parts[[length(parts) + 1]] <- dplyr::mutate(cnn, model = "CNN-LSTM")
        }
        if (!is.null(gnn)) {
            parts[[length(parts) + 1]] <- dplyr::mutate(gnn, model = "GNN-LSTM")
        }

        if (length(parts) == 0) {
            next
        }

        df <- dplyr::bind_rows(parts) %>%
            dplyr::mutate(
                rate = exp(logmxt),
                model = factor(
                    model,
                    levels = c(
                        "Observed",
                        "LC",
                        "MortFCNet",
                        "CNN-LSTM",
                        "GNN-LSTM"
                    )
                )
            )

        all_plot_dfs[[region_value]] <- df
    }

    all_rates <- unlist(lapply(all_plot_dfs, function(d) d$rate))
    shared_ylim <- range(all_rates, na.rm = TRUE)
    y_margin <- diff(shared_ylim) * 0.03
    shared_ylim <- shared_ylim + c(-y_margin, y_margin)

    # Colour scale: Observed is fixed gray; CNN/GNN pull from top_region_plot_colors.
    # Keeping Observed inside the scale (not hardcoded) is what allows ggplot2
    # to merge the colour + linetype legends into one.
    age90_colors <- c(
        "Observed" = "#999999",
        "LC" = top_region_plot_colors[["LC"]],
        "MortFCNet" = top_region_plot_colors[["MortFCNet"]],
        "CNN-LSTM" = top_region_plot_colors[["CNN-LSTM"]],
        "GNN-LSTM" = top_region_plot_colors[["GNN-LSTM"]]
    )

    for (region_value in names(all_plot_dfs)) {
        plot_df <- all_plot_dfs[[region_value]]

        p_region_age90 <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = Week, y = rate, colour = model, linetype = model)
        ) +
            # Layer 1 — Observed: thin, gray, semi-transparent
            ggplot2::geom_line(
                data = dplyr::filter(plot_df, model == "Observed"),
                linewidth = observed_linewidth * 0.55,
                alpha = 0.50
            ) +
            # Layer 2 — LC & MortFCNet: context only, de-emphasised
            ggplot2::geom_line(
                data = dplyr::filter(plot_df, model %in% c("LC", "MortFCNet")),
                linewidth = model_linewidth * 0.80,
                alpha = 0.60
            ) +
            # Layer 3 — CNN-LSTM & GNN-LSTM: focus models, drawn on top
            ggplot2::geom_line(
                data = dplyr::filter(
                    plot_df,
                    model %in% c("CNN-LSTM", "GNN-LSTM")
                ),
                linewidth = model_linewidth * 0.95, # pulled back so LC/MortFCNet show through
                alpha = 1.0
            ) +
            ggplot2::scale_colour_manual(
                name = "Model",
                values = age90_colors
            ) +
            ggplot2::scale_linetype_manual(
                name = "Model", # same name → single merged legend
                values = c(
                    "Observed" = "solid",
                    "LC" = "dashed",
                    "MortFCNet" = "dashed",
                    "CNN-LSTM" = "solid",
                    "GNN-LSTM" = "solid"
                )
            ) +
            ggplot2::scale_y_continuous(limits = shared_ylim) +
            ggplot2::labs(
                title = region_value,
                subtitle = "Age group: 90+",
                x = "Year",
                y = "Mortality Rate"
            ) +
            base_theme +
            ggplot2::theme(legend.position = "top")

        safe_region_name <- tolower(gsub("[^A-Za-z0-9]+", "_", region_value))
        safe_region_name <- gsub("(^_+|_+$)", "", safe_region_name)

        ggplot2::ggsave(
            filename = file.path(
                age_focus_plot_dir,
                paste0("region_", safe_region_name, "_age90plus.png")
            ),
            plot = p_region_age90,
            width = 12,
            height = 6,
            dpi = 300
        )
    }
}

if (nrow(embedding_plot_df) > 0) {
    cluster_colors <- c(
        "1" = "#F8766D",
        "2" = "#7CAE00",
        "3" = "#00BFC4",
        "4" = "#C77CFF"
    )

    make_embedding_cluster_plot <- function(df_model, model_name) {
        y_rng <- range(df_model$comp2, na.rm = TRUE)
        y_span <- diff(y_rng)
        if (!is.finite(y_span) || y_span <= 0) {
            y_span <- 1
        }
        nudge_val <- 0.03 * y_span

        p <- ggplot2::ggplot(
            df_model,
            ggplot2::aes(x = comp1, y = comp2)
        ) +
            ggplot2::geom_segment(
                ggplot2::aes(xend = comp1_centroid, yend = comp2_centroid),
                color = "grey70",
                alpha = 0.45,
                linewidth = 0.3
            ) +
            ggplot2::geom_point(
                ggplot2::aes(fill = cluster),
                shape = 21,
                size = 3.8,
                stroke = 0.35,
                color = "grey20",
                alpha = 0.95
            ) +
            ggplot2::scale_fill_manual(
                values = cluster_colors,
                name = "Cluster"
            ) +
            ggplot2::coord_fixed(ratio = 1) +
            ggplot2::labs(
                title = paste0(
                    model_name,
                    " Region Embedding PCA Clusters"
                ),
                x = "PC1",
                y = "PC2"
            ) +
            ggplot2::theme_bw(base_size = 13) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 17, hjust = 0.5),
                legend.position = "right",
                plot.margin = ggplot2::margin(10, 10, 10, 10),
                aspect.ratio = 1
            )

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p +
                ggrepel::geom_text_repel(
                    ggplot2::aes(label = region, color = cluster),
                    size = 3.6,
                    min.segment.length = 0,
                    seed = 2026,
                    box.padding = 0.25,
                    point.padding = 0.2,
                    max.overlaps = Inf,
                    show.legend = FALSE
                ) +
                ggplot2::scale_color_manual(
                    values = cluster_colors,
                    guide = "none"
                )
        } else {
            p <- p +
                ggplot2::geom_label(
                    ggplot2::aes(label = region),
                    size = 3.2,
                    nudge_y = nudge_val,
                    label.padding = grid::unit(0.1, "lines"),
                    label.r = grid::unit(0.08, "lines"),
                    label.size = 0.15,
                    fill = scales::alpha("white", 0.82),
                    color = "black",
                    fontface = "bold",
                    show.legend = FALSE
                )
        }

        p
    }

    make_embedding_cluster_map <- function(df_model, model_name) {
        if (!exists("shapef") || is.null(shapef) || nrow(shapef) == 0) {
            return(NULL)
        }

        map_df <- shapef %>%
            dplyr::left_join(
                df_model %>%
                    dplyr::filter(model == model_name) %>%
                    dplyr::select(region, cluster),
                by = c("NUTS_ID" = "region")
            )

        ggplot2::ggplot(map_df) +
            ggplot2::geom_sf(
                ggplot2::aes(fill = cluster),
                color = "grey40",
                size = 0.1
            ) +
            ggplot2::scale_fill_manual(
                values = cluster_colors,
                drop = FALSE,
                na.value = "white",
                name = "Cluster"
            ) +
            ggplot2::coord_sf(expand = FALSE) +
            ggplot2::labs(
                title = paste0(model_name, " Region Embedding PCA Cluster Map")
            ) +
            ggplot2::theme_void(base_size = 13) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.title = ggplot2::element_text(size = 12),
                legend.text = ggplot2::element_text(size = 11),
                plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
            )
    }

    models_to_plot <- unique(as.character(embedding_plot_df$model))
    for (mdl_name in models_to_plot) {
        df_m <- embedding_plot_df %>% dplyr::filter(model == mdl_name)
        if (nrow(df_m) == 0) {
            next
        }

        p_region_embedding_clusters <- make_embedding_cluster_plot(
            df_m,
            mdl_name
        )
        fn_model <- tolower(gsub("[^A-Za-z0-9]+", "_", mdl_name))
        fn_model <- gsub("(^_+|_+$)", "", fn_model)
        ggplot2::ggsave(
            filename = file.path(
                embedding_results_plots_dir,
                paste0("regionClusters_", fn_model, ".png")
            ),
            plot = p_region_embedding_clusters,
            width = 9,
            height = 6,
            dpi = 300
        )

        if (exists("shapef")) {
            p_region_embedding_map <- make_embedding_cluster_map(
                embedding_plot_df,
                mdl_name
            )

            if (!is.null(p_region_embedding_map)) {
                p_side_by_side <- ggpubr::ggarrange(
                    p_region_embedding_clusters,
                    p_region_embedding_map,
                    ncol = 2,
                    labels = c("A", "B"),
                    widths = c(1, 1)
                )

                ggplot2::ggsave(
                    filename = file.path(
                        embedding_results_plots_dir,
                        paste0("regionClusters_side_by_side_", fn_model, ".png")
                    ),
                    plot = p_side_by_side,
                    width = 16,
                    height = 7,
                    dpi = 300
                )
            }
        }
    }

    embedding_cluster_map_df <- embedding_plot_df %>%
        dplyr::filter(model %in% c("CNN-LSTM", "GNN-LSTM")) %>%
        dplyr::mutate(
            model = factor(model, levels = c("CNN-LSTM", "GNN-LSTM")),
            cluster = factor(cluster, levels = names(cluster_colors))
        ) %>%
        dplyr::distinct(region, model, cluster)

    if (exists("shapef") && nrow(embedding_cluster_map_df) > 0) {
        embedding_cluster_map_sf <- shapef %>%
            dplyr::left_join(
                embedding_cluster_map_df,
                by = c("NUTS_ID" = "region")
            )

        p_embedding_cluster_map <- ggplot2::ggplot(embedding_cluster_map_sf) +
            ggplot2::geom_sf(
                ggplot2::aes(fill = cluster),
                color = "grey40",
                size = 0.1
            ) +
            ggplot2::facet_wrap(~model, ncol = 2) +
            ggplot2::scale_fill_manual(
                values = cluster_colors,
                drop = FALSE,
                na.value = "white",
                name = "Cluster"
            ) +
            ggplot2::coord_sf(expand = FALSE) +
            ggplot2::labs(
                title = "Region Embedding PCA Cluster Choropleths",
                subtitle = paste0(
                    "CNN-LSTM and GNN-LSTM region embeddings clustered in PCA space"
                )
            ) +
            ggplot2::theme_void(base_size = 13) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
                plot.subtitle = ggplot2::element_text(
                    hjust = 0.5,
                    size = 11,
                    color = "grey40"
                ),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.title = ggplot2::element_text(size = 12),
                legend.text = ggplot2::element_text(size = 11),
                strip.text = ggplot2::element_text(size = 13, face = "bold"),
                plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
            )

        ggplot2::ggsave(
            filename = file.path(
                embedding_results_plots_dir,
                "region_embedding_cluster_choropleth.png"
            ),
            plot = p_embedding_cluster_map,
            width = 12,
            height = 8,
            dpi = 300
        )
    }
}

################################################################################
# Numeric formatting helpers for table exports
################################################################################
format_scientific <- function(x) {
    x <- as.numeric(x)
    ifelse(
        is.na(x),
        NA_character_,
        vapply(
            x,
            function(xx) {
                s <- sprintf("%.4e", xx)
                parts <- regmatches(
                    s,
                    regexec("^([0-9.]+)e([+-]?)(0*)([0-9]+)$", s)
                )[[1]]
                if (length(parts) != 5) {
                    return(s)
                }
                mant <- parts[2]
                sign <- parts[3]
                expo <- parts[5]
                if (sign == "+") {
                    sign <- ""
                }
                paste0(mant, "$\\times10^{", sign, expo, "}$")
            },
            character(1L)
        )
    )
}

################################################################################
# Per-age improvement table (BASELINE vs CNN-LSTM vs GNN-LSTM)
################################################################################
baseline_df <- per_age_metrics_table %>%
    dplyr::filter(model == "BASELINE") %>%
    dplyr::select(age_group = age, BASELINE = baseline_mse_mxt_test)

mortfcnet_df <- per_age_metrics_table %>%
    dplyr::filter(model == "MortFCNet") %>%
    dplyr::select(
        age_group = age,
        MortFCNet = mse_mxt_test
    )

cnn_df <- per_age_metrics_table %>%
    dplyr::filter(model == "CNN-LSTM") %>%
    dplyr::select(
        age_group = age,
        `CNN-LSTM` = mse_mxt_test
    )

gnn_df <- per_age_metrics_table %>%
    dplyr::filter(model == "GNN-LSTM") %>%
    dplyr::select(
        age_group = age,
        `GNN-LSTM` = mse_mxt_test
    )

per_age_compare <- baseline_df %>%
    dplyr::left_join(mortfcnet_df, by = "age_group") %>%
    dplyr::left_join(cnn_df, by = "age_group") %>%
    dplyr::left_join(gnn_df, by = "age_group") %>%
    dplyr::arrange(age_group)

per_age_compare_formatted <- per_age_compare %>%
    dplyr::mutate(
        dplyr::across(-age_group, format_scientific)
    )

dir.create("Results/tables", recursive = TRUE, showWarnings = FALSE)
if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(
        per_age_compare_formatted,
        file.path("Results/tables", "per_age_gnn_improvement.csv")
    )
} else {
    write.csv(
        per_age_compare_formatted,
        file = file.path("Results/tables", "per_age_gnn_improvement.csv"),
        row.names = FALSE
    )
}

################################################################################
# Tables and plots: train/test MSE summaries and region-level comparisons
################################################################################

results_tables_dir <- "Results/tables"
dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Per-model train/test loss table
model_mse_summary <- metrics_table %>%
    dplyr::select(model, mse_mxt_train, mse_mxt_test) %>%
    dplyr::arrange(model)

if (nrow(model_mse_summary) > 0) {
    if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_csv(
            model_mse_summary,
            file.path(results_tables_dir, "model_train_test_mse_summary.csv")
        )
    } else {
        write.csv(
            model_mse_summary,
            file.path(results_tables_dir, "model_train_test_mse_summary.csv"),
            row.names = FALSE
        )
    }
}

# Per-age train/test loss table
age_col <- if ("age_group" %in% names(per_age_metrics_table)) {
    "age_group"
} else {
    "age"
}

if (age_col == "age_group") {
    per_age_mse_summary <- per_age_metrics_table %>%
        dplyr::select(
            age = age_group,
            model,
            mse_mxt_train,
            mse_mxt_test
        ) %>%
        dplyr::arrange(age, model)
} else {
    per_age_mse_summary <- per_age_metrics_table %>%
        dplyr::select(age, model, mse_mxt_train, mse_mxt_test) %>%
        dplyr::arrange(age, model)
}

if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(
        per_age_mse_summary,
        file.path(results_tables_dir, "age_train_test_mse_summary.csv")
    )
} else {
    write.csv(
        per_age_mse_summary,
        file.path(results_tables_dir, "age_train_test_mse_summary.csv"),
        row.names = FALSE
    )
}

# Per-region train/test loss table
per_region_mse_summary <- per_region_metrics_table %>%
    dplyr::select(region, model, mse_mxt_train, mse_mxt_test) %>%
    dplyr::arrange(region, model)

if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(
        per_region_mse_summary,
        file.path(results_tables_dir, "region_train_test_mse_summary.csv")
    )
} else {
    write.csv(
        per_region_mse_summary,
        file.path(results_tables_dir, "region_train_test_mse_summary.csv"),
        row.names = FALSE
    )
}

################################################################################
# Plot: strongest region improvements against baseline
################################################################################

if (exists("per_region_metrics_table")) {
    baseline_mse <- per_region_metrics_table %>%
        dplyr::filter(model == "BASELINE") %>%
        dplyr::select(region, baseline_mse = mse_mxt_test)

    top_region_improvement <- per_region_metrics_table %>%
        dplyr::filter(model %in% c("CNN-LSTM", "GNN-LSTM")) %>%
        dplyr::left_join(baseline_mse, by = "region") %>%
        dplyr::mutate(
            improvement_pct = dplyr::if_else(
                is.na(baseline_mse) | baseline_mse == 0,
                NA_real_,
                (baseline_mse - mse_mxt_test) / baseline_mse * 100
            )
        ) %>%
        dplyr::filter(!is.na(improvement_pct)) %>%
        dplyr::arrange(dplyr::desc(improvement_pct)) %>%
        dplyr::slice_head(n = 3) %>%
        dplyr::mutate(best_model = as.character(model)) %>%
        dplyr::select(
            region,
            best_model,
            baseline_mse,
            model_mse = mse_mxt_test,
            improvement_pct
        )

    if (nrow(top_region_improvement) > 0) {
        if (requireNamespace("readr", quietly = TRUE)) {
            readr::write_csv(
                top_region_improvement,
                file.path(
                    results_tables_dir,
                    "top3_region_improvement_summary.csv"
                )
            )
        } else {
            write.csv(
                top_region_improvement,
                file.path(
                    results_tables_dir,
                    "top3_region_improvement_summary.csv"
                ),
                row.names = FALSE
            )
        }
        observed_all <- NULL
        if (exists("observed_test")) {
            observed_all <- observed_test %>%
                dplyr::mutate(Week = as.Date(Week))
        }

        if (!is.null(observed_all)) {
            test_end <- max(observed_all$Week, na.rm = TRUE)
            test_start_cutoff <- test_end - 365 * 5

            get_mxt_series_test <- function(wide_df, region) {
                series <- get_mxt_series(wide_df, region)
                if (is.null(series)) {
                    return(NULL)
                }
                series <- series %>%
                    dplyr::mutate(Week = as.Date(Week)) %>%
                    dplyr::filter(Week >= test_start_cutoff)
                if (nrow(series) == 0) {
                    return(NULL)
                }
                series
            }

            model_sources <- list(
                LC = list(
                    test = if (exists("baseline_logmxt_wide")) {
                        baseline_logmxt_wide
                    } else {
                        NULL
                    }
                ),
                MortFCNet = list(
                    test = resolve_model_wides("MortFCNet")$test
                ),
                `CNN-LSTM` = list(
                    test = resolve_model_wides("CNN-LSTM")$test
                ),
                `GNN-LSTM` = list(
                    test = resolve_model_wides("GNN-LSTM")$test
                )
            )

            for (region_name in top_region_improvement$region) {
                region_info <- top_region_improvement %>%
                    dplyr::filter(region == region_name) %>%
                    dplyr::slice(1)

                plot_parts <- list()

                obs_series <- get_mxt_series_test(observed_all, region_name)
                if (!is.null(obs_series)) {
                    plot_parts[[length(plot_parts) + 1L]] <- dplyr::mutate(
                        obs_series,
                        model = "Observed"
                    )
                }

                for (model_label in names(model_sources)) {
                    series_df <- build_region_model_series(
                        NULL,
                        model_sources[[model_label]]$test,
                        region_name,
                        model_label
                    )
                    if (!is.null(series_df)) {
                        series_df <- series_df %>%
                            dplyr::mutate(Week = as.Date(Week)) %>%
                            dplyr::filter(Week >= test_start_cutoff)
                        if (nrow(series_df) > 0) {
                            plot_parts[[length(plot_parts) + 1L]] <- series_df
                        }
                    }
                }

                region_plot_df <- dplyr::bind_rows(plot_parts) %>%
                    dplyr::mutate(
                        model = factor(
                            model,
                            levels = c(
                                "Observed",
                                "LC",
                                "MortFCNet",
                                "CNN-LSTM",
                                "GNN-LSTM"
                            )
                        )
                    )

                if (nrow(region_plot_df) == 0) {
                    next
                }

                # Plot: top-improvement region fit comparison for CNN-LSTM and GNN-LSTM.
                region_plot <- ggplot2::ggplot(
                    region_plot_df,
                    ggplot2::aes(
                        x = Week,
                        y = rate,
                        colour = model,
                        linetype = model
                    )
                ) +
                    ggplot2::annotate(
                        "rect",
                        xmin = test_start_cutoff,
                        xmax = test_end,
                        ymin = -Inf,
                        ymax = Inf,
                        alpha = 0.03,
                        fill = "steelblue"
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model != "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = model_linewidth,
                        alpha = 0.75
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model == "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = observed_linewidth
                    ) +
                    ggplot2::scale_colour_manual(
                        values = top_region_plot_colors
                    ) +
                    ggplot2::scale_linetype_manual(
                        values = c(
                            "Observed" = "solid",
                            "LC" = "dashed",
                            "MortFCNet" = "solid",
                            "CNN-LSTM" = "solid",
                            "GNN-LSTM" = "solid"
                        )
                    ) +
                    ggplot2::labs(
                        title = paste0(region_name),
                        subtitle = sprintf(
                            "Top improvement model: %s | %.2f%% vs baseline",
                            as.character(region_info$best_model[[1]]),
                            region_info$improvement_pct[[1]]
                        ),
                        x = "Week",
                        y = "Mortality rate",
                        colour = "Model",
                        linetype = "Model"
                    ) +
                    base_theme +
                    ggplot2::theme(
                        plot.subtitle = ggplot2::element_text(hjust = 0.5),
                        axis.text.x = ggplot2::element_text(
                            angle = 45,
                            hjust = 1
                        ),
                        legend.position = "top",
                        legend.direction = "horizontal"
                    )

                region_plot <- region_plot +
                    ggplot2::scale_x_date(
                        limits = as.Date(c(test_start_cutoff, test_end)),
                        date_breaks = "1 year",
                        date_labels = "%Y"
                    )

                safe_region_name <- tolower(gsub(
                    "[^A-Za-z0-9]+",
                    "_",
                    region_name
                ))
                safe_region_name <- gsub("(^_+|_+$)", "", safe_region_name)

                ggplot2::ggsave(
                    filename = file.path(
                        results_plots_dir,
                        paste0("cnn_top_", safe_region_name, ".png")
                    ),
                    plot = region_plot,
                    width = 12,
                    height = 6,
                    dpi = 300
                )
            }
        }
    } else {
        NULL
    }
}


################################################################################
# Plot: strongest GNN-LSTM region improvements against baseline
################################################################################

if (exists("per_region_metrics_table")) {
    baseline_mse <- per_region_metrics_table %>%
        dplyr::filter(model == "BASELINE") %>%
        dplyr::select(region, baseline_mse = mse_mxt_test)

    top_region_improvement_gnn <- per_region_metrics_table %>%
        dplyr::filter(model == "GNN-LSTM") %>%
        dplyr::left_join(baseline_mse, by = "region") %>%
        dplyr::mutate(
            improvement_pct = dplyr::if_else(
                is.na(baseline_mse) | baseline_mse == 0,
                NA_real_,
                (baseline_mse - mse_mxt_test) / baseline_mse * 100
            )
        ) %>%
        dplyr::filter(!is.na(improvement_pct)) %>%
        dplyr::arrange(dplyr::desc(improvement_pct)) %>%
        dplyr::slice_head(n = 3) %>%
        dplyr::mutate(best_model = as.character(model)) %>%
        dplyr::select(
            region,
            best_model,
            baseline_mse,
            model_mse = mse_mxt_test,
            improvement_pct
        )

    if (nrow(top_region_improvement_gnn) > 0) {
        if (requireNamespace("readr", quietly = TRUE)) {
            readr::write_csv(
                top_region_improvement_gnn,
                file.path(
                    results_tables_dir,
                    "top3_region_improvement_gnn_summary.csv"
                )
            )
        } else {
            write.csv(
                top_region_improvement_gnn,
                file = file.path(
                    results_tables_dir,
                    "top3_region_improvement_gnn_summary.csv"
                ),
                row.names = FALSE
            )
        }
        observed_all <- NULL
        if (exists("observed_test")) {
            observed_all <- observed_test %>%
                dplyr::mutate(Week = as.Date(Week))
        }

        if (!is.null(observed_all)) {
            test_end <- max(observed_all$Week, na.rm = TRUE)
            test_start_cutoff <- test_end - 365 * 5

            get_mxt_series_test <- function(wide_df, region) {
                series <- get_mxt_series(wide_df, region)
                if (is.null(series)) {
                    return(NULL)
                }
                series <- series %>%
                    dplyr::mutate(Week = as.Date(Week)) %>%
                    dplyr::filter(Week >= test_start_cutoff)
                if (nrow(series) == 0) {
                    return(NULL)
                }
                series
            }

            model_sources <- list(
                LC = list(
                    test = if (exists("baseline_logmxt_wide")) {
                        baseline_logmxt_wide
                    } else {
                        NULL
                    }
                ),
                MortFCNet = list(
                    test = resolve_model_wides("MortFCNet")$test
                ),
                `CNN-LSTM` = list(
                    test = resolve_model_wides("CNN-LSTM")$test
                ),
                `GNN-LSTM` = list(
                    test = resolve_model_wides("GNN-LSTM")$test
                )
            )

            for (region_name in top_region_improvement_gnn$region) {
                plot_parts <- list()

                obs_series <- get_mxt_series_test(observed_all, region_name)
                if (!is.null(obs_series)) {
                    plot_parts[[length(plot_parts) + 1L]] <- dplyr::mutate(
                        obs_series,
                        model = "Observed"
                    )
                }

                for (model_label in names(model_sources)) {
                    series_df <- build_region_model_series(
                        NULL,
                        model_sources[[model_label]]$test,
                        region_name,
                        model_label
                    )
                    if (!is.null(series_df)) {
                        series_df <- series_df %>%
                            dplyr::mutate(Week = as.Date(Week)) %>%
                            dplyr::filter(Week >= test_start_cutoff)
                        if (nrow(series_df) > 0) {
                            plot_parts[[length(plot_parts) + 1L]] <- series_df
                        }
                    }
                }

                region_plot_df <- dplyr::bind_rows(plot_parts) %>%
                    dplyr::mutate(
                        model = factor(
                            model,
                            levels = c(
                                "Observed",
                                "LC",
                                "MortFCNet",
                                "CNN-LSTM",
                                "GNN-LSTM"
                            )
                        )
                    )

                if (nrow(region_plot_df) == 0) {
                    next
                }

                # Plot: top-improvement region fit comparison for GNN-LSTM.
                region_plot <- ggplot2::ggplot(
                    region_plot_df,
                    ggplot2::aes(
                        x = Week,
                        y = rate,
                        colour = model,
                        linetype = model
                    )
                ) +
                    ggplot2::annotate(
                        "rect",
                        xmin = test_start_cutoff,
                        xmax = test_end,
                        ymin = -Inf,
                        ymax = Inf,
                        alpha = 0.03,
                        fill = "steelblue"
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model != "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = model_linewidth,
                        alpha = 0.75
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model == "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = observed_linewidth
                    ) +
                    ggplot2::scale_colour_manual(
                        values = top_region_plot_colors
                    ) +
                    ggplot2::scale_linetype_manual(
                        values = c(
                            "Observed" = "solid",
                            "LC" = "dashed",
                            "MortFCNet" = "solid",
                            "CNN-LSTM" = "solid",
                            "GNN-LSTM" = "solid"
                        )
                    ) +
                    ggplot2::labs(
                        title = paste0("GNN top region fit: ", region_name),
                        subtitle = sprintf(
                            "Top improvement model: %s | %.2f%% vs baseline",
                            "GNN-LSTM",
                            top_region_improvement_gnn %>%
                                dplyr::filter(region == region_name) %>%
                                dplyr::pull(improvement_pct)
                        ),
                        x = "Week",
                        y = "Mortality rate",
                        colour = "Model",
                        linetype = "Model"
                    ) +
                    ggplot2::theme_minimal(base_size = 12) +
                    ggplot2::theme(
                        plot.title = ggplot2::element_text(hjust = 0.5),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5),
                        axis.text.x = ggplot2::element_text(
                            angle = 45,
                            hjust = 1
                        ),
                        legend.position = "bottom"
                    ) +
                    ggplot2::scale_x_date(
                        limits = as.Date(c(test_start_cutoff, test_end)),
                        date_breaks = "1 year",
                        date_labels = "%Y"
                    )

                safe_region_name <- tolower(gsub(
                    "[^A-Za-z0-9]+",
                    "_",
                    region_name
                ))
                safe_region_name <- gsub("(^_+|_+$)", "", safe_region_name)

                ggplot2::ggsave(
                    filename = file.path(
                        results_plots_dir,
                        paste0("gnn_top_", safe_region_name, ".png")
                    ),
                    plot = region_plot,
                    width = 12,
                    height = 6,
                    dpi = 300
                )
            }
        }
    } else {
        NULL
    }
}

################################################################################
# Plot: age-specific best-region comparisons for ages 80 and 90+
################################################################################
if (
    exists("per_age_mse_summary") &&
        exists("observed_test") &&
        exists("baseline_logmxt_wide") &&
        exists("results_plots_dir")
) {
    focus_age_groups <- c(80L, 90L)
    focus_model_labels <- c("BASELINE", "MortFCNet", "CNN-LSTM", "GNN-LSTM")
    focus_model_colors <- c(
        "BASELINE" = lee_carter_blue,
        "MortFCNet" = mortfcnet_green,
        "CNN-LSTM" = cnn_lstm_orange,
        "GNN-LSTM" = gnn_lstm_violet
    )

    age_focus_label <- function(age_value) {
        age_value <- as.integer(age_value)
        ifelse(age_value == 90L, "90+", as.character(age_value))
    }

    extract_age_region_series <- function(wide_df, age_value, region_value) {
        if (
            is.null(wide_df) ||
                !is.data.frame(wide_df) ||
                !("Week" %in% names(wide_df))
        ) {
            return(NULL)
        }

        age_cols <- get_age_metric_cols(wide_df, age_value)
        region_cols <- get_region_metric_cols(wide_df, region_value)
        series_cols <- intersect(age_cols, region_cols)
        if (length(series_cols) == 0) {
            return(NULL)
        }

        series_df <- wide_df[, c("Week", series_cols[1]), drop = FALSE]
        names(series_df)[2] <- "logmxt"
        series_df$Week <- as.Date(series_df$Week)
        series_df$logmxt <- as.numeric(series_df$logmxt)
        series_df <- series_df[is.finite(series_df$logmxt), , drop = FALSE]
        series_df[order(series_df$Week), , drop = FALSE]
    }

    build_age_region_plot_data <- function(age_value, region_values) {
        model_wides <- list(
            "Observed" = observed_test,
            "LC" = baseline_logmxt_wide,
            "MortFCNet" = resolve_model_wides("MortFCNet")$test,
            "CNN-LSTM" = resolve_model_wides("CNN-LSTM")$test,
            "GNN-LSTM" = resolve_model_wides("GNN-LSTM")$test
        )

        plot_parts <- list()
        for (region_value in region_values) {
            for (model_label in names(model_wides)) {
                series_df <- extract_age_region_series(
                    model_wides[[model_label]],
                    age_value,
                    region_value
                )
                if (is.null(series_df) || nrow(series_df) == 0) {
                    next
                }

                plot_parts[[length(plot_parts) + 1L]] <- data.frame(
                    Week = series_df$Week,
                    rate = exp(series_df$logmxt),
                    model = model_label,
                    region = region_value,
                    stringsAsFactors = FALSE
                )
            }
        }

        dplyr::bind_rows(plot_parts)
    }

    age_focus_summary <- per_age_mse_summary %>%
        dplyr::filter(
            age %in% focus_age_groups,
            model %in% focus_model_labels
        ) %>%
        dplyr::mutate(
            age_label = factor(
                age_focus_label(age),
                levels = c("80", "90+")
            ),
            model = factor(model, levels = focus_model_labels)
        ) %>%
        dplyr::arrange(age, model)

    if (nrow(age_focus_summary) > 0) {
        age_focus_plot_dir <- file.path(results_plots_dir, "age_focus")
        dir.create(age_focus_plot_dir, recursive = TRUE, showWarnings = FALSE)

        age_focus_tables_dir <- results_tables_dir

        if (requireNamespace("readr", quietly = TRUE)) {
            readr::write_csv(
                age_focus_summary,
                file.path(
                    age_focus_tables_dir,
                    "age_80_90plus_model_summary.csv"
                )
            )
        } else {
            write.csv(
                age_focus_summary,
                file = file.path(
                    age_focus_tables_dir,
                    "age_80_90plus_model_summary.csv"
                ),
                row.names = FALSE
            )
        }

        region_candidates <- if (exists("regions")) {
            regions
        } else {
            sort(unique(per_region_metrics_table$region))
        }

        model_test_wides <- list(
            "MortFCNet" = resolve_model_wides("MortFCNet")$test,
            "CNN-LSTM" = resolve_model_wides("CNN-LSTM")$test,
            "GNN-LSTM" = resolve_model_wides("GNN-LSTM")$test
        )

        selected_region_values <- if (exists("selected_regions")) {
            as.character(selected_regions)
        } else {
            c("FRK2", "FRJ2")
        }
        selected_age_region_pairs <- data.frame(
            age = rep(90L, length(selected_region_values)),
            region = selected_region_values,
            stringsAsFactors = FALSE
        )
        selected_age_region_pairs$age_label <- age_focus_label(
            selected_age_region_pairs$age
        )

        age_region_perf_rows <- list()
        for (age_value in focus_age_groups) {
            age_label <- age_focus_label(age_value)
            for (region_value in region_candidates) {
                obs_series <- extract_age_region_series(
                    observed_test,
                    age_value,
                    region_value
                )
                baseline_series <- extract_age_region_series(
                    baseline_logmxt_wide,
                    age_value,
                    region_value
                )

                if (
                    is.null(obs_series) ||
                        is.null(baseline_series) ||
                        nrow(obs_series) == 0 ||
                        nrow(baseline_series) == 0
                ) {
                    next
                }

                baseline_metrics <- compute_metrics(
                    baseline_series,
                    obs_series
                )

                model_mse_values <- vapply(
                    focus_model_labels[-1L],
                    function(model_label) {
                        model_series <- extract_age_region_series(
                            model_test_wides[[model_label]],
                            age_value,
                            region_value
                        )
                        if (is.null(model_series) || nrow(model_series) == 0) {
                            return(NA_real_)
                        }
                        compute_metrics(model_series, obs_series)$mse_mxt
                    },
                    numeric(1)
                )
                names(model_mse_values) <- focus_model_labels[-1L]

                if (all(!is.finite(model_mse_values))) {
                    next
                }

                model_mse_for_best <- model_mse_values
                model_mse_for_best[!is.finite(model_mse_for_best)] <- Inf
                best_idx <- which.min(model_mse_for_best)
                best_model <- names(model_mse_values)[best_idx]
                best_model_mse <- model_mse_values[best_idx]

                age_region_perf_rows[[
                    length(age_region_perf_rows) + 1L
                ]] <- data.frame(
                    age = age_value,
                    age_label = age_label,
                    region = region_value,
                    baseline_mse = baseline_metrics$mse_mxt,
                    MortFCNet_mse = model_mse_values["MortFCNet"],
                    CNN_LSTM_mse = model_mse_values["CNN-LSTM"],
                    GNN_LSTM_mse = model_mse_values["GNN-LSTM"],
                    best_model = best_model,
                    best_model_mse = best_model_mse,
                    improvement_pct = ifelse(
                        is.finite(baseline_metrics$mse_mxt) &&
                            baseline_metrics$mse_mxt > 0 &&
                            is.finite(best_model_mse),
                        (baseline_metrics$mse_mxt - best_model_mse) /
                            baseline_metrics$mse_mxt *
                            100,
                        NA_real_
                    ),
                    stringsAsFactors = FALSE
                )
            }
        }

        age_region_perf_table <- dplyr::bind_rows(age_region_perf_rows) %>%
            dplyr::arrange(age, dplyr::desc(improvement_pct), region)

        if (nrow(age_region_perf_table) > 0) {
            if (requireNamespace("readr", quietly = TRUE)) {
                readr::write_csv(
                    age_region_perf_table,
                    file.path(
                        age_focus_tables_dir,
                        "age_80_90plus_region_model_summary.csv"
                    )
                )
            } else {
                write.csv(
                    age_region_perf_table,
                    file = file.path(
                        age_focus_tables_dir,
                        "age_80_90plus_region_model_summary.csv"
                    ),
                    row.names = FALSE
                )
            }

            age_region_top3 <- age_region_perf_table %>%
                dplyr::filter(is.finite(improvement_pct)) %>%
                dplyr::group_by(age, age_label) %>%
                dplyr::slice_max(improvement_pct, n = 3, with_ties = FALSE) %>%
                dplyr::ungroup() %>%
                dplyr::arrange(age, dplyr::desc(improvement_pct))

            if (nrow(age_region_top3) > 0) {
                if (requireNamespace("readr", quietly = TRUE)) {
                    readr::write_csv(
                        age_region_top3,
                        file.path(
                            age_focus_tables_dir,
                            "age_80_90plus_top3_regions.csv"
                        )
                    )
                } else {
                    write.csv(
                        age_region_top3,
                        file = file.path(
                            age_focus_tables_dir,
                            "age_80_90plus_top3_regions.csv"
                        ),
                        row.names = FALSE
                    )
                }

                top3_plot_dir <- file.path(age_focus_plot_dir, "top3_regions")
                dir.create(
                    top3_plot_dir,
                    recursive = TRUE,
                    showWarnings = FALSE
                )

                for (age_value in focus_age_groups) {
                    age_label <- age_focus_label(age_value)
                    selected_regions <- age_region_top3 %>%
                        dplyr::filter(age == age_value) %>%
                        dplyr::arrange(dplyr::desc(improvement_pct)) %>%
                        dplyr::pull(region)

                    if (length(selected_regions) == 0) {
                        next
                    }

                    top3_plot_df <- build_age_region_plot_data(
                        age_value,
                        selected_regions
                    )

                    if (nrow(top3_plot_df) == 0) {
                        next
                    }

                    top3_plot_df$region <- factor(
                        top3_plot_df$region,
                        levels = selected_regions
                    )
                    top3_plot_df$model <- factor(
                        top3_plot_df$model,
                        levels = c(
                            "Observed",
                            "LC",
                            "MortFCNet",
                            "CNN-LSTM",
                            "GNN-LSTM"
                        )
                    )

                    age_top3_plot <- ggplot2::ggplot(
                        top3_plot_df,
                        ggplot2::aes(
                            x = Week,
                            y = rate,
                            colour = model,
                            linetype = model
                        )
                    ) +
                        ggplot2::geom_line(
                            data = dplyr::filter(
                                top3_plot_df,
                                model != "Observed"
                            ),
                            linewidth = 0.7,
                            alpha = 1
                        ) +
                        ggplot2::geom_line(
                            data = dplyr::filter(
                                top3_plot_df,
                                model == "Observed"
                            ),
                            linewidth = 0.5,
                            alpha = 0.5
                        ) +
                        ggplot2::facet_wrap(
                            ~region,
                            ncol = 1,
                            scales = "free_y"
                        ) +
                        ggplot2::scale_colour_manual(
                            values = top_region_plot_colors
                        ) +
                        ggplot2::scale_linetype_manual(
                            values = c(
                                "Observed" = "solid",
                                "LC" = "dashed",
                                "MortFCNet" = "dashed",
                                "CNN-LSTM" = "solid",
                                "GNN-LSTM" = "solid"
                            )
                        ) +
                        ggplot2::scale_x_date(
                            date_breaks = "1 year",
                            date_labels = "%Y"
                        ) +
                        ggplot2::labs(
                            title = paste0(
                                "Age ",
                                age_label,
                                ": top 3 regions"
                            ),
                            subtitle = paste0(
                                "Ranked by improvement vs LC baseline using ",
                                "MortFCNet, CNN-LSTM, and GNN-LSTM"
                            ),
                            x = "Week",
                            y = "Mortality rate",
                            colour = "Model",
                            linetype = "Model"
                        ) +
                        base_theme +
                        ggplot2::theme(
                            legend.position = "bottom",
                            plot.subtitle = ggplot2::element_text(hjust = 0.5),
                            strip.text = ggplot2::element_text(face = "bold")
                        )

                    ggplot2::ggsave(
                        filename = file.path(
                            top3_plot_dir,
                            paste0(
                                "age_",
                                ifelse(
                                    age_value == 90L,
                                    "90plus",
                                    as.character(age_value)
                                ),
                                "_top3_regions.png"
                            )
                        ),
                        plot = age_top3_plot,
                        width = 12,
                        height = 10,
                        dpi = 300
                    )
                }
            }
        }

        build_age_region_metrics_table <- function(age_values, region_values) {
            rows <- list()

            for (age_value in age_values) {
                age_label <- age_focus_label(age_value)

                for (region_value in region_values) {
                    obs_series <- extract_age_region_series(
                        observed_test,
                        age_value,
                        region_value
                    )
                    baseline_series <- extract_age_region_series(
                        baseline_logmxt_wide,
                        age_value,
                        region_value
                    )

                    if (
                        is.null(obs_series) ||
                            is.null(baseline_series) ||
                            nrow(obs_series) == 0 ||
                            nrow(baseline_series) == 0
                    ) {
                        next
                    }

                    baseline_metrics <- compute_metrics(
                        baseline_series,
                        obs_series
                    )

                    model_mse_values <- vapply(
                        focus_model_labels[-1L],
                        function(model_label) {
                            model_series <- extract_age_region_series(
                                model_test_wides[[model_label]],
                                age_value,
                                region_value
                            )
                            if (
                                is.null(model_series) ||
                                    nrow(model_series) == 0
                            ) {
                                return(NA_real_)
                            }
                            compute_metrics(model_series, obs_series)$mse_mxt
                        },
                        numeric(1)
                    )
                    names(model_mse_values) <- focus_model_labels[-1L]

                    if (all(!is.finite(model_mse_values))) {
                        next
                    }

                    model_mse_for_best <- model_mse_values
                    model_mse_for_best[!is.finite(model_mse_for_best)] <- Inf
                    best_idx <- which.min(model_mse_for_best)
                    best_model <- names(model_mse_values)[best_idx]
                    best_model_mse <- model_mse_values[best_idx]

                    rows[[length(rows) + 1L]] <- data.frame(
                        age = age_value,
                        age_label = age_label,
                        region = region_value,
                        baseline_mse = baseline_metrics$mse_mxt,
                        MortFCNet_mse = model_mse_values["MortFCNet"],
                        CNN_LSTM_mse = model_mse_values["CNN-LSTM"],
                        GNN_LSTM_mse = model_mse_values["GNN-LSTM"],
                        best_model = best_model,
                        best_model_mse = best_model_mse,
                        improvement_pct = ifelse(
                            is.finite(baseline_metrics$mse_mxt) &&
                                baseline_metrics$mse_mxt > 0 &&
                                is.finite(best_model_mse),
                            (baseline_metrics$mse_mxt - best_model_mse) /
                                baseline_metrics$mse_mxt *
                                100,
                            NA_real_
                        ),
                        stringsAsFactors = FALSE
                    )
                }
            }

            dplyr::bind_rows(rows)
        }

        all_age_values <- sort(unique(stats::na.omit(
            as.integer(as.character(per_age_mse_summary$age))
        )))

        if (length(all_age_values) > 0) {
            all_age_region_perf_table <- build_age_region_metrics_table(
                all_age_values,
                region_candidates
            )

            if (nrow(all_age_region_perf_table) > 0) {
                best_region_age_rows <- all_age_region_perf_table %>%
                    dplyr::filter(is.finite(improvement_pct)) %>%
                    dplyr::group_by(region) %>%
                    dplyr::slice_max(
                        improvement_pct,
                        n = 1,
                        with_ties = FALSE
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::arrange(dplyr::desc(improvement_pct), region) %>%
                    dplyr::mutate(
                        region_age_label = paste0(
                            region,
                            " | Age ",
                            age_label
                        )
                    )

                if (nrow(best_region_age_rows) > 0) {
                    best_age_region_plot_dir <- file.path(
                        results_plots_dir,
                        "age_focus",
                        "best_age_by_region"
                    )
                    dir.create(
                        best_age_region_plot_dir,
                        recursive = TRUE,
                        showWarnings = FALSE
                    )

                    best_region_age_summary <- best_region_age_rows %>%
                        dplyr::select(
                            region,
                            age,
                            age_label,
                            best_model,
                            baseline_mse,
                            best_model_mse,
                            improvement_pct
                        )

                    if (requireNamespace("readr", quietly = TRUE)) {
                        readr::write_csv(
                            best_region_age_summary,
                            file.path(
                                results_tables_dir,
                                "best_age_by_region_summary.csv"
                            )
                        )
                    } else {
                        write.csv(
                            best_region_age_summary,
                            file = file.path(
                                results_tables_dir,
                                "best_age_by_region_summary.csv"
                            ),
                            row.names = FALSE
                        )
                    }

                    for (i in seq_len(nrow(best_region_age_rows))) {
                        region_value <- as.character(best_region_age_rows$region[
                            i
                        ])
                        age_value <- as.integer(best_region_age_rows$age[i])
                        age_label <- as.character(best_region_age_rows$age_label[
                            i
                        ])

                        region_plot_df <- build_age_region_plot_data(
                            age_value,
                            region_value
                        )
                        if (
                            is.null(region_plot_df) || nrow(region_plot_df) == 0
                        ) {
                            next
                        }

                        region_plot_df$model <- factor(
                            region_plot_df$model,
                            levels = c(
                                "Observed",
                                "LC",
                                "MortFCNet",
                                "CNN-LSTM",
                                "GNN-LSTM"
                            )
                        )

                        region_plot <- ggplot2::ggplot(
                            region_plot_df,
                            ggplot2::aes(
                                x = Week,
                                y = rate,
                                colour = model,
                                linetype = model
                            )
                        ) +
                            ggplot2::annotate(
                                "rect",
                                xmin = test_start_cutoff,
                                xmax = test_end,
                                ymin = -Inf,
                                ymax = Inf,
                                alpha = 0.03,
                                fill = "steelblue"
                            ) +
                            ggplot2::geom_line(
                                data = dplyr::filter(
                                    region_plot_df,
                                    model != "Observed"
                                ),
                                na.rm = TRUE,
                                linewidth = model_linewidth,
                                alpha = 1
                            ) +
                            ggplot2::geom_line(
                                data = dplyr::filter(
                                    region_plot_df,
                                    model == "Observed"
                                ),
                                na.rm = TRUE,
                                linewidth = 0.5,
                                alpha = 0.5
                            ) +
                            ggplot2::scale_colour_manual(
                                values = top_region_plot_colors
                            ) +
                            ggplot2::scale_linetype_manual(
                                values = c(
                                    "Observed" = "solid",
                                    "LC" = "dashed",
                                    "MortFCNet" = "dashed",
                                    "CNN-LSTM" = "solid",
                                    "GNN-LSTM" = "solid"
                                )
                            ) +
                            ggplot2::scale_x_date(
                                limits = as.Date(c(
                                    test_start_cutoff,
                                    test_end
                                )),
                                date_breaks = "1 year",
                                date_labels = "%Y"
                            ) +
                            ggplot2::labs(
                                title = paste0("Region ", region_value),
                                subtitle = sprintf(
                                    "Best age group: %s | Top improvement model: %s | %.2f%% vs baseline",
                                    age_label,
                                    as.character(best_region_age_rows$best_model[
                                        i
                                    ]),
                                    best_region_age_rows$improvement_pct[i]
                                ),
                                x = "Week",
                                y = "Mortality rate",
                                colour = "Model",
                                linetype = "Model"
                            ) +
                            base_theme +
                            ggplot2::theme(
                                plot.subtitle = ggplot2::element_text(
                                    hjust = 0.5
                                ),
                                axis.text.x = ggplot2::element_text(
                                    angle = 45,
                                    hjust = 1
                                ),
                                legend.position = "top",
                                legend.direction = "horizontal"
                            )

                        safe_region_name <- tolower(gsub(
                            "[^A-Za-z0-9]+",
                            "_",
                            region_value
                        ))
                        safe_region_name <- gsub(
                            "(^_+|_+$)",
                            "",
                            safe_region_name
                        )

                        safe_age_label <- ifelse(
                            age_label == "90+",
                            "90plus",
                            age_label
                        )
                        safe_age_label <- gsub(
                            "[^A-Za-z0-9]+",
                            "_",
                            safe_age_label
                        )
                        safe_age_label <- gsub(
                            "(^_+|_+$)",
                            "",
                            safe_age_label
                        )

                        plot_file <- file.path(
                            best_age_region_plot_dir,
                            paste0(
                                "region_",
                                safe_region_name,
                                "_age",
                                safe_age_label,
                                "_best_age_group.png"
                            )
                        )

                        ggplot2::ggsave(
                            filename = plot_file,
                            plot = region_plot,
                            width = 12,
                            height = 6,
                            dpi = 300
                        )
                    }
                }
            }
        }
    }
    if (
        exists("all_age_region_perf_table") &&
            nrow(all_age_region_perf_table) > 0
    ) {
        save_best_age_region_model_outputs <- function(
            model_label,
            model_mse_col,
            plot_subdir
        ) {
            model_best_rows <- all_age_region_perf_table %>%
                dplyr::mutate(
                    model_mse = .data[[model_mse_col]],
                    improvement_pct_model = dplyr::if_else(
                        is.finite(baseline_mse) &
                            baseline_mse > 0 &
                            is.finite(model_mse),
                        (baseline_mse - model_mse) / baseline_mse * 100,
                        NA_real_
                    )
                ) %>%
                dplyr::filter(is.finite(improvement_pct_model)) %>%
                dplyr::group_by(region) %>%
                dplyr::slice_max(
                    improvement_pct_model,
                    n = 1,
                    with_ties = FALSE
                ) %>%
                dplyr::ungroup() %>%
                dplyr::arrange(dplyr::desc(improvement_pct_model), region)

            if (nrow(model_best_rows) == 0) {
                return(invisible(NULL))
            }

            model_plot_dir <- file.path(
                results_plots_dir,
                "age_focus",
                "best_age_by_region",
                plot_subdir
            )
            dir.create(model_plot_dir, recursive = TRUE, showWarnings = FALSE)

            safe_model_name <- tolower(gsub("[^A-Za-z0-9]+", "_", model_label))
            safe_model_name <- gsub("(^_+|_+$)", "", safe_model_name)

            model_summary <- model_best_rows %>%
                dplyr::transmute(
                    model = model_label,
                    region,
                    age,
                    age_label,
                    baseline_mse,
                    best_model_mse = model_mse,
                    improvement_pct = improvement_pct_model
                )

            if (requireNamespace("readr", quietly = TRUE)) {
                readr::write_csv(
                    model_summary,
                    file.path(
                        results_tables_dir,
                        paste0(
                            "best_age_by_region_",
                            safe_model_name,
                            "_summary.csv"
                        )
                    )
                )
            } else {
                write.csv(
                    model_summary,
                    file = file.path(
                        results_tables_dir,
                        paste0(
                            "best_age_by_region_",
                            safe_model_name,
                            "_summary.csv"
                        )
                    ),
                    row.names = FALSE
                )
            }

            for (i in seq_len(nrow(model_best_rows))) {
                region_value <- as.character(model_best_rows$region[i])
                age_value <- as.integer(model_best_rows$age[i])
                age_label <- as.character(model_best_rows$age_label[i])

                region_plot_df <- build_age_region_plot_data(
                    age_value,
                    region_value
                )
                if (is.null(region_plot_df) || nrow(region_plot_df) == 0) {
                    next
                }

                region_plot_df$model <- factor(
                    region_plot_df$model,
                    levels = c(
                        "Observed",
                        "LC",
                        "MortFCNet",
                        "CNN-LSTM",
                        "GNN-LSTM"
                    )
                )

                safe_region_name <- tolower(gsub(
                    "[^A-Za-z0-9]+",
                    "_",
                    region_value
                ))
                safe_region_name <- gsub("(^_+|_+$)", "", safe_region_name)

                safe_age_label <- ifelse(
                    age_label == "90+",
                    "90plus",
                    age_label
                )
                safe_age_label <- gsub("[^A-Za-z0-9]+", "_", safe_age_label)
                safe_age_label <- gsub("(^_+|_+$)", "", safe_age_label)

                region_plot <- ggplot2::ggplot(
                    region_plot_df,
                    ggplot2::aes(
                        x = Week,
                        y = rate,
                        colour = model,
                        linetype = model
                    )
                ) +
                    ggplot2::annotate(
                        "rect",
                        xmin = test_start_cutoff,
                        xmax = test_end,
                        ymin = -Inf,
                        ymax = Inf,
                        alpha = 0.03,
                        fill = "steelblue"
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model != "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = model_linewidth,
                        alpha = 1
                    ) +
                    ggplot2::geom_line(
                        data = dplyr::filter(
                            region_plot_df,
                            model == "Observed"
                        ),
                        na.rm = TRUE,
                        linewidth = 0.5,
                        alpha = 0.5
                    ) +
                    ggplot2::scale_colour_manual(
                        values = top_region_plot_colors
                    ) +
                    ggplot2::scale_linetype_manual(
                        values = c(
                            "Observed" = "solid",
                            "LC" = "dashed",
                            "MortFCNet" = "dashed",
                            "CNN-LSTM" = "solid",
                            "GNN-LSTM" = "solid"
                        )
                    ) +
                    ggplot2::scale_x_date(
                        limits = as.Date(c(test_start_cutoff, test_end)),
                        date_breaks = "1 year",
                        date_labels = "%Y"
                    ) +
                    ggplot2::labs(
                        title = paste0(
                            "Region ",
                            region_value,
                            " - ",
                            model_label
                        ),
                        subtitle = sprintf(
                            "Best age group: %s | Improvement vs baseline: %.2f%%",
                            age_label,
                            model_best_rows$improvement_pct_model[i]
                        ),
                        x = "Week",
                        y = "Mortality rate",
                        colour = "Model",
                        linetype = "Model"
                    ) +
                    base_theme +
                    ggplot2::theme(
                        plot.subtitle = ggplot2::element_text(hjust = 0.5),
                        axis.text.x = ggplot2::element_text(
                            angle = 45,
                            hjust = 1
                        ),
                        legend.position = "top",
                        legend.direction = "horizontal"
                    )

                plot_file <- file.path(
                    model_plot_dir,
                    paste0(
                        "region_",
                        safe_region_name,
                        "_age",
                        safe_age_label,
                        safe_model_name,
                        "_best_age_group.png"
                    )
                )

                ggplot2::ggsave(
                    filename = plot_file,
                    plot = region_plot,
                    width = 12,
                    height = 6,
                    dpi = 300
                )
            }
        }

        save_best_age_region_model_outputs(
            "CNN-LSTM",
            "CNN_LSTM_mse",
            "cnn_lstm"
        )
        save_best_age_region_model_outputs(
            "GNN-LSTM",
            "GNN_LSTM_mse",
            "gnn_lstm"
        )
    }
}


# ============================================================================
# FIGURE 8: CLIMATE FEATURE IMPORTANCE
# Publication-quality visualization
# ============================================================================

shap_plot_dir <- file.path("Results", "plots", "shap")
dir.create(shap_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ─── Feature metadata ────────────────────────────────────────────────────────
feature_info <- expand.grid(
    lag = seq_len(seq_length_cnn_gnn) - 1L,
    var = measurement_cols,
    stringsAsFactors = FALSE
) %>%
    dplyr::mutate(
        feat_name = paste0("lag", lag, "_", var),
        t_step = seq_length_cnn_gnn - lag,
        f_idx = match(var, measurement_cols),
        category = dplyr::case_when(
            var %in% c("tn", "tg", "tx") ~ "Temperature",
            var %in% c("rr", "hu") ~ "Precipitation & Humidity",
            var == "fg" ~ "Wind",
            TRUE ~ "Other"
        ),
        var_label = dplyr::case_when(
            var == "tn" ~ "T_min",
            var == "tg" ~ "T_mean",
            var == "tx" ~ "T_max",
            var == "rr" ~ "Precipitation",
            var == "hu" ~ "Humidity",
            var == "fg" ~ "Wind speed",
            TRUE ~ var
        ),
        feat_label = paste0(var_label, "  (lag ", lag, ")")
    )

n_feats <- nrow(feature_info)

# ─── Colorblind-friendly palette ────────────────────────────────────────────
cat_colors <- c(
    "Temperature" = "#D55E00",
    "Precipitation & Humidity" = "#4D4D4D",
    "Wind" = "#009E73"
)

B_perm <- 30L

compute_feature_importance <- function(
    model,
    model_label,
    x_fit,
    region_fit,
    age_fit,
    age_mid_fit,
    y_scale_vec,
    y_center_vec,
    intercept_fix_vec,
    year_numeric_fit = NULL
) {
    set.seed(global_seed + if (identical(model_label, "CNN-LSTM")) 41L else 42L)

    n_eval_local <- min(60L, dim(x_fit)[1])
    n_bg_local <- min(120L, dim(x_fit)[1])

    all_idx <- seq_len(dim(x_fit)[1])
    eval_idx <- sample(all_idx, n_eval_local)
    bg_pool <- setdiff(all_idx, eval_idx)
    if (length(bg_pool) < n_bg_local) {
        bg_pool <- all_idx
    }
    bg_idx <- sample(
        bg_pool,
        n_bg_local,
        replace = length(bg_pool) < n_bg_local
    )

    eval_region <- region_fit[eval_idx, 1]
    eval_age <- age_fit[eval_idx, 1]
    eval_age_mid_z <- matrix(age_mid_fit[eval_idx, 1], ncol = 1)
    eval_col_names <- target_cols[as.integer(eval_age)]

    candidate_inputs <- list(
        climate_input = x_fit[eval_idx, , , drop = FALSE],
        region_input = matrix(as.integer(eval_region), ncol = 1),
        age_input = matrix(as.integer(eval_age), ncol = 1),
        age_numeric_input = eval_age_mid_z,
        year_numeric_input = if (is.null(year_numeric_fit)) {
            NULL
        } else {
            matrix(year_numeric_fit[eval_idx, 1], ncol = 1)
        }
    )

    baseline_inputs <- select_model_inputs(model, candidate_inputs)
    baseline_pred_z <- as.numeric(predict(model, baseline_inputs))
    baseline_pred <- baseline_pred_z *
        y_scale_vec[eval_col_names] +
        y_center_vec[eval_col_names] +
        intercept_fix_vec[eval_col_names]

    importance_mat <- matrix(NA_real_, nrow = B_perm, ncol = n_feats)
    colnames(importance_mat) <- feature_info$feat_name
    for (b in seq_len(B_perm)) {
        for (fi in seq_len(n_feats)) {
            ts <- feature_info$t_step[fi]
            fs <- feature_info$f_idx[fi]

            X_perturbed <- x_fit[eval_idx, , , drop = FALSE]
            X_perturbed[, ts, fs] <- sample(
                x_fit[bg_idx, ts, fs],
                size = n_eval_local,
                replace = TRUE
            )

            pert_candidate_inputs <- candidate_inputs
            pert_candidate_inputs$climate_input <- X_perturbed

            pert_inputs <- select_model_inputs(model, pert_candidate_inputs)
            pz <- as.numeric(predict(model, pert_inputs))
            pp <- pz *
                y_scale_vec[eval_col_names] +
                y_center_vec[eval_col_names] +
                intercept_fix_vec[eval_col_names]

            importance_mat[b, fi] <- mean(abs(pp - baseline_pred), na.rm = TRUE)
        }
    }

    imp_mean <- colMeans(importance_mat, na.rm = TRUE)
    imp_lo <- apply(importance_mat, 2, quantile, 0.025, na.rm = TRUE)
    imp_hi <- apply(importance_mat, 2, quantile, 0.975, na.rm = TRUE)

    list(
        imp_mean = imp_mean,
        imp_lo = imp_lo,
        imp_hi = imp_hi
    )
}

plot_feature_importance <- function(
    model_label,
    imp_mean,
    imp_lo,
    imp_hi,
    plot_file,
    y_max = NULL
) {
    total <- sum(imp_mean)
    pct <- function(x) x / total * 100

    plot_df <- feature_info %>%
        dplyr::mutate(
            importance = pct(imp_mean[feat_name]),
            ci_lo = pct(imp_lo[feat_name]),
            ci_hi = pct(imp_hi[feat_name])
        ) %>%
        dplyr::arrange(dplyr::desc(importance)) %>%
        dplyr::slice_head(n = 15) %>%
        dplyr::mutate(
            feat_label = factor(feat_label, levels = rev(feat_label)),
            category = factor(
                category,
                levels = c(
                    "Temperature",
                    "Precipitation & Humidity",
                    "Wind"
                )
            )
        )

    fig_plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = feat_label,
            y = importance,
            fill = category
        )
    ) +
        ggplot2::geom_col(
            width = 0.68,
            color = "white",
            linewidth = 0.25
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = ci_lo,
                ymax = ci_hi
            ),
            width = 0.18,
            linewidth = 0.55,
            color = "black"
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
            values = cat_colors,
            name = "Climate category"
        ) +
        ggplot2::scale_y_continuous(
            labels = function(x) sprintf("%.0f%%", x),
            expand = ggplot2::expansion(mult = c(0, 0.05)),
            limits = if (!is.null(y_max)) c(0, y_max) else NULL
        ) +
        ggplot2::labs(
            title = paste0("Climate Feature Importance (", model_label, ")"),
            subtitle = paste0(
                "Top 15 climate predictors ranked by permutation importance"
            ),
            x = NULL,
            y = "Relative importance (%)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = 15,
                face = "bold",
                margin = ggplot2::margin(b = 4)
            ),
            plot.subtitle = ggplot2::element_text(
                size = 11,
                color = "grey30",
                margin = ggplot2::margin(b = 12)
            ),
            axis.title.x = ggplot2::element_text(
                size = 11,
                margin = ggplot2::margin(t = 8)
            ),
            axis.text.y = ggplot2::element_text(
                size = 10,
                color = "black"
            ),
            axis.text.x = ggplot2::element_text(
                size = 9,
                color = "black"
            ),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(
                color = "grey88",
                linewidth = 0.35
            ),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.title = ggplot2::element_text(
                size = 10,
                face = "bold"
            ),
            legend.text = ggplot2::element_text(
                size = 9
            ),
            plot.caption = ggplot2::element_text(
                size = 8,
                color = "grey40",
                hjust = 0
            ),
            plot.margin = ggplot2::margin(
                t = 12,
                r = 18,
                b = 10,
                l = 12
            )
        )

    ggplot2::ggsave(
        filename = plot_file,
        plot = fig_plot,
        width = 9.5,
        height = 7,
        dpi = 600,
        bg = "white"
    )
    invisible(fig_plot)
}

fig8_plot_files <- list(
    "CNN-LSTM" = file.path(
        shap_plot_dir,
        "figure_8_cnn_lstm_climate_feature_importance.png"
    ),
    "GNN-LSTM" = file.path(
        shap_plot_dir,
        "figure_8_gnn_lstm_climate_feature_importance.png"
    )
)

# Compute importance for both models
imp_cnn <- compute_feature_importance(
    model = mdl_cnn,
    model_label = "CNN-LSTM",
    x_fit = X_fit_aa,
    region_fit = region_fit_aa,
    age_fit = age_fit_aa,
    age_mid_fit = age_mid_fit_aa,
    y_scale_vec = y_scale,
    y_center_vec = y_center,
    intercept_fix_vec = intercept_fix_per_age
)

imp_gnn <- compute_feature_importance(
    model = mdl_gnn,
    model_label = "GNN-LSTM",
    x_fit = X_fit_gnn,
    region_fit = region_idx_fit_gnn,
    age_fit = age_idx_fit_gnn,
    age_mid_fit = age_mid_fit_gnn,
    y_scale_vec = y_scale,
    y_center_vec = y_center,
    intercept_fix_vec = intercept_fix_gnn,
    year_numeric_fit = year_numeric_fit_gnn
)

# Compute common maximum y-value across both models
total_cnn <- sum(imp_cnn$imp_mean)
total_gnn <- sum(imp_gnn$imp_mean)
pct_cnn <- function(x) x / total_cnn * 100
pct_gnn <- function(x) x / total_gnn * 100

max_y_cnn <- max(pct_cnn(imp_cnn$imp_hi), na.rm = TRUE)
max_y_gnn <- max(pct_gnn(imp_gnn$imp_hi), na.rm = TRUE)
common_max_y <- max(max_y_cnn, max_y_gnn) * 1.1

# Plot both models with common y-axis
fig8_plots <- list(
    "CNN-LSTM" = plot_feature_importance(
        model_label = "CNN-LSTM",
        imp_mean = imp_cnn$imp_mean,
        imp_lo = imp_cnn$imp_lo,
        imp_hi = imp_cnn$imp_hi,
        plot_file = fig8_plot_files[["CNN-LSTM"]],
        y_max = common_max_y
    ),
    "GNN-LSTM" = plot_feature_importance(
        model_label = "GNN-LSTM",
        imp_mean = imp_gnn$imp_mean,
        imp_lo = imp_gnn$imp_lo,
        imp_hi = imp_gnn$imp_hi,
        plot_file = fig8_plot_files[["GNN-LSTM"]],
        y_max = common_max_y
    )
)


# ============================================================================
# PHASE 10: RELATIVE RISK (FIGURE 10)
# ============================================================================
phase10_plot_dir <- file.path("Results", "plots", "relative_risk")
dir.create(phase10_plot_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(global_seed + 52L)

B_rr <- as.integer(Sys.getenv("RR_BOOT_B", "1000"))

normalize_key <- function(x) {
    x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
    gsub("[^a-z0-9]+", "", x)
}

region_lookup <- shapef %>%
    dplyr::transmute(
        region = as.character(NUTS_ID),
        region_name = as.character(NUTS_NAME),
        region_key = normalize_key(region_name),
        region_label = iconv(region_name, to = "ASCII//TRANSLIT")
    ) %>%
    dplyr::distinct(region, .keep_all = TRUE)

target_region_names <- c(
    "Lorraine",
    "Ile-de-France",
    "Auvergne"
)

target_region_keys <- normalize_key(target_region_names)
target_region_ids <- character(length(target_region_keys))
target_region_labels <- character(length(target_region_keys))

for (i in seq_along(target_region_keys)) {
    target_key <- target_region_keys[i]
    used_region_ids <- target_region_ids[seq_len(max(0L, i - 1L))]

    exact_idx <- which(
        region_lookup$region_key == target_key &
            !region_lookup$region %in% used_region_ids
    )

    if (length(exact_idx) >= 1L) {
        target_region_ids[i] <- region_lookup$region[exact_idx[1L]]
        target_region_labels[i] <- region_lookup$region_label[exact_idx[1L]]
        next
    }

    partial_idx <- grep(target_key, region_lookup$region_key, fixed = TRUE)
    partial_idx <- partial_idx[
        !region_lookup$region[partial_idx] %in% used_region_ids
    ]

    if (length(partial_idx) >= 1L) {
        target_region_ids[i] <- region_lookup$region[partial_idx[1L]]
        target_region_labels[i] <- region_lookup$region_label[partial_idx[1L]]
        next
    }

    fuzzy_idx <- agrep(
        target_key,
        region_lookup$region_key,
        max.distance = 0.35
    )

    fuzzy_idx <- fuzzy_idx[
        !region_lookup$region[fuzzy_idx] %in% used_region_ids
    ]

    if (length(fuzzy_idx) >= 1L) {
        target_region_ids[i] <- region_lookup$region[fuzzy_idx[1L]]
        target_region_labels[i] <- region_lookup$region_label[fuzzy_idx[1L]]
        next
    }

    target_region_ids[i] <- NA_character_
    target_region_labels[i] <- NA_character_
}

if (anyNA(target_region_ids)) {
    target_region_ids <- head(region_names, 3L)

    target_region_labels <- region_lookup$region_label[match(
        target_region_ids,
        region_lookup$region
    )]
}
# ----------------------------------------------------------------------------
# Phase 10: temperature–mortality relative-risk curves by selected region and age.
# ----------------------------------------------------------------------------
selected_age_groups <- c(70L, 80L, 90L)
selected_age_labels <- c("70-74", "80-84", "90+")

mortality_temp_df <- DfM_all %>%
    dplyr::mutate(
        region = as.character(region),
        age = as.integer(age),
        ISODate = as.Date(ISODate),
        week_of_year = lubridate::isoweek(ISODate)
    ) %>%
    dplyr::left_join(
        climate_weekly_region %>%
            dplyr::transmute(
                region = as.character(Region),
                ISODate = as.Date(Week),
                temp_weekly = tg
            ),
        by = c("region", "ISODate")
    ) %>%
    dplyr::filter(
        region %in% target_region_ids,
        age %in% selected_age_groups,
        is.finite(temp_weekly),
        is.finite(deaths),
        is.finite(expo),
        expo > 0
    )

compute_rr_bundle <- function(fit, temp_grid, ref_temp, boot_draws) {
    newdata_grid <- data.frame(
        temp_weekly = temp_grid,
        week_of_year = 26,
        expo = 1
    )

    newdata_ref <- data.frame(
        temp_weekly = rep(ref_temp, length(temp_grid)),
        week_of_year = 26,
        expo = 1
    )

    design_grid <- stats::model.matrix(
        stats::delete.response(stats::terms(fit)),
        newdata_grid
    )

    design_ref <- stats::model.matrix(
        stats::delete.response(stats::terms(fit)),
        newdata_ref[1, , drop = FALSE]
    )

    beta_hat <- stats::coef(fit)

    eta_hat <- as.numeric(design_grid %*% beta_hat)
    eta_ref_hat <- as.numeric(design_ref %*% beta_hat)

    rr_hat <- exp(eta_hat - eta_ref_hat)

    vc <- tryCatch(
        stats::vcov(fit),
        error = function(e) NULL
    )

    if (is.null(vc) || any(!is.finite(vc))) {
        return(list(
            rr = rr_hat,
            rr_lo = rr_hat,
            rr_hi = rr_hat
        ))
    }

    vc <- as.matrix(vc)
    vc <- vc + diag(1e-8, nrow(vc))

    beta_draws <- tryCatch(
        MASS::mvrnorm(
            boot_draws,
            mu = beta_hat,
            Sigma = vc
        ),
        error = function(e) {
            matrix(
                rep(beta_hat, boot_draws),
                nrow = boot_draws,
                byrow = TRUE
            )
        }
    )

    rr_draws <- vapply(
        seq_len(boot_draws),
        function(b) {
            eta_draw <- as.numeric(
                design_grid %*% beta_draws[b, ]
            )

            eta_ref_draw <- as.numeric(
                design_ref %*% beta_draws[b, ]
            )

            exp(eta_draw - eta_ref_draw)
        },
        numeric(length(temp_grid))
    )

    list(
        rr = rr_hat,
        rr_lo = apply(
            rr_draws,
            1,
            stats::quantile,
            0.025,
            na.rm = TRUE
        ),
        rr_hi = apply(
            rr_draws,
            1,
            stats::quantile,
            0.975,
            na.rm = TRUE
        )
    )
}

rr_curve_list <- list()
rr_mmt_list <- list()

for (region_id in target_region_ids) {
    region_data <- mortality_temp_df %>%
        dplyr::filter(region == region_id)

    if (nrow(region_data) < 50L) {
        next
    }

    region_label <- region_lookup$region_label[match(
        region_id,
        region_lookup$region
    )]

    if (is.na(region_label) || !nzchar(region_label)) {
        region_label <- region_id
    }

    pooled_data <- region_data %>%
        dplyr::mutate(
            age_factor = factor(
                age,
                levels = selected_age_groups
            )
        )

    pooled_fit <- tryCatch(
        stats::glm(
            deaths ~ age_factor +
                splines::ns(temp_weekly, df = 4) +
                splines::ns(week_of_year, df = 6),
            family = stats::quasipoisson(),
            offset = log(expo),
            data = pooled_data
        ),
        error = function(e) NULL
    )

    if (is.null(pooled_fit)) {
        next
    }

    pooled_temp_grid <- seq(
        stats::quantile(region_data$temp_weekly, 0.02, na.rm = TRUE),
        stats::quantile(region_data$temp_weekly, 0.98, na.rm = TRUE),
        length.out = 200L
    )

    pooled_newdata <- data.frame(
        age_factor = factor(
            selected_age_groups[1],
            levels = selected_age_groups
        ),
        temp_weekly = pooled_temp_grid,
        week_of_year = 26,
        expo = 1
    )

    pooled_pred <- as.numeric(
        stats::predict(
            pooled_fit,
            newdata = pooled_newdata,
            type = "link"
        )
    )

    ref_temp <- pooled_temp_grid[which.min(pooled_pred)]

    rr_mmt_list[[length(rr_mmt_list) + 1L]] <- data.frame(
        region = region_id,
        region_label = region_label,
        mmt_temp = ref_temp,
        stringsAsFactors = FALSE
    )

    for (age_idx in seq_along(selected_age_groups)) {
        age_value <- selected_age_groups[age_idx]
        age_label <- selected_age_labels[age_idx]

        age_data <- region_data %>%
            dplyr::filter(age == age_value)

        if (nrow(age_data) < 30L) {
            next
        }

        age_fit <- tryCatch(
            stats::glm(
                deaths ~
                    splines::ns(temp_weekly, df = 4) +
                    splines::ns(week_of_year, df = 6),
                family = stats::quasipoisson(),
                offset = log(expo),
                data = age_data
            ),
            error = function(e) NULL
        )

        if (is.null(age_fit)) {
            next
        }

        temp_grid <- seq(
            stats::quantile(age_data$temp_weekly, 0.02, na.rm = TRUE),
            stats::quantile(age_data$temp_weekly, 0.98, na.rm = TRUE),
            length.out = 200L
        )

        rr_bundle <- compute_rr_bundle(
            age_fit,
            temp_grid,
            ref_temp,
            B_rr
        )

        rr_curve_list[[length(rr_curve_list) + 1L]] <- data.frame(
            region = region_id,
            region_label = region_label,
            age = age_value,
            age_label = age_label,
            temp_weekly = temp_grid,
            rr = rr_bundle$rr,
            rr_lo = rr_bundle$rr_lo,
            rr_hi = rr_bundle$rr_hi,
            mmt_temp = ref_temp,
            stringsAsFactors = FALSE
        )
    }
}

rr_plot_df <- dplyr::bind_rows(rr_curve_list)
rr_mmt_df <- dplyr::bind_rows(rr_mmt_list)

if (nrow(rr_plot_df) > 0L) {
    rr_plot_df <- rr_plot_df %>%
        dplyr::mutate(
            region_label = factor(
                region_label,
                levels = target_region_labels
            ),
            age_label = factor(
                age_label,
                levels = selected_age_labels
            )
        )

    rr_ylim_upper <- max(rr_plot_df$rr_hi, na.rm = TRUE)

    if (!is.finite(rr_ylim_upper)) {
        rr_ylim_upper <- 4
    } else {
        rr_ylim_upper <- min(rr_ylim_upper * 1.05, 8)
    }

    write.csv(
        rr_plot_df,
        file = file.path(
            phase10_plot_dir,
            "figure_10_relative_risk_phase10.csv"
        ),
        row.names = FALSE
    )

    saved_files <- character(0)

    for (region_lab in levels(rr_plot_df$region_label)) {
        region_df <- rr_plot_df %>%
            dplyr::filter(region_label == region_lab)

        if (nrow(region_df) == 0L) {
            next
        }

        # Plot: age-specific temperature response relative risk for each selected region.
        p <- ggplot2::ggplot(
            region_df,
            ggplot2::aes(
                x = temp_weekly,
                y = rr,
                color = age_label,
                fill = age_label
            )
        ) +
            ggplot2::geom_line(
                linewidth = 1.3
            ) +

            ggplot2::geom_hline(
                yintercept = 1,
                linetype = "dashed",
                linewidth = 0.6,
                color = "grey40"
            ) +

            ggplot2::scale_color_manual(
                values = c(
                    "70-74" = "#1b9e77",
                    "80-84" = "#d95f02",
                    "90+" = "#7570b3"
                )
            ) +
            ggplot2::scale_fill_manual(
                values = c(
                    "70-74" = "#1b9e77",
                    "80-84" = "#d95f02",
                    "90+" = "#7570b3"
                )
            ) +
            ggplot2::scale_y_continuous(
                breaks = scales::pretty_breaks(n = 6)
            ) +
            ggplot2::labs(
                title = paste0(
                    "Temperature–Mortality Relative Risk: ",
                    region_lab
                ),
                subtitle = NULL,
                x = "Weekly mean temperature",
                y = "Relative Risk",
                color = "Age group",
                fill = "Age group"
            ) +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    face = "bold",
                    size = 16,
                    hjust = 0.5
                ),
                plot.subtitle = ggplot2::element_text(
                    size = 11,
                    hjust = 0.5,
                    color = "grey30"
                ),
                legend.position = "top",
                legend.title = ggplot2::element_text(
                    face = "bold"
                ),
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_line(
                    color = "grey90"
                ),
                panel.grid.major.y = ggplot2::element_line(
                    color = "grey90"
                ),
                axis.title = ggplot2::element_text(
                    face = "bold"
                ),
                axis.text = ggplot2::element_text(
                    color = "black"
                ),
                plot.margin = ggplot2::margin(
                    15,
                    18,
                    15,
                    15
                )
            )

        safe_label <- gsub(
            "[ \\+]+",
            "_",
            as.character(region_lab)
        )

        safe_label <- gsub(
            "[^A-Za-z0-9_-]",
            "",
            safe_label
        )

        out_file <- file.path(
            phase10_plot_dir,
            sprintf(
                "figure_10_relative_risk_region_%s.png",
                safe_label
            )
        )

        ggplot2::ggsave(
            out_file,
            plot = p,
            width = 10,
            height = 6.5,
            dpi = 400
        )
    }
}

# ============================================================================
# Plot: temperature response curves for CNN-LSTM and GNN-LSTM
# ============================================================================

# Plot data: build age-90 temperature response series
compute_temp_response <- function(test_pred_df, model_label) {
    temp_response_df <- test_pred_df %>%
        dplyr::filter(Age == "90") %>%
        dplyr::left_join(
            climate_weekly_region %>%
                dplyr::select(Week, Region, tg),
            by = c("Week", "Region")
        ) %>%
        dplyr::left_join(
            residual_long %>%
                dplyr::rename(Age = AgeNumeric, observed_residual = residual),
            by = c("Week", "Region", "Age")
        ) %>%
        dplyr::filter(
            is.finite(tg),
            is.finite(residual_pred),
            is.finite(observed_residual)
        ) %>%
        dplyr::mutate(
            effect_pred = residual_pred - median(residual_pred, na.rm = TRUE),
            effect_obs = observed_residual -
                median(observed_residual, na.rm = TRUE),
            model = model_label
        )

    return(temp_response_df)
}

# Plot data: bin temperature response curves for smoother lines
create_temp_summary <- function(temp_response_df, model_label) {
    temp_summary <- temp_response_df %>%
        dplyr::mutate(
            temp_bin = cut(
                tg,
                breaks = seq(
                    min(tg, na.rm = TRUE),
                    max(tg, na.rm = TRUE),
                    length.out = 45
                ),
                include.lowest = TRUE
            )
        ) %>%
        dplyr::group_by(temp_bin) %>%
        dplyr::summarise(
            tg = mean(tg, na.rm = TRUE),
            effect_pred = mean(effect_pred, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::arrange(tg) %>%
        dplyr::mutate(model = model_label)

    return(temp_summary)
}

# Compute temperature response for both models
cnn_temp_response_df <- compute_temp_response(test_pred_df_cnn, "CNN-LSTM")
gnn_temp_response_df <- compute_temp_response(test_pred_df_gnn, "GNN-LSTM")

# Create summaries
cnn_temp_summary <- create_temp_summary(cnn_temp_response_df, "CNN-LSTM")
gnn_temp_summary <- create_temp_summary(gnn_temp_response_df, "GNN-LSTM")

# Create median observed value summaries for overlay
cnn_obs_summary <- cnn_temp_response_df %>%
    dplyr::mutate(
        temp_bin = cut(
            tg,
            breaks = seq(
                min(tg, na.rm = TRUE),
                max(tg, na.rm = TRUE),
                length.out = 45
            ),
            include.lowest = TRUE
        )
    ) %>%
    dplyr::group_by(temp_bin) %>%
    dplyr::summarise(
        tg = mean(tg, na.rm = TRUE),
        effect_obs = mean(effect_obs, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::arrange(tg)

gnn_obs_summary <- gnn_temp_response_df %>%
    dplyr::mutate(
        temp_bin = cut(
            tg,
            breaks = seq(
                min(tg, na.rm = TRUE),
                max(tg, na.rm = TRUE),
                length.out = 45
            ),
            include.lowest = TRUE
        )
    ) %>%
    dplyr::group_by(temp_bin) %>%
    dplyr::summarise(
        tg = mean(tg, na.rm = TRUE),
        effect_obs = mean(effect_obs, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::arrange(tg)

# Find common y-axis range (including both predicted and observed values)
y_min_common <- min(
    c(
        cnn_temp_summary$effect_pred,
        gnn_temp_summary$effect_pred,
        cnn_obs_summary$effect_obs,
        gnn_obs_summary$effect_obs
    ),
    na.rm = TRUE
)
y_max_common <- max(
    c(
        cnn_temp_summary$effect_pred,
        gnn_temp_summary$effect_pred,
        cnn_obs_summary$effect_obs,
        gnn_obs_summary$effect_obs
    ),
    na.rm = TRUE
)
y_margin <- (y_max_common - y_min_common) * 0.1
y_limits <- c(y_min_common - y_margin, y_max_common + y_margin)

# Create individual plots with synchronized y-axis
if (nrow(cnn_temp_response_df) > 0L) {
    temp_quantiles_cnn <- stats::quantile(
        cnn_temp_response_df$tg,
        probs = c(0.05, 0.95),
        na.rm = TRUE
    )

    cnn_temp_response_plot <- ggplot2::ggplot(
        cnn_temp_summary,
        ggplot2::aes(x = tg, y = effect_pred)
    ) +
        ggplot2::geom_line(
            data = cnn_obs_summary,
            inherit.aes = FALSE,
            ggplot2::aes(x = tg, y = effect_obs),
            color = "#CCCCCC",
            linewidth = 1.0,
            alpha = 0.7
        ) +
        ggplot2::geom_line(color = "#e08214", linewidth = 1.2) +
        ggplot2::geom_vline(
            xintercept = temp_quantiles_cnn,
            linetype = "dashed",
            color = "grey40",
            linewidth = 0.5
        ) +
        ggplot2::coord_cartesian(ylim = y_limits) +
        ggplot2::labs(
            title = "CNN-LSTM Predicted Excess Mortality vs Weekly Mean Temperature",
            x = "Weekly Mean Temperature",
            y = "Excess Mortality"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            legend.position = "none",
            panel.grid.minor = ggplot2::element_blank()
        )

    cnn_temp_plot_dir <- file.path(
        "Results",
        "plots",
        "cnn_temperature_response"
    )
    dir.create(cnn_temp_plot_dir, recursive = TRUE, showWarnings = FALSE)

    cnn_temp_plot_file <- file.path(
        cnn_temp_plot_dir,
        "figure_cnn_lstm_temperature_response_age_90.png"
    )

    ggplot2::ggsave(
        filename = cnn_temp_plot_file,
        plot = cnn_temp_response_plot,
        width = 10,
        height = 6,
        dpi = 300,
        bg = "white"
    )
}

if (nrow(gnn_temp_response_df) > 0L) {
    temp_quantiles_gnn <- stats::quantile(
        gnn_temp_response_df$tg,
        probs = c(0.05, 0.95),
        na.rm = TRUE
    )

    gnn_temp_response_plot <- ggplot2::ggplot(
        gnn_temp_summary,
        ggplot2::aes(x = tg, y = effect_pred)
    ) +
        ggplot2::geom_line(
            data = gnn_obs_summary,
            inherit.aes = FALSE,
            ggplot2::aes(x = tg, y = effect_obs),
            color = "#CCCCCC",
            linewidth = 1.0,
            alpha = 0.7
        ) +
        ggplot2::geom_line(color = "#8b5fbf", linewidth = 1.2) +
        ggplot2::geom_vline(
            xintercept = temp_quantiles_gnn,
            linetype = "dashed",
            color = "grey40",
            linewidth = 0.5
        ) +
        ggplot2::coord_cartesian(ylim = y_limits) +
        ggplot2::labs(
            title = "GNN-LSTM Predicted Excess Mortality vs Weekly Mean Temperature",
            x = "Weekly Mean Temperature",
            y = "Excess Mortality"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14),
            legend.position = "none",
            panel.grid.minor = ggplot2::element_blank()
        )

    gnn_temp_plot_dir <- file.path(
        "Results",
        "plots",
        "gnn_temperature_response"
    )
    dir.create(gnn_temp_plot_dir, recursive = TRUE, showWarnings = FALSE)

    gnn_temp_plot_file <- file.path(
        gnn_temp_plot_dir,
        "figure_gnn_lstm_temperature_response_age_90.png"
    )

    ggplot2::ggsave(
        filename = gnn_temp_plot_file,
        plot = gnn_temp_response_plot,
        width = 10,
        height = 6,
        dpi = 300,
        bg = "white"
    )
}

# ============================================================================
# Plot: regional temperature-response overlays for CNN-LSTM vs GNN-LSTM
# ============================================================================

# Get unique regions from data
regions <- sort(unique(c(
    cnn_temp_response_df$Region,
    gnn_temp_response_df$Region
)))

# Create model colors mapping
model_colors <- c("CNN-LSTM" = "#e08214", "GNN-LSTM" = "#8b5fbf")

# Create directory for regional overlays
overlay_plot_dir <- file.path(
    "Results",
    "plots",
    "temperature_response_overlay"
)
dir.create(overlay_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Create overlay plot for each region
for (region_name in regions) {
    cnn_region_raw <- cnn_temp_response_df %>%
        dplyr::filter(Region == region_name)

    cnn_region_data <- cnn_region_raw %>%
        dplyr::mutate(
            temp_bin = cut(
                tg,
                breaks = seq(
                    min(tg, na.rm = TRUE),
                    max(tg, na.rm = TRUE),
                    length.out = 45
                ),
                include.lowest = TRUE
            )
        ) %>%
        dplyr::group_by(temp_bin) %>%
        dplyr::summarise(
            tg = mean(tg, na.rm = TRUE),
            effect_pred = mean(effect_pred, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::arrange(tg) %>%
        dplyr::mutate(model = "CNN-LSTM")

    gnn_region_raw <- gnn_temp_response_df %>%
        dplyr::filter(Region == region_name)

    gnn_region_data <- gnn_region_raw %>%
        dplyr::mutate(
            temp_bin = cut(
                tg,
                breaks = seq(
                    min(tg, na.rm = TRUE),
                    max(tg, na.rm = TRUE),
                    length.out = 45
                ),
                include.lowest = TRUE
            )
        ) %>%
        dplyr::group_by(temp_bin) %>%
        dplyr::summarise(
            tg = mean(tg, na.rm = TRUE),
            effect_pred = mean(effect_pred, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::arrange(tg) %>%
        dplyr::mutate(model = "GNN-LSTM")

    region_combined <- dplyr::bind_rows(cnn_region_data, gnn_region_data)
    region_observed <- dplyr::bind_rows(cnn_region_raw, gnn_region_raw) %>%
        dplyr::mutate(
            temp_bin = cut(
                tg,
                breaks = seq(
                    min(tg, na.rm = TRUE),
                    max(tg, na.rm = TRUE),
                    length.out = 45
                ),
                include.lowest = TRUE
            )
        ) %>%
        dplyr::group_by(temp_bin) %>%
        dplyr::summarise(
            tg = mean(tg, na.rm = TRUE),
            effect_obs = mean(effect_obs, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::arrange(tg)

    if (nrow(region_combined) > 0L) {
        region_overlay_plot <- ggplot2::ggplot(
            region_combined,
            ggplot2::aes(x = tg, y = effect_pred, color = model)
        ) +
            ggplot2::geom_line(
                data = region_observed,
                inherit.aes = FALSE,
                ggplot2::aes(x = tg, y = effect_obs),
                color = "#CCCCCC",
                linewidth = 1.0,
                alpha = 0.7
            ) +
            ggplot2::geom_line(linewidth = 1.2, linetype = "solid") +
            ggplot2::scale_color_manual(values = model_colors) +
            ggplot2::coord_cartesian(ylim = y_limits) +
            ggplot2::labs(
                title = paste(
                    "Temperature Response Comparison -",
                    region_name,
                    "(Age 90+)"
                ),
                x = "Weekly Mean Temperature",
                y = "Excess Mortality",
                color = "Model"
            ) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 14),
                legend.position = "top",
                legend.direction = "horizontal",
                legend.title = ggplot2::element_text(size = 10, face = "bold"),
                legend.text = ggplot2::element_text(size = 10),
                panel.grid.minor = ggplot2::element_blank()
            )

        overlay_plot_file <- file.path(
            overlay_plot_dir,
            paste0(
                "temperature_response_overlay_",
                tolower(region_name),
                ".png"
            )
        )

        ggplot2::ggsave(
            filename = overlay_plot_file,
            plot = region_overlay_plot,
            width = 10,
            height = 6,
            dpi = 300,
            bg = "white"
        )
    }
}
