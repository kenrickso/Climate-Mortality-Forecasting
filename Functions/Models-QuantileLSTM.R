# ------------------------------------------------------------------------------
# Function: cnn_mort_fc_net
# Description: CNN-GRU/LSTM framework for forecasting mortality residuals, incorporating temporal, age, and regional embeddings.
# Inputs:
#   - seq_length: Length of the input sequence.
#   - input_size: Number of features per timestep.
#   - cnn_filters: Vector of CNN filter sizes.
#   - hidden_size_1: Units for the first dense layer.
#   - hidden_size_2: Units for the second dense layer.
#   - hidden_size_3: Units for the third dense layer.
#   - output_size: Size of the output projection.
#   - cell_type: Type of recurrent cell ("lstm" or "gru").
#   - use_cnn: Logical to enable/disable CNN layers.
#   - use_dense: Logical to enable/disable dense layers.
#   - use_age_embedding: Logical to use post-LSTM age embeddings.
#   - use_age_embedding_pre_lstm: Logical to use pre-LSTM age embeddings.
#   - use_region_embedding: Logical to use region embeddings.
#   - n_regions: Total number of regions (if using region embeddings).
#   - region_embedding_dim: Dimension for region embeddings.
#   - zero_embeddings: Logical to zero out embeddings if testing without them.
#   - lstm_units: Explicit number of LSTM units (overrides auto-derived value).
#   - dropout_rate: Configurable dropout rate.
# Output:
#   - A compiled Keras model predicting mortality residuals.
# ------------------------------------------------------------------------------
cnn_mort_fc_net <- function(
    seq_length,
    input_size,
    cnn_filters = c(64, 32, 16),
    hidden_size_1 = 80,
    hidden_size_2 = 40,
    hidden_size_3 = 20,
    output_size = 1,
    cell_type = c("lstm", "gru"),
    use_cnn = TRUE,
    use_dense = TRUE,
    use_age_embedding = FALSE,
    # New options for pre-LSTM embeddings
    use_age_embedding_pre_lstm = FALSE,
    use_region_embedding = FALSE,
    n_regions = NULL,
    region_embedding_dim = 8,
    zero_embeddings = FALSE,
    # Explicit LSTM units (overrides default auto-derived value)
    lstm_units = NULL,
    # Configurable dropout
    dropout_rate = 0.05
) {
    cell_type <- match.arg(cell_type)
    # Primary temporal input: (seq_length, input_size)
    inputs <- layer_input(
        shape = c(seq_length, input_size),
        name = "temporal_input"
    )

    x <- inputs
    if (use_cnn) {
        # CNN layers for dimension reduction
        for (f in cnn_filters) {
            x <- x %>%
                layer_conv_1d(
                    filters = f,
                    kernel_size = 5,
                    padding = 'same',
                    activation = 'relu'
                ) %>%
                layer_batch_normalization() %>%
                layer_spatial_dropout_1d(rate = dropout_rate)
        }
        # feature size per timestep after CNN
        feat_size <- as.integer(tail(cnn_filters, 1))
        default_lstm_units <- if (use_dense) feat_size else output_size
    } else {
        # Skip CNN layers
        feat_size <- as.integer(input_size)
        default_lstm_units <- if (use_dense) feat_size else output_size
    }

    # Allow explicit override of LSTM units
    lstm_units_local <- if (!is.null(lstm_units)) {
        as.integer(lstm_units)
    } else {
        as.integer(default_lstm_units)
    }

    # Optional region embedding (per-sample, tiled across timesteps and concatenated)
    region_input <- NULL
    if (isTRUE(use_region_embedding)) {
        if (is.null(n_regions)) {
            stop("n_regions must be provided when use_region_embedding = TRUE")
        }
        region_input <- layer_input(
            shape = c(1),
            dtype = 'int32',
            name = "region_input"
        )
        region_emb <- region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions) + 1L,
                output_dim = as.integer(max(1, region_embedding_dim)),
                input_length = 1,
                name = "region_embedding"
            ) %>%
            layer_flatten()
        # Optionally zero-out embeddings
        if (isTRUE(zero_embeddings)) {
            region_emb <- region_emb %>%
                layer_lambda(
                    f = function(x) x * 0,
                    output_shape = function(sh) sh
                )
        }
        # Repeat to match timesteps: (batch, seq_length, region_embedding_dim)
        region_emb_seq <- region_emb %>%
            layer_repeat_vector(as.integer(seq_length))
        feat_size <- feat_size + as.integer(region_embedding_dim)
    }

    # Optional age embedding tiled per timestep (pre-LSTM)
    age_tiled_layer <- NULL
    if (
        isTRUE(use_age_embedding_pre_lstm) &&
            exists("age_embeddings", envir = .GlobalEnv)
    ) {
        age_dim <- as.integer(ncol(age_embeddings))
        if (isTRUE(zero_embeddings)) {
            age_tiled_layer <- layer_lambda(
                f = function(seq_x) {
                    batch_size <- tensorflow::tf$shape(seq_x)[1]
                    tensorflow::tf$zeros(
                        shape = c(
                            batch_size,
                            as.integer(seq_length),
                            as.integer(age_dim)
                        )
                    )
                },
                output_shape = c(as.integer(seq_length), as.integer(age_dim))
            )
        } else {
            age_tiled_layer <- layer_lambda(
                f = function(seq_x) {
                    batch_size <- tensorflow::tf$shape(seq_x)[1]
                    # age_embeddings assumed to be a 1 x age_dim numeric matrix in global env
                    age_embed_const <- tensorflow::tf$constant(
                        age_embeddings,
                        dtype = tensorflow::tf$float32
                    )
                    # reshape to (1,1,age_dim) then tile to (batch, seq_length, age_dim)
                    age_embeds_reshaped <- tensorflow::tf$reshape(
                        age_embed_const,
                        shape = c(1L, 1L, as.integer(age_dim))
                    )
                    tensorflow::tf$tile(
                        age_embeds_reshaped,
                        c(batch_size, as.integer(seq_length), 1L)
                    )
                },
                output_shape = c(as.integer(seq_length), as.integer(age_dim))
            )
        }
        feat_size <- feat_size + as.integer(age_dim)
    }

    # If either region or age pre-LSTM embeddings exist, concatenate them to timestep features
    if (!is.null(region_input) || !is.null(age_tiled_layer)) {
        concat_list <- list(x)
        if (!is.null(region_input)) {
            concat_list <- c(concat_list, list(region_emb_seq))
        }
        if (!is.null(age_tiled_layer)) {
            concat_list <- c(concat_list, list(age_tiled_layer(x)))
        }
        # Note: if age_tiled_layer is present we call it with x to get a tiled tensor
        # Use layer_concatenate to merge across the feature axis
        x <- layer_concatenate(concat_list, axis = 2L)
    }

    # Recurrent layer: choose LSTM or GRU by `cell_type` (processes combined timestep features)
    if (cell_type == "gru") {
        x <- x %>%
            layer_gru(units = lstm_units_local, return_sequences = FALSE)
    } else {
        x <- x %>%
            layer_lstm(units = lstm_units_local, return_sequences = FALSE)
    }

    # Add age embeddings after LSTM if requested (old behaviour)
    if (
        use_age_embedding &&
            use_dense &&
            exists("age_embeddings", envir = .GlobalEnv)
    ) {
        # Concatenate standardized age values with temporal features
        concat_output_size <- as.integer(lstm_units_local + output_size)

        add_age_embeddings <- layer_lambda(
            f = function(temporal_features) {
                batch_size <- tensorflow::tf$shape(temporal_features)[1]
                # Get age embeddings from global environment
                age_embed_const <- tensorflow::tf$constant(
                    age_embeddings,
                    dtype = tensorflow::tf$float32
                )
                # Tile to match batch size
                age_embeds_tiled <- tensorflow::tf$tile(
                    age_embed_const,
                    c(batch_size, 1L)
                )
                # Concatenate
                tensorflow::tf$concat(
                    list(temporal_features, age_embeds_tiled),
                    axis = 1L
                )
            },
            output_shape = list(concat_output_size)
        )

        x <- add_age_embeddings(x)
    }

    if (use_dense) {
        # Simplified dense head to reduce gradient dampening
        # FC1 -> LeakyReLU -> Dropout (removed LayerNorm to improve gradient flow)
        x <- x %>%
            layer_dense(
                units = hidden_size_1,
                use_bias = TRUE
            ) %>%
            layer_activation_leaky_relu(alpha = 0.1) %>%
            layer_dropout(rate = dropout_rate)

        # FC2 -> LeakyReLU -> Dropout
        x <- x %>%
            layer_dense(
                units = hidden_size_2,
                use_bias = TRUE
            ) %>%
            layer_activation_leaky_relu(alpha = 0.1) %>%
            layer_dropout(rate = dropout_rate)

        # FC3 -> LeakyReLU -> Dropout
        x <- x %>%
            layer_dense(
                units = hidden_size_3,
                use_bias = TRUE
            ) %>%
            layer_activation_leaky_relu(alpha = 0.1) %>%
            layer_dropout(rate = dropout_rate)

        # FC4 (output) - linear activation for regression
        outputs <- x %>%
            layer_dense(units = output_size, use_bias = TRUE)
    } else {
        # No dense layers, output is the recurrent output
        outputs <- x
    }

    # Build model: include region input if used
    model_inputs <- if (!is.null(region_input)) {
        list(inputs, region_input)
    } else {
        inputs
    }

    model <- keras_model(inputs = model_inputs, outputs = outputs)
    return(model)
}

# ------------------------------------------------------------------------------
# Function: gnn_lstm_scalar_output
# Description: A scalar-output GNN-LSTM model that takes spatial climate sequences, age index, and region index to output a single mortality residual.
# Inputs:
#   - seq_length: Sequence length for temporal input.
#   - n_regions: Number of regions in the spatial graph.
#   - features_per_region: Number of climate features per region.
#   - adjacency_matrix: Spatial adjacency matrix for the GNN.
#   - n_ages: Number of age groups.
#   - gnn_hidden_units: Vector defining GNN hidden layer dimensions.
#   - gnn_activation: Activation function for GNN layers.
#   - lstm_units: Number of units in the LSTM layer.
#   - hidden_size_1: Units for the first dense layer.
#   - hidden_size_2: Units for the second dense layer.
#   - hidden_size_3: Units for the third dense layer.
#   - age_embedding_dim: Dimension for the age embedding.
#   - region_embedding_dim: Dimension for the region embedding.
#   - dropout_rate: Dropout rate for regularization.
# Output:
#   - A compiled Keras model outputting a single scalar residual.
# ------------------------------------------------------------------------------
gnn_lstm_scalar_output_region_preserving <- function(
    seq_length,
    n_regions,
    features_per_region,
    adjacency_matrix,
    n_ages,
    gnn_hidden_units = c(64, 32, 16),
    gnn_activation = "relu",
    lstm_units = 64,
    hidden_size_1 = 64,
    hidden_size_2 = 32,
    hidden_size_3 = 16,
    age_embedding_dim = 8,
    region_embedding_dim = 8,
    dropout_rate = 0.1,
    use_age_embedding = TRUE,
    use_region_embedding = TRUE,
    use_age_numeric = FALSE,
    use_year_numeric = FALSE,
    age_numeric_dim = NULL
) {
    climate_input <- layer_input(
        shape = c(seq_length, n_regions * features_per_region),
        name = "climate_input"
    )
    age_input <- layer_input(shape = c(1), dtype = "int32", name = "age_input")
    region_input <- layer_input(
        shape = c(1),
        dtype = "int32",
        name = "region_input"
    )

    adj_tf <- tensorflow::tf$constant(
        adjacency_matrix,
        dtype = tensorflow::tf$float32
    )

    gnn_output_size <- as.integer(tail(gnn_hidden_units, 1))
    seq_len_int <- as.integer(seq_length)
    n_regions_int <- as.integer(n_regions)
    features_per_region_int <- as.integer(features_per_region)
    combined_region_feature_size <- as.integer(2L * gnn_output_size)

    gnn_layers <- list()
    for (i in seq_along(gnn_hidden_units)) {
        gnn_layers[[i]] <- layer_dense(
            units = as.integer(gnn_hidden_units[i]),
            activation = gnn_activation,
            name = paste0("gnn_scalar_keep_regions_layer_", i),
            kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
        )
    }

    # Share one age embedding across the graph block and the scalar head so
    # the same categorical table is not trained twice.
    if (isTRUE(use_age_embedding)) {
        age_emb <- age_input %>%
            layer_embedding(
                input_dim = as.integer(n_ages + 1),
                output_dim = as.integer(age_embedding_dim),
                name = "gnn_scalar_age_embedding"
            ) %>%
            layer_flatten()
    } else {
        age_emb <- age_input %>%
            layer_embedding(
                input_dim = as.integer(n_ages + 1),
                output_dim = as.integer(age_embedding_dim),
                name = "gnn_scalar_age_embedding",
                embeddings_initializer = "zeros",
                trainable = FALSE
            ) %>%
            layer_flatten()
    }

    age_graph_emb <- age_emb %>%
        layer_dense(
            units = as.integer(gnn_output_size),
            activation = "relu",
            name = "gnn_graph_age_projection",
            kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
        )

    region_graph_emb <- if (isTRUE(use_region_embedding)) {
        region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(region_embedding_dim),
                name = "gnn_graph_region_embedding"
            ) %>%
            layer_flatten() %>%
            layer_dense(
                units = as.integer(gnn_output_size),
                activation = "relu",
                name = "gnn_graph_region_projection",
                kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
            )
    } else {
        region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(region_embedding_dim),
                name = "gnn_graph_region_embedding",
                embeddings_initializer = "zeros",
                trainable = FALSE
            ) %>%
            layer_flatten() %>%
            layer_dense(
                units = as.integer(gnn_output_size),
                activation = "relu",
                name = "gnn_graph_region_projection",
                kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
            )
    }

    graph_context <- layer_concatenate(
        list(age_graph_emb, region_graph_emb),
        axis = -1
    ) %>%
        layer_dense(
            units = as.integer(features_per_region_int),
            activation = "relu",
            name = "gnn_graph_context_projection",
            kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
        )

    graph_context_tiled <- graph_context %>%
        layer_repeat_vector(seq_len_int)

    # Optional scalar numeric age input (midpoint like 67.5). When provided,
    # project via a non-negative dense layer, tile across timesteps, and
    # include both a tiled version for pre-LSTM concatenation and a flat
    # version for the post-LSTM skip connection.
    age_numeric_input <- NULL
    age_numeric_proj <- NULL
    age_numeric_tiled <- NULL
    age_numeric_flat <- NULL
    if (isTRUE(use_age_numeric)) {
        age_numeric_input <- layer_input(
            shape = c(1),
            dtype = 'float32',
            name = 'age_numeric_input'
        )
        if (is.null(age_numeric_dim)) {
            age_numeric_dim_local <- as.integer(age_embedding_dim)
        } else {
            age_numeric_dim_local <- as.integer(age_numeric_dim)
        }
        age_numeric_proj <- age_numeric_input %>%
            layer_dense(
                units = age_numeric_dim_local,
                activation = 'relu',
                name = 'gnn_age_numeric_projection'
            )
        age_numeric_tiled <- age_numeric_proj %>%
            layer_repeat_vector(seq_len_int)
        age_numeric_flat <- age_numeric_input %>%
            layer_dense(
                units = age_numeric_dim_local,
                activation = 'relu',
                name = 'gnn_age_numeric_flat'
            ) %>%
            layer_flatten()
    }

    year_numeric_input <- NULL
    year_numeric_tiled <- NULL
    year_numeric_flat <- NULL
    if (isTRUE(use_year_numeric)) {
        year_numeric_input <- layer_input(
            shape = c(1),
            dtype = 'float32',
            name = 'year_numeric_input'
        )
        year_numeric_tiled <- year_numeric_input %>%
            layer_repeat_vector(seq_len_int)
        year_numeric_flat <- year_numeric_input %>%
            layer_flatten()
    }

    # Keep the entire message-passing path inside a custom Keras layer so the
    # Dense sublayers are registered as part of the graph instead of being
    # captured through a Lambda closure.
    GNNMessagePassingLayer <- new_layer_class(
        "GNNMessagePassingLayer",
        initialize = function(
            gnn_layers,
            adj_tf,
            seq_len_int,
            n_regions_int,
            features_per_region_int,
            gnn_output_size
        ) {
            super$initialize()
            self$gnn_layers <- gnn_layers
            self$adj_tf <- adj_tf
            self$seq_len_int <- as.integer(seq_len_int)
            self$n_regions_int <- as.integer(n_regions_int)
            self$features_per_region_int <- as.integer(features_per_region_int)
            self$gnn_output_size <- as.integer(gnn_output_size)
        },
        call = function(inputs) {
            x <- inputs[[1]]
            graph_context_seq <- inputs[[2]]

            batch_size <- tensorflow::tf$shape(x)[1]

            x_flat <- tensorflow::tf$reshape(
                x,
                list(-1L, self$n_regions_int, self$features_per_region_int)
            )
            ctx_flat <- tensorflow::tf$reshape(
                graph_context_seq,
                list(-1L, 1L, self$features_per_region_int)
            )
            ctx_tiled <- tensorflow::tf$tile(
                ctx_flat,
                c(1L, self$n_regions_int, 1L)
            )

            h <- tensorflow::tf$concat(list(x_flat, ctx_tiled), axis = -1L)
            for (gnn_layer in self$gnn_layers) {
                h_agg <- tensorflow::tf$matmul(self$adj_tf, h)
                h <- gnn_layer(h_agg)
            }

            tensorflow::tf$reshape(
                h,
                tensorflow::tf$stack(list(
                    batch_size,
                    self$seq_len_int,
                    self$n_regions_int,
                    self$gnn_output_size
                ))
            )
        }
    )

    x_gnn_nodes <- GNNMessagePassingLayer(
        gnn_layers = gnn_layers,
        adj_tf = adj_tf,
        seq_len_int = seq_len_int,
        n_regions_int = n_regions_int,
        features_per_region_int = features_per_region_int,
        gnn_output_size = gnn_output_size
    )(list(climate_input, graph_context_tiled))

    # Learned attention over regions preserves the spatial structure until the
    # model has produced a graph-aware regional sequence.
    region_attention_logits <- x_gnn_nodes %>%
        layer_reshape(
            target_shape = c(seq_len_int * n_regions_int, gnn_output_size)
        ) %>%
        layer_dense(
            units = 1L,
            activation = NULL,
            name = "gnn_region_attention_logits",
            kernel_regularizer = keras$regularizers$L2(l2 = 1e-5)
        ) %>%
        layer_reshape(
            target_shape = c(seq_len_int, n_regions_int, 1L)
        )

    x_gnn_with_attention <- layer_concatenate(
        list(x_gnn_nodes, region_attention_logits),
        axis = -1
    )

    TargetAttentionExtractionLayer <- new_layer_class(
        "TargetAttentionExtractionLayer",
        initialize = function(seq_len_int, n_regions_int, gnn_output_size) {
            super$initialize()
            self$seq_len_int <- as.integer(seq_len_int)
            self$n_regions_int <- as.integer(n_regions_int)
            self$gnn_output_size <- as.integer(gnn_output_size)
        },
        call = function(inputs) {
            gnn_tensor <- inputs[[1]]
            region_tensor <- inputs[[2]]

            region_idx <- tensorflow::tf$cast(
                region_tensor,
                tensorflow::tf$int32
            )
            region_idx <- tensorflow::tf$squeeze(region_idx, axis = 1L) - 1L

            split_parts <- tensorflow::tf$split(
                gnn_tensor,
                num_or_size_splits = c(self$gnn_output_size, 1L),
                axis = -1L
            )
            node_tensor <- split_parts[[1]]
            attention_logits <- split_parts[[2]]

            attention_weights <- tensorflow::tf$nn$softmax(
                attention_logits,
                axis = 2L
            )
            global_region_seq <- tensorflow::tf$reduce_sum(
                node_tensor * attention_weights,
                axis = 2L
            )

            region_one_hot <- tensorflow::tf$one_hot(
                region_idx,
                depth = self$n_regions_int,
                dtype = tensorflow::tf$float32
            )
            region_one_hot <- tensorflow::tf$expand_dims(
                region_one_hot,
                axis = 1L
            )
            region_one_hot <- tensorflow::tf$expand_dims(
                region_one_hot,
                axis = -1L
            )

            target_region_seq <- tensorflow::tf$reduce_sum(
                node_tensor * region_one_hot,
                axis = 2L
            )

            tensorflow::tf$concat(
                list(target_region_seq, global_region_seq),
                axis = -1L
            )
        }
    )

    FocalRegionRawExtractionLayer <- new_layer_class(
        "FocalRegionRawExtractionLayer",
        initialize = function(
            seq_len_int,
            n_regions_int,
            features_per_region_int
        ) {
            super$initialize()
            self$seq_len_int <- as.integer(seq_len_int)
            self$n_regions_int <- as.integer(n_regions_int)
            self$features_per_region_int <- as.integer(features_per_region_int)
        },
        call = function(inputs) {
            climate_tensor <- inputs[[1]]
            region_tensor <- inputs[[2]]

            batch_size <- tensorflow::tf$shape(climate_tensor)[1]
            climate_tensor_4d <- tensorflow::tf$reshape(
                climate_tensor,
                tensorflow::tf$stack(list(
                    batch_size,
                    self$seq_len_int,
                    self$n_regions_int,
                    self$features_per_region_int
                ))
            )

            region_idx <- tensorflow::tf$cast(
                region_tensor,
                tensorflow::tf$int32
            )
            region_idx <- tensorflow::tf$squeeze(region_idx, axis = 1L) - 1L
            region_one_hot <- tensorflow::tf$one_hot(
                region_idx,
                depth = self$n_regions_int,
                dtype = tensorflow::tf$float32
            )
            region_one_hot <- tensorflow::tf$expand_dims(
                region_one_hot,
                axis = 1L
            )
            region_one_hot <- tensorflow::tf$expand_dims(
                region_one_hot,
                axis = -1L
            )

            tensorflow::tf$reduce_sum(
                climate_tensor_4d * region_one_hot,
                axis = 2L
            )
        }
    )

    x_gnn_seq <- TargetAttentionExtractionLayer(
        seq_len_int = seq_len_int,
        n_regions_int = n_regions_int,
        gnn_output_size = gnn_output_size
    )(list(x_gnn_with_attention, region_input))

    focal_region_raw_seq <- FocalRegionRawExtractionLayer(
        seq_len_int = seq_len_int,
        n_regions_int = n_regions_int,
        features_per_region_int = features_per_region_int
    )(list(climate_input, region_input))

    region_emb <- if (isTRUE(use_region_embedding)) {
        region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(region_embedding_dim),
                name = "gnn_scalar_region_embedding"
            ) %>%
            layer_flatten()
    } else {
        region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(region_embedding_dim),
                name = "gnn_scalar_region_embedding",
                embeddings_initializer = "zeros",
                trainable = FALSE
            ) %>%
            layer_flatten()
    }

    age_emb_tiled <- age_emb %>% layer_repeat_vector(seq_len_int)
    region_emb_tiled <- region_emb %>% layer_repeat_vector(seq_len_int)

    concat_inputs <- list(
        x_gnn_seq,
        focal_region_raw_seq,
        age_emb_tiled,
        region_emb_tiled
    )
    if (!is.null(age_numeric_tiled)) {
        concat_inputs <- c(concat_inputs, list(age_numeric_tiled))
    }
    if (!is.null(year_numeric_tiled)) {
        concat_inputs <- c(concat_inputs, list(year_numeric_tiled))
    }
    x <- layer_concatenate(
        concat_inputs,
        axis = -1
    ) %>%
        layer_layer_normalization()

    x <- x %>%
        layer_lstm(units = as.integer(lstm_units), return_sequences = FALSE) %>%
        layer_dropout(rate = dropout_rate)

    final_skip <- list(x, age_emb, region_emb)
    if (!is.null(age_numeric_flat)) {
        final_skip <- c(final_skip, list(age_numeric_flat))
    }
    if (!is.null(year_numeric_flat)) {
        final_skip <- c(final_skip, list(year_numeric_flat))
    }
    x <- layer_concatenate(final_skip, axis = -1)

    x <- x %>%
        layer_dense(
            units = as.integer(hidden_size_1),
            activation = "relu"
        ) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(
            units = as.integer(hidden_size_2),
            activation = "relu"
        ) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(
            units = as.integer(hidden_size_3),
            activation = "relu"
        ) %>%
        layer_dropout(rate = dropout_rate)

    output <- x %>% layer_dense(units = 1L, name = "residual_output")

    model_inputs <- list(climate_input, age_input, region_input)
    if (!is.null(age_numeric_input)) {
        model_inputs <- c(model_inputs, list(age_numeric_input))
    }
    if (!is.null(year_numeric_input)) {
        model_inputs <- c(model_inputs, list(year_numeric_input))
    }

    keras_model(
        inputs = model_inputs,
        outputs = output
    )
}

# ------------------------------------------------------------------------------
# Function: gnn_lstm_scalar_output_quantile
# Description: Quantile-regression variant of gnn_lstm_scalar_output_region_preserving.
gnn_lstm_scalar_output_quantile <- function(
    seq_length,
    n_regions,
    features_per_region,
    adjacency_matrix,
    n_ages,
    gnn_hidden_units = c(64, 32, 16),
    gnn_activation = "relu",
    lstm_units = 64,
    hidden_size_1 = 64,
    hidden_size_2 = 32,
    hidden_size_3 = 16,
    age_embedding_dim = 8,
    region_embedding_dim = 8,
    dropout_rate = 0.1,
    use_age_embedding = TRUE,
    use_region_embedding = TRUE,
    use_age_numeric = FALSE,
    use_year_numeric = FALSE,
    age_numeric_dim = NULL,
    n_quantiles = 7L
) {
    base_model <- gnn_lstm_scalar_output_region_preserving(
        seq_length = seq_length,
        n_regions = n_regions,
        features_per_region = features_per_region,
        adjacency_matrix = adjacency_matrix,
        n_ages = n_ages,
        gnn_hidden_units = gnn_hidden_units,
        gnn_activation = gnn_activation,
        lstm_units = lstm_units,
        hidden_size_1 = hidden_size_1,
        hidden_size_2 = hidden_size_2,
        hidden_size_3 = hidden_size_3,
        age_embedding_dim = age_embedding_dim,
        region_embedding_dim = region_embedding_dim,
        dropout_rate = dropout_rate,
        use_age_embedding = use_age_embedding,
        use_region_embedding = use_region_embedding,
        use_age_numeric = use_age_numeric,
        use_year_numeric = use_year_numeric,
        age_numeric_dim = age_numeric_dim
    )

    hidden_output <- base_model$layers[[length(base_model$layers) - 1]]$output
    output <- hidden_output %>%
        layer_dense(units = as.integer(n_quantiles), name = "quantile_output")

    model <- keras_model(
        inputs = base_model$inputs,
        outputs = output,
        name = "gnn_lstm_quantile"
    )

    return(model)
}

# ------------------------------------------------------------------------------
# Function: make_sequences_age_aware
# Description: <What it does in 1-2 sentences>
# Inputs:
#   - df_key: <description>
#   - # data.frame with Week: <description>
#   - Region columns
# X_climate: <description>
#   - # matrix of climate features (nrow: <description>
# Output:
#   - <What is returned>
# ------------------------------------------------------------------------------
make_sequences_age_aware <- function(
    df_key, # data.frame with Week, Region columns
    X_climate, # matrix of climate features (nrow = nrow(df_key))
    y_mat, # matrix of targets, columns are ages
    seq_length,
    region_levels, # factor levels for regions
    age_cols # column names/indices for ages
) {
    n_ages <- length(age_cols)
    n_regions <- length(region_levels)

    # First, create per-region climate sequences (same as before)
    regions <- sort(unique(df_key$Region))

    # Storage
    X_list <- list()
    y_list <- list()
    region_list <- list()
    age_list <- list()
    week_list <- list()
    region_name_list <- list()
    age_name_list <- list()
    age_mid_list <- list()

    sample_idx <- 1

    for (r in regions) {
        idx <- which(df_key$Region == r)
        idx <- idx[order(df_key$Week[idx])]

        n_r <- length(idx)
        if (n_r < seq_length) {
            next
        }

        x_sub <- X_climate[idx, , drop = FALSE]
        y_sub <- y_mat[idx, , drop = FALSE]
        wk_sub <- df_key$Week[idx]

        n_seqs <- n_r - seq_length + 1

        # Get region index
        region_idx <- which(region_levels == r)

        for (i in seq_len(n_seqs)) {
            start <- i
            end <- i + seq_length - 1
            target_row <- end

            x_seq <- x_sub[start:end, , drop = FALSE]

            # For each age, create a separate sample
            for (a in seq_len(n_ages)) {
                X_list[[sample_idx]] <- x_seq
                y_list[[sample_idx]] <- y_sub[target_row, a]
                region_list[[sample_idx]] <- region_idx
                age_list[[sample_idx]] <- a
                week_list[[sample_idx]] <- wk_sub[target_row]
                region_name_list[[sample_idx]] <- r
                age_name_list[[sample_idx]] <- age_cols[a]
                # compute midpoint for age group (e.g. 65 -> 67.5, 90 -> 92.5)
                age_label <- age_cols[a]
                age_num <- suppressWarnings(as.numeric(as.character(age_label)))
                if (is.na(age_num)) {
                    m <- regmatches(
                        as.character(age_label),
                        regexpr("\\d+", as.character(age_label))
                    )
                    age_num <- ifelse(length(m) > 0, as.numeric(m), NA_real_)
                }
                age_mid_val <- ifelse(
                    is.finite(age_num),
                    age_num + 2.5,
                    NA_real_
                )
                age_mid_list[[sample_idx]] <- age_mid_val
                sample_idx <- sample_idx + 1
            }
        }
    }

    n_total <- sample_idx - 1
    n_features <- ncol(X_climate)

    # Convert to arrays
    X_arr <- array(0, dim = c(n_total, seq_length, n_features))
    for (i in seq_len(n_total)) {
        X_arr[i, , ] <- X_list[[i]]
    }

    y_vec <- unlist(y_list)
    region_vec <- unlist(region_list)
    age_vec <- unlist(age_list)
    week_vec <- as.Date(unlist(week_list), origin = "1970-01-01")
    region_name_vec <- unlist(region_name_list)
    age_name_vec <- unlist(age_name_list)
    age_mid_vec <- as.numeric(unlist(age_mid_list))

    list(
        X = X_arr,
        y = y_vec,
        region_idx = region_vec,
        age_idx = age_vec,
        target_week = week_vec,
        Region = region_name_vec,
        Age = age_name_vec,
        age_mid = age_mid_vec
    )
}

# ------------------------------------------------------------------------------
# Function: make_sequences_gnn_scalar
# Description: Creates temporal sequences for the scalar GNN-LSTM model, generating one sample per week, age, and region.
# Inputs:
#   - df_key: data.frame with Week and Region columns.
#   - X_mat: Matrix of climate features corresponding to df_key.
#   - y_mat: Matrix of scalar targets, where columns are ages.
#   - seq_length: Length of the temporal sequence.
#   - regions_list: List of region identifiers.
#   - age_cols: Column names corresponding to different ages.
# Output:
#   - A list containing X sequences, y targets, and metadata arrays.
# ------------------------------------------------------------------------------
make_sequences_gnn_scalar <- function(
    df_key, # data.frame with Week, Region columns
    X_mat, # matrix of climate features (nrow = nrow(df_key))
    y_mat, # matrix of targets, columns are ages
    seq_length,
    regions_list,
    age_cols # column names for ages
) {
    n_regions <- length(regions_list)
    n_ages <- length(age_cols)
    n_features_per_region <- ncol(X_mat)

    weeks_unique <- sort(unique(df_key$Week))
    n_weeks <- length(weeks_unique)

    # Reshape X_mat to 3D: (n_weeks, n_regions, n_features_per_region)
    X_3d <- array(0, dim = c(n_weeks, n_regions, n_features_per_region))

    for (w_idx in seq_along(weeks_unique)) {
        w <- weeks_unique[w_idx]
        for (r_idx in seq_along(regions_list)) {
            r <- regions_list[r_idx]
            idx <- which(df_key$Week == w & df_key$Region == r)
            if (length(idx) > 0) {
                X_3d[w_idx, r_idx, ] <- X_mat[idx[1], ]
            }
        }
    }

    # Also reshape y_mat to 3D: (n_weeks, n_regions, n_ages)
    y_3d <- array(0, dim = c(n_weeks, n_regions, n_ages))

    for (w_idx in seq_along(weeks_unique)) {
        w <- weeks_unique[w_idx]
        for (r_idx in seq_along(regions_list)) {
            r <- regions_list[r_idx]
            idx <- which(df_key$Week == w & df_key$Region == r)
            if (length(idx) > 0) {
                y_3d[w_idx, r_idx, ] <- y_mat[idx[1], ]
            }
        }
    }

    # Number of temporal sequences
    n_seqs <- n_weeks - seq_length + 1

    # Total samples: n_seqs * n_ages * n_regions
    n_total <- n_seqs * n_ages * n_regions

    # Allocate arrays
    X_seq <- array(
        0,
        dim = c(n_total, seq_length, n_regions * n_features_per_region)
    )
    y_vec <- numeric(n_total) # SCALAR target
    age_idx_vec <- integer(n_total)
    region_idx_vec <- integer(n_total)
    target_week_vec <- as.Date(rep(NA, n_total))
    age_name_vec <- character(n_total)
    region_name_vec <- character(n_total)

    sample_idx <- 1

    for (i in seq_len(n_seqs)) {
        seq_weeks_idx <- i:(i + seq_length - 1)
        target_week_idx <- i + seq_length - 1

        # Build spatial-temporal input: (seq_length, n_regions * n_features)
        x_seq_i <- X_3d[seq_weeks_idx, , , drop = FALSE]
        x_flat <- matrix(
            x_seq_i,
            nrow = seq_length,
            ncol = n_regions * n_features_per_region
        )

        # For each age and region, create a separate sample
        for (a in seq_len(n_ages)) {
            for (r in seq_len(n_regions)) {
                X_seq[sample_idx, , ] <- x_flat
                # Target: SINGLE residual for this age and region
                y_vec[sample_idx] <- y_3d[target_week_idx, r, a]
                age_idx_vec[sample_idx] <- a
                region_idx_vec[sample_idx] <- r
                target_week_vec[sample_idx] <- weeks_unique[target_week_idx]
                age_name_vec[sample_idx] <- age_cols[a]
                region_name_vec[sample_idx] <- regions_list[r]
                sample_idx <- sample_idx + 1
            }
        }
    }

    # Replace any NA targets with 0
    y_vec[!is.finite(y_vec)] <- 0

    list(
        X = X_seq,
        y = y_vec,
        age_idx = age_idx_vec,
        region_idx = region_idx_vec,
        target_week = target_week_vec,
        Age = age_name_vec,
        Region = region_name_vec,
        n_regions = n_regions,
        n_ages = n_ages,
        features_per_region = n_features_per_region,
        regions = regions_list,
        age_cols = age_cols
    )
}

# ------------------------------------------------------------------------------
# Function: cnn_lstm_age_aware
# Description: <What it does in 1-2 sentences>
# Inputs:
#   - seq_length: <description>
#   - n_climate_features: <description>
#   - n_regions: <description>
#   - n_ages: <description>
#   - cnn_filters: <description>
#   - 16: <description>
# Output:
#   - <What is returned>
# ------------------------------------------------------------------------------
cnn_lstm_age_aware <- function(
    seq_length,
    n_climate_features,
    n_regions,
    n_ages,
    cnn_filters = c(32, 16),
    lstm_units = 64,
    hidden_size_1 = 64,
    hidden_size_2 = 32,
    hidden_size_3 = 16,
    hidden_size_4 = 8,
    hidden_size_5 = 4,
    region_embedding_dim = 8,
    age_embedding_dim = 8,
    use_age_embedding = TRUE,
    use_region_embedding = TRUE,
    zero_embeddings = FALSE,
    dropout_rate = 0.2,
    # Optional numeric age feature: when TRUE the model accepts an
    # additional scalar input `age_numeric_input` (e.g. midpoint like 67.5)
    # which is projected and tiled similarly to embeddings.
    use_age_numeric = FALSE,
    age_numeric_dim = NULL,
    cell_type = c("lstm", "gru")
) {
    cell_type <- match.arg(cell_type)

    # Input 1: Climate sequence (seq_length, n_climate_features)
    climate_input <- layer_input(
        shape = c(seq_length, n_climate_features),
        name = "climate_input"
    )

    # Input 2: Region index (integer)
    region_input <- layer_input(
        shape = c(1),
        dtype = "int32",
        name = "region_input"
    )

    # Input 3: Age index (integer)
    age_input <- layer_input(
        shape = c(1),
        dtype = "int32",
        name = "age_input"
    )

    # Optional numeric age input (e.g. 67.5, 72.5). Only created when
    # use_age_numeric == TRUE.
    age_numeric_input <- NULL
    age_numeric_tiled <- NULL
    age_numeric_flat <- NULL
    if (isTRUE(use_age_numeric)) {
        age_numeric_input <- layer_input(
            shape = c(1),
            dtype = 'float32',
            name = "age_numeric_input"
        )
        # determine projection dim
        if (is.null(age_numeric_dim)) {
            age_numeric_dim_local <- as.integer(age_embedding_dim)
        } else {
            age_numeric_dim_local <- as.integer(age_numeric_dim)
        }
        # project scalar to vector and repeat across timesteps for pre-LSTM concat
        # Use a non-negative kernel constraint so the projection is monotone in the input age
        age_numeric_proj <- age_numeric_input %>%
            layer_dense(
                units = age_numeric_dim_local,
                activation = 'linear',
                name = "age_numeric_proj",
                kernel_constraint = constraint_nonneg()
            )
        age_numeric_tiled <- age_numeric_proj %>%
            layer_repeat_vector(as.integer(seq_length))
        # flat version for skip-connection after LSTM (also constrained to be monotone)
        age_numeric_flat <- age_numeric_input %>%
            layer_dense(
                units = age_numeric_dim_local,
                activation = 'linear',
                name = "age_numeric_flat",
                kernel_constraint = constraint_nonneg()
            ) %>%
            layer_flatten()
    }

    # Process climate through CNN
    # LayerNorm instead of BatchNorm: stable with any batch size and sample weights
    x <- climate_input
    for (f in cnn_filters) {
        x <- x %>%
            layer_conv_1d(
                filters = f,
                kernel_size = 3,
                padding = "same"
            ) %>%
            layer_layer_normalization() %>%
            layer_activation_leaky_relu(alpha = 0.1) %>%
            layer_spatial_dropout_1d(rate = dropout_rate)
    }

    # Region embedding
    region_embedding_enabled <- isTRUE(use_region_embedding) &&
        !isTRUE(zero_embeddings)
    if (region_embedding_enabled) {
        region_emb <- region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(max(1, region_embedding_dim)),
                name = "region_embedding"
            ) %>%
            layer_flatten()
    } else {
        region_emb <- region_input %>%
            layer_embedding(
                input_dim = as.integer(n_regions + 1),
                output_dim = as.integer(max(1, region_embedding_dim)),
                name = "region_embedding",
                embeddings_initializer = "zeros",
                trainable = FALSE
            ) %>%
            layer_flatten()
    }

    # Age embedding
    age_embedding_enabled <- isTRUE(use_age_embedding) &&
        !isTRUE(zero_embeddings)
    if (age_embedding_enabled) {
        age_emb <- age_input %>%
            layer_embedding(
                input_dim = as.integer(n_ages + 1),
                output_dim = as.integer(max(1, age_embedding_dim)),
                name = "age_embedding"
            ) %>%
            layer_flatten()
    } else {
        age_emb <- age_input %>%
            layer_embedding(
                input_dim = as.integer(n_ages + 1),
                output_dim = as.integer(max(1, age_embedding_dim)),
                name = "age_embedding",
                embeddings_initializer = "zeros",
                trainable = FALSE
            ) %>%
            layer_flatten()
    }

    # Tile embeddings across timesteps and concatenate with climate features
    region_emb_tiled <- region_emb %>%
        layer_repeat_vector(as.integer(seq_length))

    age_emb_tiled <- age_emb %>%
        layer_repeat_vector(as.integer(seq_length))

    # Concatenate: (seq_length, cnn_output + region_emb + age_emb)
    # Note: we intentionally do NOT include the numeric-age tile here. Keeping the
    # numeric age out of the pre-LSTM concatenation makes the age effect act more
    # additively (via the skip connection below) which helps preserve global
    # monotonicity in age when combined with the (potentially interacting) LSTM/climate branch.
    concat_list <- list(x)
    if (!is.null(region_emb_tiled)) {
        concat_list <- c(concat_list, list(region_emb_tiled))
    }
    if (!is.null(age_emb_tiled)) {
        concat_list <- c(concat_list, list(age_emb_tiled))
    }
    x <- layer_concatenate(concat_list, axis = -1)

    # LSTM/GRU layer
    if (cell_type == "gru") {
        x <- x %>%
            layer_gru(units = lstm_units, return_sequences = FALSE)
    } else {
        x <- x %>%
            layer_lstm(units = lstm_units, return_sequences = FALSE)
    }

    x <- x %>% layer_layer_normalization()

    # Also concatenate embeddings (and optional numeric age) directly to LSTM output for skip connection
    # Keep region embedding (and numeric age) in the skip connection, but move the
    # explicit `age_emb` concatenation to the final output stage so the head is
    # explicitly conditioned on age identity.
    skip_list <- list(x, region_emb)
    if (!is.null(age_numeric_flat)) {
        skip_list <- c(skip_list, list(age_numeric_flat))
    }
    x <- layer_concatenate(skip_list, axis = -1)

    # Dense head: Dense → LayerNorm → LeakyReLU → Dropout
    # (following MortFCNet best practice; avoids dying ReLU with small signals)
    x <- x %>%
        layer_dense(units = as.integer(hidden_size_1)) %>%
        layer_layer_normalization() %>%
        layer_activation_leaky_relu(alpha = 0.1) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(units = as.integer(hidden_size_2)) %>%
        layer_layer_normalization() %>%
        layer_activation_leaky_relu(alpha = 0.1) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(units = as.integer(hidden_size_3)) %>%
        layer_layer_normalization() %>%
        layer_activation_leaky_relu(alpha = 0.1) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(units = as.integer(hidden_size_4)) %>%
        layer_layer_normalization() %>%
        layer_activation_leaky_relu(alpha = 0.1) %>%
        layer_dropout(rate = dropout_rate) %>%
        layer_dense(units = as.integer(hidden_size_5)) %>%
        layer_layer_normalization() %>%
        layer_activation_leaky_relu(alpha = 0.1) %>%
        layer_dropout(rate = dropout_rate)

    # Age-conditioning: concatenate the age embedding right before the final output
    # so the final projection must use the explicit age identity.
    x <- layer_concatenate(list(x, age_emb), axis = -1)

    # Output: single residual prediction
    output <- x %>%
        layer_dense(units = 1, name = "residual_output")

    # Build inputs list conditionally including the numeric age input
    model_inputs <- list(climate_input, region_input, age_input)
    if (!is.null(age_numeric_input)) {
        model_inputs <- c(model_inputs, list(age_numeric_input))
    }

    model <- keras_model(
        inputs = model_inputs,
        outputs = output
    )

    return(model)
}

# ------------------------------------------------------------------------------
# Function: cnn_lstm_age_aware_quantile
# Description: Quantile-regression variant of cnn_lstm_age_aware.
cnn_lstm_age_aware_quantile <- function(
    seq_length,
    n_climate_features,
    n_regions,
    n_ages,
    cnn_filters = c(32, 16),
    lstm_units = 64,
    hidden_size_1 = 64,
    hidden_size_2 = 32,
    hidden_size_3 = 16,
    hidden_size_4 = 8,
    hidden_size_5 = 4,
    region_embedding_dim = 8,
    age_embedding_dim = 8,
    use_age_embedding = TRUE,
    use_region_embedding = TRUE,
    zero_embeddings = FALSE,
    dropout_rate = 0.2,
    use_age_numeric = FALSE,
    age_numeric_dim = NULL,
    cell_type = c("lstm", "gru"),
    n_quantiles = 7L
) {
    base_model <- cnn_lstm_age_aware(
        seq_length = seq_length,
        n_climate_features = n_climate_features,
        n_regions = n_regions,
        n_ages = n_ages,
        cnn_filters = cnn_filters,
        lstm_units = lstm_units,
        hidden_size_1 = hidden_size_1,
        hidden_size_2 = hidden_size_2,
        hidden_size_3 = hidden_size_3,
        hidden_size_4 = hidden_size_4,
        hidden_size_5 = hidden_size_5,
        region_embedding_dim = region_embedding_dim,
        age_embedding_dim = age_embedding_dim,
        use_age_embedding = use_age_embedding,
        use_region_embedding = use_region_embedding,
        zero_embeddings = zero_embeddings,
        dropout_rate = dropout_rate,
        use_age_numeric = use_age_numeric,
        age_numeric_dim = age_numeric_dim,
        cell_type = cell_type
    )

    hidden_output <- base_model$layers[[length(base_model$layers) - 1]]$output
    output <- hidden_output %>%
        layer_dense(units = as.integer(n_quantiles), name = "quantile_output")

    model <- keras_model(
        inputs = base_model$inputs,
        outputs = output,
        name = "cnn_lstm_quantile"
    )

    return(model)
}
