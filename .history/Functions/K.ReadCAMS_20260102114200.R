# R6 class wrapping file reading & aggregation utilities for datasets

# null-coalesce helper
`%||%` <- function(x, y) if (is.null(x)) y else x

if (!requireNamespace("R6", quietly = TRUE)) stop("Please install the 'R6' package to use DatasetAggregator")

DatasetAggregator <- R6::R6Class(
  "DatasetAggregator",
  public = list(
    dir = NULL,
    pattern = NULL,
    date_col = NULL,
    region_col = NULL,
    agg_fun = NULL,
    na.rm = NULL,
    read_fun = NULL,
    files = NULL,
    aggregated = NULL,

    initialize = function(
      dir = "C:\1 My Code\FullDataset\DataData/CAMS",
      pattern = "\\.txt$",
      date_col = "Date",
      region_col = "Region",
      agg_fun = mean, # default provided here (local var)
      na.rm = TRUE,
      read_fun = NULL
    ) {
      # assign basic scalar defaults to public fields (no locked-binding issue)
      self$dir <- dir
      self$pattern <- pattern
      self$date_col <- date_col
      self$region_col <- region_col

      # set function and logical defaults inside initialize
      # accept either a function or a string name for agg_fun
      if (is.null(agg_fun)) {
        self$agg_fun <- mean
      } else if (is.character(agg_fun) || is.symbol(agg_fun)) {
        self$agg_fun <- match.fun(agg_fun)
      } else {
        self$agg_fun <- agg_fun
      }

      self$na.rm <- na.rm

      if (is.null(read_fun)) {
        self$read_fun <- if (requireNamespace("data.table", quietly = TRUE)) {
          data.table::fread
        } else {
          utils::read.csv
        }
      } else {
        self$read_fun <- read_fun
      }
      self$refresh_files()
    },

    refresh_files = function() {
      self$files <- list.files(
        self$dir,
        pattern = self$pattern,
        full.names = TRUE
      )
      invisible(self$files)
    },

    head_all_files = function(n = 6, read_fun = NULL, ...) {
      read_fun <- read_fun %||% self$read_fun
      files <- self$files
      if (length(files) == 0) {
        stop("No files to show. Call refresh_files() or check dir.")
      }
      out <- list()
      for (f in files) {
        cat("\n---", basename(f), "---\n")
        df_head <- tryCatch(
          read_fun(f, nrows = n, stringsAsFactors = FALSE, ...),
          error = function(e) {
            message("Could not read ", basename(f), ": ", conditionMessage(e))
            return(NULL)
          }
        )
        print(df_head)
        out[[basename(f)]] <- df_head
      }
      invisible(out)
    },

    head_all_files_fread = function(n = 6, ...) {
      if (!requireNamespace("data.table", quietly = TRUE)) {
        return(self$head_all_files(n = n, ...))
      }
      files <- self$files
      if (length(files) == 0) {
        stop("No files to show. Call refresh_files() or check dir.")
      }
      out <- list()
      for (f in files) {
        cat("\n---", basename(f), "---\n")
        df_head <- tryCatch(
          data.table::fread(f, nrows = n, showProgress = FALSE, ...),
          error = function(e) {
            message("Could not fread ", basename(f), ": ", conditionMessage(e))
            return(NULL)
          }
        )
        print(df_head)
        out[[basename(f)]] <- df_head
      }
      invisible(out)
    },

    read_and_agg_file = function(
      file,
      date_col = NULL,
      region_col = NULL,
      agg_fun = NULL,
      na.rm = NULL,
      read_fun = NULL,
      ...
    ) {
      date_col <- date_col %||% self$date_col
      region_col <- region_col %||% self$region_col
      agg_fun <- agg_fun %||% self$agg_fun
      na.rm <- na.rm %||% self$na.rm
      read_fun <- read_fun %||% self$read_fun

      df <- tryCatch(
        read_fun(file, stringsAsFactors = FALSE, ...),
        error = function(e) {
          stop("Read failed for ", basename(file), ": ", conditionMessage(e))
        }
      )
      df <- as.data.frame(df, stringsAsFactors = FALSE)

      if (!date_col %in% names(df)) {
        stop("Date column '", date_col, "' not found in ", basename(file))
      }
      parsed <- as.POSIXct(
        df[[date_col]],
        tz = "UTC",
        tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d")
      )
      if (all(is.na(parsed))) {
        parsed <- as.Date(df[[date_col]])
      }
      df[[date_col]] <- as.Date(parsed)

      if (!(region_col %in% names(df))) {
        df[[region_col]] <- NA_character_
      } else {
        df[[region_col]] <- as.character(df[[region_col]])
      }

      num_cols <- setdiff(
        names(df)[sapply(df, is.numeric)],
        c(date_col, region_col)
      )
      if (length(num_cols) == 0) {
        stop("No numeric columns to aggregate in ", basename(file))
      }

      # use the actual column names provided by the user for grouping (no hard-coded "Region"/"Date")
      agg_by <- setNames(
        list(df[[region_col]], df[[date_col]]),
        c(region_col, date_col)
      )
      agg_df <- aggregate(df[num_cols], by = agg_by, FUN = function(x) {
        agg_fun(x, na.rm = na.rm)
      })

      prefix <- tools::file_path_sans_ext(basename(file))
      new_names <- c(region_col, date_col, paste0(prefix, "_", num_cols))

      names(agg_df) <- new_names
      agg_df[[date_col]] <- as.Date(agg_df[[date_col]])
      agg_df[[region_col]] <- as.character(agg_df[[region_col]])
      as.data.frame(agg_df, stringsAsFactors = FALSE)
    },

    aggregate_files = function(
      files = NULL,
      date_col = NULL,
      region_col = NULL,
      agg_fun = NULL,
      na.rm = NULL,
      read_fun = NULL,
      ...
    ) {
      date_col <- date_col %||% self$date_col
      region_col <- region_col %||% self$region_col
      agg_fun <- agg_fun %||% self$agg_fun
      na.rm <- na.rm %||% self$na.rm
      read_fun <- read_fun %||% self$read_fun
      files <- files %||% self$files
      if (length(files) == 0) {
        stop("No files matching pattern in dir: ", self$dir)
      }

      lst <- lapply(files, function(f) {
        tryCatch(
          self$read_and_agg_file(
            f,
            date_col = date_col,
            region_col = region_col,
            agg_fun = agg_fun,
            na.rm = na.rm,
            read_fun = read_fun,
            ...
          ),
          error = function(e) {
            message("Skipping ", basename(f), ": ", conditionMessage(e))
            NULL
          }
        )
      })
      lst <- Filter(Negate(is.null), lst)
      if (length(lst) == 0) {
        stop("No files could be read and aggregated.")
      }

      lst <- lapply(lst, function(x) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        if (!(region_col %in% names(x))) {
          x[[region_col]] <- NA_character_
        }
        x[[region_col]] <- as.character(x[[region_col]])
        x[[date_col]] <- as.Date(x[[date_col]])
        x
      })

      by_cols <- c(region_col, date_col)
      aggregated_df <- Reduce(
        function(a, b) merge(a, b, by = by_cols, all = TRUE),
        lst
      ) # generic name

      aggregated_df[[region_col]] <- as.character(aggregated_df[[region_col]])
      aggregated_df[[date_col]] <- as.Date(aggregated_df[[date_col]])
      aggregated_df <- aggregated_df[
        order(aggregated_df[[region_col]], aggregated_df[[date_col]]),
        ,
        drop = FALSE
      ]
      rownames(aggregated_df) <- NULL
      self$aggregated <- as.data.frame(aggregated_df, stringsAsFactors = FALSE)
      invisible(self$aggregated)
    },

    get_aggregated = function() {
      self$aggregated
    },

    select_avg_cols = function(pattern = "avg", ignore.case = TRUE) {
      df <- self$get_aggregated()
      if (is.null(df)) {
        return(NULL)
      }
      keep <- c(
        self$region_col,
        self$date_col,
        grep(pattern, names(df), value = TRUE, ignore.case = ignore.case)
      )
      df[, intersect(keep, names(df)), drop = FALSE]
    }
  )
)