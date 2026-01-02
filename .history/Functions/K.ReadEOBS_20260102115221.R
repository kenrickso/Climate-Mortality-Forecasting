# Utility to list and read E-OBS region files like: FR_NUTS2_xx_daily_2013-2024.rds

# List available variables (the "xx") for files matching COUNTRY_NUTS_XX_daily_2013-2024.rds
list_eobs_files <- function(path = "Data/E-OBS", country = "FR", nuts = "NUTS2") {
  pattern <- paste0("^", country, "_", nuts, "_(.*)_daily_2013-2024\\.rds$")
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(character(0))
  vars <- sub(pattern, "\\1", basename(files))
  setNames(files, vars)
}

# Read one or more E-OBS region files.
# - vars: character vector of variable codes (e.g. "hu","tg"); if NULL reads all available.
# - combine: if TRUE returns a single data.frame with an added "variable" column; else returns a named list of data.frames.
read_eobs_files <- function(path = "../../Data/E-OBS", country = "FR", nuts = "NUTS2",
                            vars = NULL, combine = FALSE) {
  files_named <- list_eobs_files(path = path, country = country, nuts = nuts)
  if (length(files_named) == 0) {
    warning("No matching E-OBS files found in ", path)
    return(if (combine) data.frame() else list())
  }
  
  if (!is.null(vars)) {
    # keep only requested vars that exist
    missing_vars <- setdiff(vars, names(files_named))
    if (length(missing_vars)) {
      warning("Requested vars not found and will be ignored: ", paste(missing_vars, collapse = ", "))
    }
    files_named <- files_named[intersect(vars, names(files_named))]
  }
  
  # read files and normalize each data.frame: single climate column -> "value", ensure Date is Date, add variable column
  read_list <- lapply(names(files_named), function(v) {
    dat_raw <- readRDS(files_named[[v]])
    # coerce to data.frame for consistent handling
    dat <- as.data.frame(dat_raw)
    
    # normalize Date and Region column names (case-insensitive common variants)
    nm <- names(dat)
    date_idx <- which(tolower(nm) %in% c("date", "day"))
    region_idx <- which(tolower(nm) %in% c("region", "nuts_id", "nuts", "region_id"))
    if (length(date_idx) == 1) names(dat)[date_idx] <- "Date"
    if (length(region_idx) == 1) names(dat)[region_idx] <- "Region"
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
        numeric_cols <- climate_cols[sapply(dat[climate_cols], is.numeric)]
        if (length(numeric_cols) > 0) {
          na_counts <- sapply(dat[numeric_cols], function(col) sum(is.na(col)))
          chosen <- numeric_cols[which.min(na_counts)]
        } else {
          # fallback: choose the first climate col
          chosen <- climate_cols[1]
        }
      }
      dat$value <- dat[[chosen]]
      # drop the original climate columns (except the new 'value')
      to_keep <- c("Date", "Region", "variable", "value", setdiff(names(dat), c(climate_cols)))
      dat <- dat[to_keep]
    }
    
    # ensure mandatory columns exist
    if (!"Region" %in% names(dat)) dat$Region <- NA
    if (!"Date" %in% names(dat)) dat$Date <- NA
    dat
  })
  names(read_list) <- names(files_named)
  
  if (combine) {
    # keep only data.frames (drop unexpected objects)
    dfs <- Filter(function(x) is.data.frame(x), read_list)
    if (length(dfs) == 0) return(data.frame())
    
    # ensure each df has the core columns and align order
    dfs_core <- lapply(dfs, function(df) {
      # ensure Date, Region, value, variable exist
      if (!"Date" %in% names(df)) df$Date <- NA
      if (!"Region" %in% names(df)) df$Region <- NA
      if (!"value" %in% names(df)) {
        # try to detect a single remaining climate-like column
        other_cols <- setdiff(names(df), c("Date","Region","variable"))
        if (length(other_cols) >= 1) df$value <- df[[other_cols[1]]] else df$value <- NA
      }
      if (!"variable" %in% names(df)) df$variable <- NA
      # keep only the core long-format columns
      df[, intersect(c("Date","Region","variable","value"), names(df))]
    })
    
    # bind into a single long table
    long <- do.call(rbind, dfs_core)
    rownames(long) <- NULL
    # ensure Date is Date
    if ("Date" %in% names(long)) long$Date <- as.Date(long$Date)
    
    # pivot to wide: variables become column names
    if (requireNamespace("tidyr", quietly = TRUE)) {
      wide <- tidyr::pivot_wider(long, names_from = "variable", values_from = "value")
    } else {
      # base R fallback using reshape
      wide <- reshape(long[, c("Date","Region","variable","value")],
                      idvar = c("Date","Region"), timevar = "variable", direction = "wide")
      # reshape creates column names like value.<var>; clean them
      names(wide) <- sub("^value\\.", "", names(wide))
    }
    # return wide data.frame (Date, Region, <vars...>)
    return(wide)
  }
  read_list
}