#import pacakges
library(data.table)
library(stringr)
library(rlang)
library(httr2)

# load .env if CENSUS_API_KEY is not already set
if (Sys.getenv("CENSUS_API_KEY") == "") {
  script_dir <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)), error = function(e) ".")
  candidates <- c(file.path(script_dir, ".env"), ".env")
  env_file   <- Filter(file.exists, candidates)[1]
  if (!is.na(env_file)) {
    lines <- readLines(env_file, warn = FALSE)
    for (line in lines) {
      if (grepl("^CENSUS_API_KEY=", line)) {
        key <- gsub("[^A-Za-z0-9]", "", sub("^CENSUS_API_KEY=", "", line))
        Sys.setenv(CENSUS_API_KEY = key)
      }
    }
  }
}
if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("CENSUS_API_KEY is not set — add it to .env or set it in your environment")
}

# functions
#get_vars function
getvars <- function(formula) {
  constr <- gsub("\\(|\\)", "", formula)
  constr <- gsub("\\* 100", "", constr)
  vars   <- unlist(strsplit(constr, "[\\+]|[\\-]|[\\/]|[\\*]"))
  vars   <- gsub("[[:space:]]", "", vars)
  vars   <- toupper(vars)
  vars   <- vars[!duplicated(vars)]
  return(vars)
}
#Parse funciton 
parse_formula <- function(formula_str) {
  # 1. detect and strip "* 100"
  multiply_100 <- grepl("\\* 100", formula_str)
  constr <- gsub("\\* 100", "", formula_str)

  # 2. strip parentheses, trim whitespace
  constr <- gsub("\\(|\\)", "", constr)
  constr <- trimws(constr)

  # 3. split on "/" into numerator and denominator
  has_division <- grepl("/", constr)
  if (has_division) {
    parts <- strsplit(constr, "/")[[1]]
    if (length(parts) != 2) {
      stop(paste0("Formula must have exactly one /: ", formula_str))
    }
    num_str <- trimws(parts[1])
    den_str <- trimws(parts[2])
  } else {
    num_str <- constr
    den_str <- NULL
  }

  # 4. helper: split one side into variables + operators
  parse_side <- function(s) {
    # capture the + and - operators in order
    op_matches <- gregexpr("[+\\-]", s)[[1]]
    if (op_matches[1] == -1) {
      operators <- character(0)
    } else {
      operators <- substring(s, op_matches, op_matches)
    }
    # split into variable names
    vars <- strsplit(s, "[+\\-]")[[1]]
    vars <- trimws(vars)
    vars <- toupper(vars)
    vars <- vars[vars != ""]
    list(vars = vars, ops = operators)
  }

  # 5. parse each side
  num_parsed <- parse_side(num_str)
  den_parsed <- if (!is.null(den_str)) parse_side(den_str) else list(vars = character(0), ops = character(0))

  # 6. return the structured result
  list(
    num_vars     = num_parsed$vars,
    num_ops      = num_parsed$ops,
    den_vars     = den_parsed$vars,
    den_ops      = den_parsed$ops,
    multiply_100 = multiply_100,
    has_division = has_division
  )
}
#Computing MOE engine. 
compute_var_vector <- function(df, est_cols, one_zero = TRUE, z = 1.645) {
  if (length(est_cols) == 0) return(rep(0, nrow(df)))

  # find the matching MOE columns (swap trailing E for M)
  moe_cols <- str_replace(est_cols, "E$", "M")

  # extract estimates and MOEs as matrices
  E_mat <- as.matrix(df[, ..est_cols])
  V_mat <- as.matrix(df[, ..moe_cols])

  # replace Census "not available" special value with NA
  V_mat[V_mat %in% c(-555555555)] <- NA
  E_mat[E_mat %in% c(-555555555)] <- NA

  # convert MOE to variance
  V_mat <- (V_mat / z)^2

  # simple case: no one-zero rule, or only one variable
  if (!one_zero || length(est_cols) == 1) {
    return(rowSums(V_mat, na.rm = TRUE))
  }

  # one-zero: sum variances where estimate is NOT zero
  nonzero_sum <- rowSums(V_mat * (E_mat != 0), na.rm = TRUE)

  # one-zero: max variance among the zero-estimate components
  zero_vars_mat <- V_mat * (E_mat == 0)
  zero_vars_mat[is.na(zero_vars_mat)] <- 0
  zero_max <- exec(pmax, !!!as.data.frame(zero_vars_mat))

  nonzero_sum + zero_max
}

moe_proprat <- function(df, parsed, method, one_zero = TRUE, z = 1.645) {
  # column names: add E to each variable to get the estimate columns
  num_cols <- paste0(parsed$num_vars, "E")
  den_cols <- paste0(parsed$den_vars, "E")

  # numerator estimate = sum of the numerator columns
  numerator   <- rowSums(as.matrix(df[, ..num_cols]), na.rm = TRUE)
  # denominator estimate = sum of the denominator columns
  denominator <- rowSums(as.matrix(df[, ..den_cols]), na.rm = TRUE)

  # the derived value
  p <- numerator / denominator

  # variances of numerator and denominator (one-zero aware)
  var_num <- compute_var_vector(df, num_cols, one_zero, z)
  var_den <- compute_var_vector(df, den_cols, one_zero, z)

  # proportion vs ratio: only the sign under the root differs
  if (method %in% c("proportion", "prop")) {
    var_prop  <- (var_num - p^2 * var_den) / (denominator^2)
    var_ratio <- (var_num + p^2 * var_den) / (denominator^2)
    # use proportion formula, fall back to ratio if it goes negative
    var_final <- ifelse(!is.na(var_prop) & var_prop >= 0, var_prop, var_ratio)
  } else {
    var_final <- (var_num + p^2 * var_den) / (denominator^2)
  }

  est <- p
  moe <- sqrt(pmax(var_final, 0)) * z

  # scale to percent if the formula had "* 100"
  if (parsed$multiply_100) {
    est <- est * 100
    moe <- moe * 100
  }

  # zero denominator -> undefined
  est <- ifelse(denominator == 0, NA_real_, est)
  moe <- ifelse(denominator == 0, NA_real_, moe)

  data.table(est = est, moe = moe)
}
moe_agg <- function(df, parsed, one_zero = TRUE, z = 1.645) {
  # estimate columns
  est_cols <- paste0(parsed$num_vars, "E")

  # estimate = sum of the columns
  est <- rowSums(as.matrix(df[, ..est_cols]), na.rm = TRUE)

  # variance from the engine, then convert to MOE
  variance <- compute_var_vector(df, est_cols, one_zero, z)
  moe <- sqrt(variance) * z

  data.table(est = est, moe = moe)
}

##Traffic cop fucntion, main loop reads it row by row and would need a whole lot of branching logic, so this compute_one_variable function decides which method to compute instead of all of this branching logic being in the main loop. 
#main loop only has a simple for(eachrow){compute_one_variable(df, parsed, method)}
compute_one_variable <- function(df, parsed, method, one_zero = TRUE, z = 1.645) {
  method <- tolower(method)

  if (method == "variable") {
    # just grab the one column's estimate and MOE directly
    est_col <- paste0(parsed$num_vars[1], "E")
    moe_col <- paste0(parsed$num_vars[1], "M")
    return(data.table(est = df[[est_col]], moe = df[[moe_col]]))
  }

  if (method %in% c("aggregation", "agg")) {
    return(moe_agg(df, parsed, one_zero, z))
  }

  if (method %in% c("proportion", "prop", "ratio")) {
    return(moe_proprat(df, parsed, method, one_zero, z))
  }

  stop(paste0("Unknown method: ", method))
}

# ── geo_lookup ────────────────────────────────────────────────────────────────
# Maps level names to Census API geography parameters and GEO_ID prefix.
# needs_county: TRUE means the "in" clause must be "state:XX+county:*"
geo_lookup <- data.table(
  level = c(
    "county", "county.subdivision", "tract", "block.group",
    "place", "school.district.unified", "school.district.elementary",
    "school.district.secondary"
  ),
  for_param = c(
    "county:*", "county subdivision:*", "tract:*", "block group:*",
    "place:*", "unified school district:*",
    "elementary school district:*", "secondary school district:*"
  ),
  needs_county = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  geo_id_prefix = c(
    "0500000US", "0600000US", "1400000US", "1500000US",
    "1600000US", "9700000US", "9500000US", "9600000US"
  )
)

# ── fetch_acs ─────────────────────────────────────────────────────────────────
# Fetches ACS 5-year estimates from the Census API.
# variables: bare variable names without E/M suffix (e.g. "B17001_001")
# year:      ACS release year (e.g. 2021)
# for_param: Census "for" clause  (e.g. "county:*")
# in_param:  Census "in"  clause  (e.g. "state:55" or "state:55+county:*")
# Returns a data.table with GEO_ID, NAME, and {var}E / {var}M columns.
fetch_acs <- function(variables, year, for_param, in_param,
                      api_key = Sys.getenv("CENSUS_API_KEY")) {
  base_url <- sprintf("https://api.census.gov/data/%d/acs/acs5", year)
  sentinels <- c(-555555555, -666666666, -999999999, -888888888)

  # split into chunks of 24 vars (= 48 data cols) to stay under the 50-col limit
  chunks <- split(variables, ceiling(seq_along(variables) / 24))

  fetch_chunk <- function(vars) {
    cols <- c("GEO_ID", "NAME",
              as.vector(rbind(paste0(vars, "E"), paste0(vars, "M"))))

    resp <- request(base_url) |>
      req_url_query(
        get   = paste(cols, collapse = ","),
        `for` = for_param,
        `in`  = in_param,
        key   = api_key
      ) |>
      req_perform()

    ct <- resp_content_type(resp)
    if (!grepl("json", ct)) {
      body <- resp_body_string(resp)
      # extract text from HTML error page if present
      msg <- regmatches(body, regexpr("(?<=<p>)[^<]+", body, perl = TRUE))
      stop("Census API error: ", if (length(msg)) trimws(msg) else body)
    }

    json    <- resp_body_json(resp, simplifyVector = FALSE)
    headers <- unlist(json[[1]])
    mat     <- do.call(rbind, lapply(json[-1], unlist))
    dt      <- as.data.table(mat)
    setnames(dt, headers)
    dt
  }

  batches <- lapply(chunks, fetch_chunk)

  # join all batches on GEO_ID + NAME
  out <- Reduce(
    function(a, b) merge(a, b, by = c("GEO_ID", "NAME"), all = TRUE),
    batches
  )

  # convert data columns to numeric and null out Census sentinel values
  data_cols <- setdiff(names(out), c("GEO_ID", "NAME"))
  out[, (data_cols) := lapply(.SD, as.numeric), .SDcols = data_cols]
  for (col in data_cols) {
    set(out, which(out[[col]] %in% sentinels), col, NA_real_)
  }

  out
}

# ── extract_fips ──────────────────────────────────────────────────────────────
# Parses the Census GEO_ID column into FIPS component columns.
# GEO_ID format: "{prefix}US{digits}" e.g. "0500000US55025"
# Adds st, and level-specific columns (cnty, tract, bg, cousub, place, sdid).
extract_fips <- function(df) {
  dt <- as.data.table(df)

  prefix  <- substr(dt$GEO_ID[1], 1, 9)
  geo_row <- geo_lookup[geo_id_prefix == prefix]
  if (nrow(geo_row) == 0) stop(paste("Unknown GEO_ID prefix:", prefix))

  lvl    <- geo_row$level
  suffix <- sub(".*US", "", dt$GEO_ID)

  dt[, st := substr(suffix, 1, 2)]

  if (lvl == "county") {
    dt[, cnty := substr(suffix, 3, 5)]

  } else if (lvl == "county.subdivision") {
    dt[, cnty   := substr(suffix, 3, 5)]
    dt[, cousub := substr(suffix, 6, 10)]

  } else if (lvl == "tract") {
    dt[, cnty  := substr(suffix, 3, 5)]
    dt[, tract := substr(suffix, 6, 11)]

  } else if (lvl == "block.group") {
    dt[, cnty  := substr(suffix, 3, 5)]
    dt[, tract := substr(suffix, 6, 11)]
    dt[, bg    := substr(suffix, 12, 12)]

  } else if (lvl == "place") {
    dt[, place := substr(suffix, 3, 7)]

  } else {
    dt[, sdid := substr(suffix, 3, 7)]
  }

  dt
}

# ── acsdata ───────────────────────────────────────────────────────────────────
# Fetches all ACS variables needed by a set of formula specs for one geography
# level, with saveRDS caching keyed on (level, year, state).
# formulas:  character vector of formula strings
# level:     one of geo_lookup$level
# year:      ACS release year
# state:     2-digit state FIPS as character (e.g. "55")
# cache_dir: directory for .rds cache files (default: tempdir())
acsdata <- function(formulas, level, year, state,
                    api_key   = Sys.getenv("CENSUS_API_KEY"),
                    cache_dir = tempdir()) {

  all_vars <- unique(unlist(lapply(formulas, getvars)))

  # fetch one level (the original single-level logic, pulled into a helper)
  fetch_one_level <- function(lvl) {
    cache_file <- file.path(
      cache_dir,
      sprintf("acs_%s_%s_%s.rds", lvl, year, state)
    )

    if (file.exists(cache_file)) {
      message("Loading from cache: ", cache_file)
      return(readRDS(cache_file))
    }

    geo_row <- geo_lookup[level == lvl]
    if (nrow(geo_row) == 0) stop(paste("Unknown level:", lvl))

    in_param <- if (geo_row$needs_county) {
      sprintf("state:%s+county:*", state)
    } else {
      sprintf("state:%s", state)
    }

    dt <- fetch_acs(all_vars, year, geo_row$for_param, in_param, api_key)
    dt <- extract_fips(dt)

    saveRDS(dt, cache_file)
    message("Cached to: ", cache_file)
    dt
  }

  # loop over every requested level, return a named list
  result <- setNames(lapply(level, fetch_one_level), level)
  result
}
# ── sumacs ────────────────────────────────────────────────────────────────────
# Main computation loop. For each row in spec, parses the formula and calls
# compute_one_variable, then assembles results into a wide data.table.
# spec:     data.table with columns: varname, formula, method
# data:     data.table from acsdata()
# one_zero: apply Census one-zero MOE rule (default TRUE)
# z:        z-score for MOE (1.645 = 90%, 1.96 = 95%)
sumacs <- function(spec, data, one_zero = TRUE, z = 1.645) {

  # compute all derived variables for a single level's table
  compute_level <- function(level_data) {
    geo_cols <- intersect(
      c("GEO_ID", "NAME", "st", "cnty", "cousub", "tract", "bg", "place", "sdid"),
      names(level_data)
    )

    results <- vector("list", nrow(spec))
    for (i in seq_len(nrow(spec))) {
      varname <- spec$varname[i]
      parsed  <- parse_formula(spec$formula[i])
      out     <- compute_one_variable(level_data, parsed, spec$method[i], one_zero, z)
      setnames(out, c("est", "moe"),
               c(paste0(varname, "_est"), paste0(varname, "_moe")))
      results[[i]] <- out
    }

    cbind(level_data[, ..geo_cols], do.call(cbind, results))
  }

  # if data is a single data.table, wrap it so the loop works uniformly
  if (is.data.table(data)) {
    data <- list(data)
  }

  # compute each level, then stack them with fill=TRUE
  # (fill handles levels having different FIPS columns, e.g. tract has 'tract', county doesn't)
  per_level <- lapply(data, compute_level)
  rbindlist(per_level, fill = TRUE, idcol = "level")
}