#import packages
library(data.table)
library(stringr)
library(rlang)
library(httr2)

# load .Renviron if CENSUS_API_KEY is not already set
# recommended: store key in ~/.Renviron as CENSUS_API_KEY=your_key_here
if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("CENSUS_API_KEY is not set — add it to ~/.Renviron or set it with Sys.setenv()")
}

# ── getvars ───────────────────────────────────────────────────────────────────
# Extracts all ACS variable names from a formula string.
getvars <- function(formula) {
  constr <- gsub("\\(|\\)", "", formula)
  constr <- gsub("\\* 100", "", constr)
  vars   <- unlist(strsplit(constr, "[\\+]|[\\-]|[\\/]|[\\*]"))
  vars   <- gsub("[[:space:]]", "", vars)
  vars   <- toupper(vars)
  vars   <- vars[!duplicated(vars)]
  return(vars)
}

# ── parse_formula ─────────────────────────────────────────────────────────────
# Breaks a formula string into numerator/denominator variable lists and flags * 100.
parse_formula <- function(formula_str) {
  multiply_100 <- grepl("\\* 100", formula_str)
  constr <- gsub("\\* 100", "", formula_str)
  constr <- gsub("\\(|\\)", "", constr)
  constr <- trimws(constr)

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

  parse_side <- function(s) {
    op_matches <- gregexpr("[+\\-]", s)[[1]]
    if (op_matches[1] == -1) {
      operators <- character(0)
    } else {
      operators <- substring(s, op_matches, op_matches)
    }
    vars <- strsplit(s, "[+\\-]")[[1]]
    vars <- trimws(vars)
    vars <- toupper(vars)
    vars <- vars[vars != ""]
    list(vars = vars, ops = operators)
  }

  num_parsed <- parse_side(num_str)
  den_parsed <- if (!is.null(den_str)) parse_side(den_str) else list(vars = character(0), ops = character(0))

  list(
    num_vars     = num_parsed$vars,
    num_ops      = num_parsed$ops,
    den_vars     = den_parsed$vars,
    den_ops      = den_parsed$ops,
    multiply_100 = multiply_100,
    has_division = has_division
  )
}

# ── compute_var_vector ────────────────────────────────────────────────────────
# Combines a set of MOEs into one variance number, with the one-zero rule.
compute_var_vector <- function(df, est_cols, one_zero = TRUE, z = 1.645) {
  if (length(est_cols) == 0) return(rep(0, nrow(df)))

  moe_cols <- str_replace(est_cols, "E$", "M")

  E_mat <- as.matrix(df[, ..est_cols])
  V_mat <- as.matrix(df[, ..moe_cols])

  V_mat[V_mat %in% c(-555555555)] <- NA
  E_mat[E_mat %in% c(-555555555)] <- NA

  V_mat <- (V_mat / z)^2

  if (!one_zero || length(est_cols) == 1) {
    return(rowSums(V_mat, na.rm = TRUE))
  }

  nonzero_sum <- rowSums(V_mat * (E_mat != 0), na.rm = TRUE)

  zero_vars_mat <- V_mat * (E_mat == 0)
  zero_vars_mat[is.na(zero_vars_mat)] <- 0
  zero_max <- exec(pmax, !!!as.data.frame(zero_vars_mat))

  nonzero_sum + zero_max
}

# ── moe_agg ───────────────────────────────────────────────────────────────────
# Aggregation: sum estimates, combine MOEs via sum of variances.
moe_agg <- function(df, parsed, one_zero = TRUE, z = 1.645) {
  est_cols <- paste0(parsed$num_vars, "E")
  est      <- rowSums(as.matrix(df[, ..est_cols]), na.rm = TRUE)
  variance <- compute_var_vector(df, est_cols, one_zero, z)
  moe      <- sqrt(variance) * z
  data.table(est = est, moe = moe)
}

# ── moe_proprat ───────────────────────────────────────────────────────────────
# Proportion or ratio: divide numerator by denominator, combine MOEs correctly.
moe_proprat <- function(df, parsed, method, one_zero = TRUE, z = 1.645) {
  num_cols <- paste0(parsed$num_vars, "E")
  den_cols <- paste0(parsed$den_vars, "E")

  numerator   <- rowSums(as.matrix(df[, ..num_cols]), na.rm = TRUE)
  denominator <- rowSums(as.matrix(df[, ..den_cols]), na.rm = TRUE)

  p <- numerator / denominator

  var_num <- compute_var_vector(df, num_cols, one_zero, z)
  var_den <- compute_var_vector(df, den_cols, one_zero, z)

  if (method %in% c("proportion", "prop")) {
    var_prop  <- (var_num - p^2 * var_den) / (denominator^2)
    var_ratio <- (var_num + p^2 * var_den) / (denominator^2)
    var_final <- ifelse(!is.na(var_prop) & var_prop >= 0, var_prop, var_ratio)
  } else {
    var_final <- (var_num + p^2 * var_den) / (denominator^2)
  }

  est <- p
  moe <- sqrt(pmax(var_final, 0)) * z

  if (parsed$multiply_100) {
    est <- est * 100
    moe <- moe * 100
  }

  est <- ifelse(denominator == 0, NA_real_, est)
  moe <- ifelse(denominator == 0, NA_real_, moe)

  data.table(est = est, moe = moe)
}

# ── compute_one_variable ──────────────────────────────────────────────────────
# Traffic cop: routes each formula to the right MOE function.
compute_one_variable <- function(df, parsed, method, one_zero = TRUE, z = 1.645) {
  method <- tolower(method)

  if (method == "variable") {
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
# needs_state: FALSE means no "in" clause needed (e.g. nation, state)
geo_lookup <- data.table(
  level = c(
    "nation", "state",
    "county", "county.subdivision", "tract", "block.group",
    "place", "school.district.unified", "school.district.elementary",
    "school.district.secondary"
  ),
  for_param = c(
    "us:1", "state:*",
    "county:*", "county subdivision:*", "tract:*", "block group:*",
    "place:*", "unified school district:*",
    "elementary school district:*", "secondary school district:*"
  ),
  needs_state  = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  needs_county = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  geo_id_prefix = c(
    "0100000US", "0400000US",
    "0500000US", "0600000US", "1400000US", "1500000US",
    "1600000US", "9700000US", "9500000US", "9600000US"
  )
)

# ── fetch_acs ─────────────────────────────────────────────────────────────────
# Fetches ACS 5-year estimates from the Census API.
fetch_acs <- function(variables, year, for_param, in_param,
                      api_key = Sys.getenv("CENSUS_API_KEY")) {
  base_url  <- sprintf("https://api.census.gov/data/%d/acs/acs5", year)
  sentinels <- c(-555555555, -666666666, -999999999, -888888888)

  chunks <- split(variables, ceiling(seq_along(variables) / 24))

  fetch_chunk <- function(vars) {
    cols <- c("GEO_ID", "NAME",
              as.vector(rbind(paste0(vars, "E"), paste0(vars, "M"))))

    req <- request(base_url) |>
      req_url_query(
        get = paste(cols, collapse = ","),
        `for` = for_param,
        key = api_key
      )

    # only add "in" clause if in_param is not empty
    if (!is.null(in_param) && nchar(in_param) > 0) {
      req <- req |> req_url_query(`in` = in_param)
    }

    resp <- req |> req_perform()

    ct <- resp_content_type(resp)
    if (!grepl("json", ct)) {
      body <- resp_body_string(resp)
      msg  <- regmatches(body, regexpr("(?<=<p>)[^<]+", body, perl = TRUE))
      stop("Census API error: ", if (length(msg)) trimws(msg) else body)
    }

    json    <- resp_body_json(resp, simplifyVector = FALSE)
    headers <- unlist(json[[1]])
    mat     <- do.call(rbind, lapply(json[-1], unlist))
    dt      <- as.data.table(mat)
    setnames(dt, headers)

    keep <- intersect(names(dt), c("GEO_ID", "NAME", cols))
    dt[, ..keep]
  }

  batches <- lapply(chunks, fetch_chunk)

  out <- Reduce(
    function(a, b) merge(a, b, by = c("GEO_ID", "NAME"), all = TRUE),
    batches
  )

  data_cols <- setdiff(names(out), c("GEO_ID", "NAME"))
  out[, (data_cols) := lapply(.SD, as.numeric), .SDcols = data_cols]
  for (col in data_cols) {
    set(out, which(out[[col]] %in% sentinels), col, NA_real_)
  }

  out
}

# ── extract_fips ──────────────────────────────────────────────────────────────
# Parses the Census GEO_ID column into FIPS component columns.
extract_fips <- function(df) {
  dt <- as.data.table(df)

  prefix  <- substr(dt$GEO_ID[1], 1, 9)
  geo_row <- geo_lookup[geo_id_prefix == prefix]
  if (nrow(geo_row) == 0) stop(paste("Unknown GEO_ID prefix:", prefix))

  lvl    <- geo_row$level
  suffix <- sub(".*US", "", dt$GEO_ID)

  if (lvl == "nation") {
    # no FIPS columns needed

  } else if (lvl == "state") {
    dt[, st := substr(suffix, 1, 2)]

  } else if (lvl == "county") {
    dt[, st   := substr(suffix, 1, 2)]
    dt[, cnty := substr(suffix, 3, 5)]

  } else if (lvl == "county.subdivision") {
    dt[, st     := substr(suffix, 1, 2)]
    dt[, cnty   := substr(suffix, 3, 5)]
    dt[, cousub := substr(suffix, 6, 10)]

  } else if (lvl == "tract") {
    dt[, st    := substr(suffix, 1, 2)]
    dt[, cnty  := substr(suffix, 3, 5)]
    dt[, tract := substr(suffix, 6, 11)]

  } else if (lvl == "block.group") {
    dt[, st    := substr(suffix, 1, 2)]
    dt[, cnty  := substr(suffix, 3, 5)]
    dt[, tract := substr(suffix, 6, 11)]
    dt[, bg    := substr(suffix, 12, 12)]

  } else if (lvl == "place") {
    dt[, st    := substr(suffix, 1, 2)]
    dt[, place := substr(suffix, 3, 7)]

  } else {
    dt[, st   := substr(suffix, 1, 2)]
    dt[, sdid := substr(suffix, 3, 7)]
  }

  dt
}

# ── acsdata ───────────────────────────────────────────────────────────────────
# Fetches all ACS variables for one or more geography levels, with RDS caching.
# state: 2-digit FIPS string or vector of strings e.g. c("55", "36")
#        not required for level = "nation"
acsdata <- function(formulas, level, year, state = NULL,
                    api_key   = Sys.getenv("CENSUS_API_KEY"),
                    cache_dir = tempdir()) {

  all_vars <- unique(unlist(lapply(formulas, getvars)))

  # collapse multiple states into one string for cache key
  state_key <- if (is.null(state)) "all" else paste(sort(state), collapse = "-")

  fetch_one_level <- function(lvl) {
    cache_file <- file.path(
      cache_dir,
      sprintf("acs_%s_%s_%s.rds", lvl, year, state_key)
    )

    if (file.exists(cache_file)) {
      message("Loading from cache: ", cache_file)
      return(readRDS(cache_file))
    }

    geo_row <- geo_lookup[level == lvl]
    if (nrow(geo_row) == 0) stop(paste("Unknown level:", lvl))

    # build in_param based on what the geography needs
    in_param <- if (!geo_row$needs_state) {
      ""  # nation level — no in clause
    } else if (geo_row$needs_county) {
      sprintf("state:%s+county:*", paste(state, collapse = ","))
    } else {
      sprintf("state:%s", paste(state, collapse = ","))
    }

    dt <- fetch_acs(all_vars, year, geo_row$for_param, in_param, api_key)
    dt <- extract_fips(dt)

    saveRDS(dt, cache_file)
    message("Cached to: ", cache_file)
    dt
  }

  setNames(lapply(level, fetch_one_level), level)
}

# ── sumacs ────────────────────────────────────────────────────────────────────
# Main loop: computes derived estimates and MOEs for all spec rows and levels.
sumacs <- function(spec, data, one_zero = TRUE, z = 1.645, file = NULL) {

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

  if (is.data.table(data)) {
    data <- list(data)
  }

  per_level <- lapply(data, compute_level)
  result    <- rbindlist(per_level, fill = TRUE, idcol = "level")

  if (!is.null(file)) fwrite(result, file)
  result
}

# ── read_spec ─────────────────────────────────────────────────────────────────
# Reads and validates the variable spec CSV file.
read_spec <- function(path) {
  sheet <- fread(file = path)

  required <- c("formula", "myfield", "type")
  missing  <- setdiff(required, names(sheet))
  if (length(missing) > 0)
    stop("Spec CSV missing required columns: ", paste(missing, collapse = ", "))

  sheet[, type := fcase(
    tolower(type) %in% c("proportion", "prop"),  "prop",
    tolower(type) %in% c("ratio"),               "ratio",
    tolower(type) %in% c("aggregation", "agg"),  "agg",
    tolower(type) %in% c("variable", "var"),     "variable",
    default = NA_character_
  )]

  bad <- sheet[is.na(type), myfield]
  if (length(bad) > 0)
    warning("Unrecognized method for variable(s): ", paste(bad, collapse = ", "))

  setnames(sheet, c("myfield", "type"), c("varname", "method"))

  sheet
}