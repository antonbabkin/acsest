# ── Test 1: moe_agg ───────────────────────────────────────────────────────────
# MOEs 10, 8, 6 → sqrt(10² + 8² + 6²) = 14.14
df_test <- data.table(
  AE = 100, AM = 10,
  BE = 50,  BM = 8,
  CE = 25,  CM = 6
)
parsed_test <- parse_formula("A + B + C")
result      <- moe_agg(df_test, parsed_test)
stopifnot(round(result$moe, 2) == 14.14)
message("moe_agg: PASS")

# ── Test 2: moe_proprat ───────────────────────────────────────────────────────
# numerator 50 (MOE 5), denominator 200 (MOE 10) → est 0.25, moe ≈ 0.046
df_test2 <- data.table(
  AE = 50,  AM = 5,
  BE = 200, BM = 10
)
parsed_test2 <- parse_formula("A / B")
result2      <- moe_proprat(df_test2, parsed_test2, method = "prop")
stopifnot(round(result2$est, 2) == 0.25)
stopifnot(round(result2$moe, 3) == 0.022)
message("moe_proprat: PASS")

# ── Test 3: compute_one_variable routing ─────────────────────────────────────
result3 <- compute_one_variable(df_test,  parse_formula("A + B + C"), "agg")
result4 <- compute_one_variable(df_test2, parse_formula("A / B"),     "prop")
result5 <- compute_one_variable(df_test,  parse_formula("A"),         "variable")
stopifnot(!is.null(result3$moe))
stopifnot(!is.null(result4$moe))
stopifnot(result5$est == 100)
message("compute_one_variable routing: PASS")

# ── Test 4: extract_fips ──────────────────────────────────────────────────────
df_geo <- data.table(GEO_ID = "0500000US55025", NAME = "Dane County, WI")
result6 <- extract_fips(df_geo)
stopifnot(result6$st   == "55")
stopifnot(result6$cnty == "025")
message("extract_fips: PASS")



# -- BLOCK 2 -----------------------------------------------------------------------
# ── Test 5: fetch_acs live ────────────────────────────────────────────────────
raw <- fetch_acs("B17001_001", year = 2022, 
                 for_param = "county:*", in_param = "state:55")
stopifnot(nrow(raw) == 72)   # 72 WI counties
stopifnot("B17001_001E" %in% names(raw))
message("fetch_acs: PASS — ", nrow(raw), " counties returned")

# ── Test 6: acsdata multi-level ───────────────────────────────────────────────
wi_acs <- acsdata("B17001_002 / B17001_001", 
                  level = c("county", "tract"), year = 2022, state = "55")
stopifnot(length(wi_acs) == 2)
stopifnot("county" %in% names(wi_acs))
stopifnot("tract"  %in% names(wi_acs))
message("acsdata multi-level: PASS")

# ── Test 7: end-to-end sumacs ─────────────────────────────────────────────────
spec_test <- data.table(
  varname = "pct_poverty",
  formula = "B17001_002 / B17001_001",
  method  = "prop"
)
out <- sumacs(spec = spec_test, data = wi_acs)
stopifnot("pct_poverty_est" %in% names(out))
stopifnot("pct_poverty_moe" %in% names(out))
stopifnot("level" %in% names(out))
message("sumacs end-to-end: PASS — ", nrow(out), " rows")
print(out[level == "county"][order(-pct_poverty_est)][1:5])  # top 5 highest poverty counties