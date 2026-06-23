# acsest

A self-contained R script for computing derived American Community Survey (ACS)
estimates and margins of error (MOEs) from a batch specification file.

Built at the [Applied Population Laboratory](https://apl.wisc.edu/), UW–Madison
to replace the archived [`acsr`](https://github.com/sdaza/acsr) package.

---

## What it does

- Reads a CSV file of ACS variable formulas (proportions, aggregations, ratios)
- Downloads raw estimates and MOEs from the Census Bureau REST API
- Computes derived statistics with correct MOE propagation
- Implements the Census Bureau one-zero variance correction
- Supports multiple geography levels in one run (county, tract, block group, school districts, etc.)
- Outputs a wide CSV with `{varname}_est` and `{varname}_moe` columns

---

## Quick start

```r
# Set your Census API key (get one at https://api.census.gov/data/key_signup.html)
Sys.setenv(CENSUS_API_KEY = "your_key_here")

# Load the script
source("acsest.R")

# Read your spec file
sheet <- read_spec("RR_vars_22ACS.csv")

# Download raw ACS data
wi_acs <- acsdata(sheet[, formula], level = "county", year = 2022, state = "55")

# Compute derived estimates and MOEs
out <- sumacs(spec = sheet, data = wi_acs, file = "output.csv")
```

---

## Spec CSV format

One row per derived variable. Required columns:

| Column | Description | Example |
|--------|-------------|---------|
| `formula` | ACS variable codes with arithmetic operators | `(B17001_004 + B17001_018) / B17001_001 * 100` |
| `myfield` | Output column name | `poverty_rate` |
| `type` | `variable`, `agg`, `prop`, or `ratio` (case-insensitive) | `prop` |

---

## Supported geographies

| Level string | Geography |
|-------------|-----------|
| `county` | Counties |
| `tract` | Census tracts |
| `block.group` | Block groups |
| `county.subdivision` | County subdivisions |
| `place` | Places |
| `school.district.unified` | Unified school districts |
| `school.district.elementary` | Elementary school districts |
| `school.district.secondary` | Secondary school districts |

---

## Dependencies

```r
install.packages(c("data.table", "httr2", "stringr", "rlang"))
```

---

## API key setup

Option 1 — set in your R session:
```r
Sys.setenv(CENSUS_API_KEY = "your_key_here")
```

Option 2 — create a `.env` file in the same directory as `acsest.R`:
```
CENSUS_API_KEY=your_key_here
```

The script loads it automatically on startup.

---

## Differences from ACSR

| | ACSR | acsest.R |
|-|------|----------|
| API dependency | `acs` (archived) | `httr2` (maintained) |
| Year argument | `endyear=` | `year=` |
| State argument | Implicit via `geo.make()` | Explicit `state="55"` |
| Spec loading | Manual column extraction | `read_spec()` |
| One-zero correction | `one.zero=TRUE` | `one_zero=TRUE` |

---

## Projects using this script

- [UW Food Security Mapping](https://foodsecurity.wisc.edu/)
- [APL Risk & Reach](https://riskandreach.apl.wisc.edu/)