
# ACS API key is needed
if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("CENSUS_API_KEY environmental variable is not set")
}

# functions
