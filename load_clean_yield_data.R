setwd("/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts")

# install.packages("rnassqs")   # uncomment if needed
library(rnassqs)
library(dplyr)
library(stringr)
library(purrr)

# 1) Authenticate once per session (get a free key at: quickstats.nass.usda.gov)
# Either set env var beforehand: Sys.setenv(NASSQS_TOKEN="YOUR_KEY")
nassqs_auth(key = Sys.getenv("NASSQS_TOKEN"))

# 2) Map your crop names to Quick Stats commodity descriptors
commod_map <- tibble::tribble(
  ~name,        ~commodity_desc,         ~class_desc,
  "sugarbeets", "SUGARBEETS",             "ALL CLASSES",
  "sunflower",  "SUNFLOWER",              "ALL CLASSES",
  "safflower",  "SAFFLOWER",              "ALL CLASSES",
  "canola",     "CANOLA",                 "ALL CLASSES",
  "dry beans",  "BEANS, DRY EDIBLE",      "ALL CLASSES",
  "green peas", "PEAS, GREEN",            "ALL CLASSES"
)

# 3) Helper to fetch yields
get_yield <- function(commodity_desc,
                      class_desc = "ALL CLASSES",
                      state = NULL,            # e.g. "NORTH DAKOTA"; NULL = all states
                      agg_level = "STATE",     # "STATE" or "COUNTY"
                      year_start = 1980,
                      year_end   = as.numeric(format(Sys.Date(), "%Y"))) {
  
  params <- list(
    source_desc         = "SURVEY",
    sector_desc         = "CROPS",
    statisticcat_desc   = "YIELD",
    commodity_desc      = commodity_desc,
    class_desc          = class_desc,
    prodn_practice_desc = "ALL PRODUCTION PRACTICES",
    util_practice_desc  = "ALL UTILIZATION PRACTICES",
    agg_level_desc      = toupper(agg_level),
    year__GE            = year_start,
    year__LE            = year_end
  )
  if (!is.null(state)) params$state_name <- toupper(state)
  
  out <- nassqs(params)
  
  # Clean common fields
  out %>%
    transmute(
      year        = as.integer(year),
      state       = dplyr::coalesce(state_name, state_alpha),
      county      = dplyr::coalesce(county_name, NA_character_),
      commodity   = commodity_desc,
      class       = class_desc,
      unit        = unit_desc,                      # e.g., "TONS / ACRE", "LB / ACRE", "CWT / ACRE", "BU / ACRE"
      value       = as.numeric(str_replace_all(Value, ",", "")),
      agg_level   = agg_level_desc
    ) %>%
    arrange(year, state, county)
}

# 4) Pull all six crops at STATE level (all states), 1980–present
yield_state_all <- map_dfr(commod_map$commodity_desc, ~get_yield(.x))

# 5) Example: restrict to a state (North Dakota) and bind the six crops
yield_nd <- map2_dfr(commod_map$commodity_desc,
                     commod_map$class_desc,
                     ~get_yield(.x, class_desc = .y, state = "NORTH DAKOTA"))

# 6) Example: COUNTY level for one crop (be mindful of row counts/rate limits)
yield_nd_canola_county <- get_yield("CANOLA", state = "NORTH DAKOTA", agg_level = "COUNTY", year_start = 2000)

# --- Optional: quick pivots or checks ---------------------------------------
# Units vary by crop (don’t assume bushels). Quick glance:
dplyr::count(yield_state_all, commodity, unit) %>% arrange(commodity, unit)


