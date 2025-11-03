# ============================================================
# Two-season CV using MONTHLY temperature exposure bins
# (data.frame-only; sanity CV first, then full grid)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)  # only for fread/setnames; we convert to data.frame immediately
  library(fixest)
  library(ggplot2)
})

# -----------------------------
# PATHS / SETTINGS
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"

yield_csv <- file.path(root, "sugarbeet_yields.csv")

# Quiet fixest notes
if ("setFixest_notes" %in% ls(getNamespace("fixest"))) fixest::setFixest_notes(FALSE) else options(fixest_notes = FALSE)

# -----------------------------
# LOAD & STANDARDIZE (then force data.frame)
# -----------------------------
tempbinsdata_dt <- fread(bins_csv)
weatherdata_dt  <- fread(weath_csv)
yld_dt          <- fread(yield_csv)
# --- Normalize sugarbeet column names (spaces -> dots), then build fips ---
# Run this RIGHT AFTER you read yld_dt

# 1) Rename if needed
if ("State ANSI"  %in% names(yld_dt)) data.table::setnames(yld_dt, "State ANSI",  "State.ANSI")
if ("County ANSI" %in% names(yld_dt)) data.table::setnames(yld_dt, "County ANSI", "County.ANSI")
if ("Year"        %in% names(yld_dt) && !("year" %in% names(yld_dt))) data.table::setnames(yld_dt, "Year", "year")
if ("VALUE"       %in% names(yld_dt)) data.table::setnames(yld_dt, "VALUE", "Value")
if ("value"       %in% names(yld_dt)) data.table::setnames(yld_dt, "value", "Value")

# 2) Make Value numeric (idempotent) and drop rollups w/ missing county ANSI
yld_dt[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
yld_dt <- yld_dt[is.finite(Value) & Value > 0]

# 3) Build 5-digit FIPS and keep all states
stopifnot(all(c("State.ANSI","County.ANSI") %in% names(yld_dt)))
yld_dt[, fips := 1000L * as.integer(State.ANSI) + as.integer(County.ANSI)]
yld_dt <- yld_dt[is.finite(fips)]
yld_dt[, fips := sprintf("%05d", fips)]
yld_dt[, year := as.integer(year)]
yld_dt <- yld_dt[year %in% years_use]
# --- (optional) if years_use wasn't defined earlier, fall back to data span
if (!exists("years_use")) years_use <- sort(unique(yld_dt$year))

# 5) Map (only states with sugarbeets; no county outlines)
suppressPackageStartupMessages({ library(sf); library(tigris); library(ggplot2); library(scales); library(dplyr) })
options(tigris_use_cache = TRUE, tigris_class = "sf")

# Geometries
cnty_all <- tigris::counties(year = 2020, cb = TRUE) |>
  sf::st_transform(5070) |>
  mutate(fips = GEOID)

states_all <- tigris::states(year = 2020, cb = TRUE) |>
  sf::st_transform(5070)

# Join county geometries to sugarbeet averages
mdat <- cnty_all |> left_join(beet_avg, by = "fips")

# Identify states that actually have sugarbeet observations (non-NA avg_yield)
states_with_beets <- mdat |>
  filter(is.finite(avg_yield)) |>
  st_drop_geometry() |>
  distinct(STATEFP) |>
  pull(STATEFP)

# Filter to those states (both counties and state borders)
cnty   <- mdat       |> filter(STATEFP %in% states_with_beets)
states <- states_all |> filter(STATEFP %in% states_with_beets)

# (Optional) crop to the union of those states (removes oceans / AK/HI if absent)
clip_poly <- st_union(states)
cnty  <- st_intersection(cnty,  clip_poly) |> st_make_valid()
# states <- st_intersection(states, clip_poly)  # usually identical to 'states'

# Nice breaks & map (NO county outlines)
brks <- pretty(range(cnty$avg_yield, na.rm = TRUE), n = 6)

p_beet_all <- ggplot() +
  geom_sf(data = cnty,   aes(fill = avg_yield), color = NA) +         # <- no outlines
  geom_sf(data = states, fill = NA, color = "grey25", size = 0.25) +  # keep state borders
  coord_sf(crs = st_crs(5070)) +
  scale_fill_viridis_c(name = "tons/acre",
                       labels = comma,
                       breaks = brks,
                       limits = range(brks, na.rm = TRUE)) +
  labs(title    = "Average Sugarbeet Yield by County, 1981–2022",
       subtitle = "Only states with observed sugarbeet yields are shown",
       caption  = "Arithmetic mean of NASS county yields (tons/acre)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11, margin = margin(b = 6)))

print(p_beet_all)
# ggsave("sugarbeet_avg_yield_1981_2022_states-with-beets.png", p_beet_all, width = 13, height = 8.5, dpi = 350)



## post 2010 to ensure focus on currently relevant locations:
# ================================
# 0) Restrict to post-2010
# ================================
years_use <- sort(unique(yld_dt$year))   # if not already set
years_use <- years_use[years_use >= 2011]
yld_dt     <- yld_dt[year %in% years_use]

# ================================
# 1) County average yield (2011+)
# ================================
# arithmetic mean on original scale (Value should be tons/acre)
beet_avg <- yld_dt[, .(
  avg_yield = mean(Value, na.rm = TRUE),
  n_years   = uniqueN(year)
), by = .(fips)]

# ================================
# 2) State ranking (choose one)
#    (A) With acreage -> rank by avg annual PRODUCTION (tons)
#    (B) Without acreage -> rank by avg annual county YIELD
# ================================

# ---- (A) If you have acreage, compute production shares ----
ac_csv <- file.path(root, "sugarbeet_acreage.csv")
if (file.exists(ac_csv)) {
  ac_dt <- data.table::fread(ac_csv)
  # normalize names
  setnames(ac_dt,
           old = c("State ANSI","County ANSI","Year","VALUE","value"),
           new = c("State.ANSI","County.ANSI","year","Value","Value"),
           skip_absent = TRUE)
  ac_dt[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
  ac_dt <- ac_dt[is.finite(Value) & Value > 0]
  ac_dt[, `:=`(
    year = as.integer(year),
    fips = sprintf("%05d", 1000L*as.integer(State.ANSI) + as.integer(County.ANSI))
  )]
  ac_dt <- ac_dt[year %in% years_use]
  
  # county-year join: acres × yield (tons/acre) -> tons
  # yld_dt already restricted to post-2010
  cy <- merge(ac_dt[, .(fips, year, acres = Value, State.ANSI)],
              yld_dt[, .(fips, year, yld_tpa = Value, State.ANSI)],
              by = c("fips","year"),
              all = FALSE)
  # if State.ANSI missing on either side, fill from the non-missing
  cy[is.na(State.ANSI.x) & !is.na(State.ANSI.y), State.ANSI.x := State.ANSI.y]
  cy[is.na(State.ANSI.y) & !is.na(State.ANSI.x), State.ANSI.y := State.ANSI.x]
  cy[, State.ANSI := fifelse(is.na(State.ANSI.x), State.ANSI.y, State.ANSI.x)]
  cy[, statefp := sprintf("%02d", as.integer(State.ANSI))]
  
  cy[, prod_tons := acres * yld_tpa]
  
  # state-year totals, then average across 2011+ years
  st_year <- cy[, .(state_prod_tons = sum(prod_tons, na.rm = TRUE)), by = .(statefp, year)]
  st_avg  <- st_year[, .(avg_annual_prod_tons = mean(state_prod_tons, na.rm = TRUE)), by = statefp]
  
  # shares + ranking
  st_avg[, total := sum(avg_annual_prod_tons, na.rm = TRUE)]
  st_avg[, share := avg_annual_prod_tons / total]
  setorder(st_avg, -share)
  st_avg[, cum_share := cumsum(share)]
  
  cat("\n=== State ranking by average annual PRODUCTION (tons), 2011+ ===\n")
  print(st_avg[, .(statefp, avg_annual_prod_tons, share, cum_share)])
} else {
  # ---- (B) No acreage: rank states by average annual county YIELD ----
  # build statefp from fips (first two digits)
  yld_dt[, statefp := substr(fips, 1, 2)]
  st_avg_yield <- yld_dt[, .(avg_yield = mean(Value, na.rm = TRUE)), by = statefp]
  setorder(st_avg_yield, -avg_yield)
  cat("\n=== State ranking by average YIELD (tons/acre), 2011+ ===\n")
  print(st_avg_yield)
}

# ================================
# 3) Map (2011+ only; only states with data; no county outlines)
# ================================
suppressPackageStartupMessages({ library(sf); library(tigris); library(ggplot2); library(scales); library(dplyr) })
options(tigris_use_cache = TRUE, tigris_class = "sf")

cnty_all <- tigris::counties(year = 2020, cb = TRUE) |>
  sf::st_transform(5070) |>
  dplyr::mutate(fips = GEOID)

states_all <- tigris::states(year = 2020, cb = TRUE) |>
  sf::st_transform(5070)

mdat <- cnty_all |>
  dplyr::left_join(beet_avg, by = "fips")

states_with_beets <- mdat |>
  dplyr::filter(is.finite(avg_yield)) |>
  sf::st_drop_geometry() |>
  dplyr::distinct(STATEFP) |>
  dplyr::pull(STATEFP)

cnty   <- mdat       |> dplyr::filter(STATEFP %in% states_with_beets)
states <- states_all |> dplyr::filter(STATEFP %in% states_with_beets)

clip_poly <- sf::st_union(states)
cnty  <- sf::st_intersection(cnty, clip_poly) |> sf::st_make_valid()

brks <- pretty(range(cnty$avg_yield, na.rm = TRUE), n = 6)

p_beet_post2010 <- ggplot() +
  geom_sf(data = cnty,   aes(fill = avg_yield), color = NA) +
  geom_sf(data = states, fill = NA, color = "grey25", size = 0.25) +
  coord_sf(crs = sf::st_crs(5070)) +
  scale_fill_viridis_c(name = "tons/acre",
                       labels = scales::comma,
                       breaks = brks,
                       limits = range(brks, na.rm = TRUE)) +
  labs(
    title    = "Average Sugarbeet Yield by County, 2011–2022",
    subtitle = "Only states with observed sugarbeet yields (post-2010)",
    caption  = "Arithmetic mean of NASS county yields (tons/acre); 2011–2022"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11, margin = margin(b = 6)))

print(p_beet_post2010)
