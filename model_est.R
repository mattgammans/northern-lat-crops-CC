# ============================================================
# Northern-lat crops: GS selection, Bin model, Avg-T/P model,
# and Merel–Gammans climate penalty (squared weather–climate gaps)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fixest)   # FE estimation + clustering
  library(purrr)
  library(stringr)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
years_use        <- 1981:2022
gs_candidates_s  <- 4:7     # candidate GS starts (Apr–Jul)
gs_candidates_e  <- 8:10    # candidate GS ends   (Aug–Oct)
climate_baseline <- 1981:2018 # normals window for MG penalty (pick explicitly)
# Expecting NASS-style columns: Year, Value, State.ANSI, County.ANSI, (maybe fips)
yld <- fread("sugarbeet_yields.csv")
yld$State.ANSI = yld$"State ANSI"
yld$County.ANSI = yld$"County ANSI"
yld[, Value := as.numeric(gsub("[^0-9.\\-]", "", Value))]
yld$yield = log(as.numeric(yld$Value))
yld$fips = 1000*as.integer(yld$"State.ANSI") + as.integer(yld$"County.ANSI")
yld <- yld[yld$State %in% c("MICHIGAN","MINNESOTA","MONTANA","NEBRASKA","NORTH DAKOTA"),]
# Year as numeric
if (!("year" %in% names(yld))) {
  if ("Year" %in% names(yld)) yld[, year := as.integer(Year)] else stop("No 'year' or 'Year' column in yield file.")
} else {
  yld[, year := as.integer(year)]
}

# Keep analysis years overlapping weather
yld <- yld[year %in% years]

# -----------------------------
# HELPERS
# -----------------------------
mk_weather_agg <- function(weatherdata, tempbinsdata, years, gstart, gstop) {
  stopifnot(gstop > gstart)
  gseas <- gstart:gstop
  
  # Precip sum over GS
  precip <- aggregate(
    weatherdata[weatherdata$month %in% gseas & weatherdata$year %in% years, "ppt"],
    list(
      weatherdata[weatherdata$month %in% gseas & weatherdata$year %in% years, ]$fips,
      weatherdata[weatherdata$month %in% gseas & weatherdata$year %in% years, ]$year
    ),
    sum
  )
  colnames(precip) <- c("fips","year","ppt")
  
  # Temperature bins (sum hours over GS)
  Bins <- c(paste0("bin_", -15:-1), paste0("bin_", 0:50))
  Bins <- Bins[Bins %in% names(tempbinsdata)]
  
  tempbins <- aggregate(
    tempbinsdata[tempbinsdata$month %in% gseas & tempbinsdata$year %in% years, Bins, drop = FALSE],
    list(
      tempbinsdata[tempbinsdata$month %in% gseas & tempbinsdata$year %in% years, ]$fips,
      tempbinsdata[tempbinsdata$month %in% gseas & tempbinsdata$year %in% years, ]$year
    ),
    sum
  )
  colnames(tempbins) <- c("fips","year",Bins)
  
  dw <- merge(precip, tempbins, by = c("fips","year"))
  dw <- dw[order(dw$year, dw$fips), ]
  
  # 5°C bins (shares of hours; divide by 24 to match your convention)
  dw$bin5_freeze <- rowSums(dw[, paste0("bin_", -15:-1), drop = FALSE]) / 24
  mkb <- function(df, L, U) rowSums(df[, paste0("bin_", L:U), drop = FALSE]) / 24
  for (L in seq(0, 30, by = 5)) {
    nm <- sprintf("bin5_%d_%d", L, L+5)
    dw[[nm]] <- mkb(dw, L, L+4)
  }
  dw$bin5_30_35 <- mkb(dw, 30, 34)
  hi <- intersect(paste0("bin_", 35:50), names(dw))
  dw$bin5_35_50 <- if (length(hi)) rowSums(dw[, hi, drop = FALSE]) / 24 else 0
  
  # Precip scaling and square
  dw$ppt  <- dw$ppt / 100
  dw$ppt2 <- dw$ppt^2
  
  # Average GS temperature from 5°C bins (approx.; freeze bin uses -5 °C midpoint)
  bin5_names <- c("bin5_freeze","bin5_0_5","bin5_5_10","bin5_10_15","bin5_15_20",
                  "bin5_20_25","bin5_25_30","bin5_30_35","bin5_35_50")
  mids <- c(-5, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 40)  # last is pooled tail midpoint
  tot_hours <- 24 * rowSums(dw[, bin5_names, drop = FALSE])
  w_hours   <- sweep(24*dw[, bin5_names, drop = FALSE], 2, mids, `*`)
  gs_tbar   <- rowSums(w_hours) / pmax(tot_hours, 1e-9)
  dw$gs_tavg <- as.numeric(gs_tbar)
  
  dw$gs_start <- gstart
  dw$gs_stop  <- gstop
  dw
}

build_df <- function(yld, dw, years) {
  y <- data.table(yld)[year %in% years]
  df <- merge(y, dw, by = c("fips","year"))
  if (!"yield" %in% names(df)) df$yield <- log(df$Value)
  if (!"State.ANSI" %in% names(df) && "State ANSI" %in% names(df)) {
    df$State.ANSI <- df$`State ANSI`
  }
  df$fips <- as.integer(df$fips)
  df$year <- as.integer(df$year)
  as.data.frame(df)
}

rhs_bins <- function(df) {
  omit <- c("bin5_10_15","bin5_15_20")
  rhs <- c("bin5_freeze","bin5_0_5","bin5_5_10","bin5_10_15","bin5_15_20",
           "bin5_20_25","bin5_25_30","bin5_30_35","bin5_35_50","ppt","ppt2")
  setdiff(rhs, omit)
}

# LOYO CV R^2 for a given df + RHS
# --- replace this whole function ---
cv_r2 <- function(df, rhs) {
  yrs <- sort(unique(df$year))
  r2s <- c()
  for (yv in yrs) {
    tr <- df[df$year != yv, , drop = FALSE]
    te <- df[df$year == yv, , drop = FALSE]
    
    # drop zero-variance regressors within the training fold
    vrs  <- sapply(tr[, rhs, drop = FALSE], function(x) stats::var(x, na.rm = TRUE))
    keep <- names(vrs)[is.finite(vrs) & vrs > 0]
    if (!length(keep)) next
    
    fml <- as.formula(paste0("yield ~ ", paste(keep, collapse = " + "),
                             " | fips + year"))
    fit <- tryCatch(fixest::feols(fml, data = tr), error = function(e) NULL)
    if (is.null(fit)) next
    
    yhat <- predict(fit, newdata = te)
    r2   <- 1 - sum((te$yield - yhat)^2, na.rm = TRUE) /
      sum((te$yield - mean(tr$yield))^2, na.rm = TRUE)
    r2s <- c(r2s, r2)
  }
  mean(r2s, na.rm = TRUE)
}

# Select GS by maximizing LOYO CV R^2 of the bin model
select_gs <- function(weatherdata, tempbinsdata, yld, years, gstarts, gstops) {
  cand <- expand.grid(gsstart = gstarts, gsstop = gstops)
  cand <- cand[cand$gsstop > cand$gsstart, ]
  scores <- purrr::pmap_dbl(cand, function(gsstart, gsstop) {
    dw  <- mk_weather_agg(weatherdata, tempbinsdata, years, gsstart, gsstop)
    df0 <- build_df(yld, dw, years)
    rhs <- rhs_bins(df0)
    cv_r2(df0, rhs)
  })
  best <- which.max(scores)
  list(gsstart = cand$gsstart[best], gsstop = cand$gsstop[best], cv_R2 = scores[best])
}

# -----------------------------
# 1) Endogenous GS selection
# -----------------------------
sel <- select_gs(weatherdata, tempbinsdata, yld, years_use, gs_candidates_s, gs_candidates_e)
cat(sprintf("Selected GS: %s–%s (CV R^2 = %.3f)\n", sel$gsstart, sel$gsstop, sel$cv_R2))

# Build final weather aggregates + merged DF with chosen GS
dataweather <- mk_weather_agg(weatherdata, tempbinsdata, years_use, sel$gsstart, sel$gsstop)
df          <- build_df(yld, dataweather, years_use)

# -----------------------------
# 2) Temperature-bin model
# -----------------------------
# (2) Temperature-bin model  --------
fml_bins <- yield ~ bin5_freeze + bin5_0_5 + bin5_5_10 +
  bin5_20_25 + bin5_25_30 + bin5_30_35 + bin5_35_50 +
  ppt + ppt2 | fips + year
mod_bins <- fixest::feols(fml_bins, data = df, cluster = ~ State.ANSI)
print(summary(mod_bins))

# -----------------------------
# 3) Average GS temperature + precip model
#     (yield ~ gs_tavg + ppt + ppt^2 with FE)
# -----------------------------
# (3) Avg-GS temp + precip ----------
fml_avg  <- yield ~ gs_tavg + ppt + ppt2 | fips + year
mod_avg  <- fixest::feols(fml_avg,  data = df,    cluster = ~ State.ANSI)

cat("\n=== Avg-T + Precip model ===\n")
print(summary(mod_avg))

# -----------------------------
# 4) Merel–Gammans climate penalty
#     Add squared deviations from county climate normals:
#     pen_tempL2 = sum_b (s_b - \bar{s}_b)^2  over 5°C bins
#     pen_ppt    = (ppt - \bar{ppt})^2
# -----------------------------
# Compute county normals on the SAME GS aggregation
# (re-aggregate for baseline window to avoid mismatches)
dw_norm <- mk_weather_agg(weatherdata, tempbinsdata,
                          years = climate_baseline, 
                          gstart = sel$gsstart, gstop = sel$gsstop)

bin_cols <- grep("^bin5_(freeze|\\d+_\\d+)$", names(dw_norm), value = TRUE)

# --- MG normals on the SAME GS; compute climate averages, not bin L2 ---
dw_norm <- mk_weather_agg(
  weatherdata, tempbinsdata,
  years  = climate_baseline,           # e.g., 1981–2010
  gstart = sel$gsstart, gstop = sel$gsstop
)

# county climate normals: GS-average temperature and GS precip
clim_norms <- dw_norm |>
  dplyr::group_by(fips) |>
  dplyr::summarise(
    clim_gs_tavg = mean(gs_tavg, na.rm = TRUE),
    clim_ppt     = mean(ppt,    na.rm = TRUE)
  ) |>
  as.data.frame()

df_mg <- merge(df, clim_norms, by = "fips", all.x = TRUE)

# --- Merel–Gammans penalties: squared deviations from climate averages ---
df_mg$pen_temp <- (df_mg$gs_tavg - df_mg$clim_gs_tavg)^2
df_mg$pen_ppt  <- (df_mg$ppt    - df_mg$clim_ppt    )^2

# penalties
df_mg$gs_tavg2 <- df_mg$gs_tavg^2



# MG model = avg-T/P model + penalties
fml_mg <- yield ~ gs_tavg +gs_tavg2+ ppt  + ppt2 +pen_temp + pen_ppt +year | fips   
mod_mg <- fixest::feols(fml_mg, data = df_mg,  cluster = ~ State.ANSI)

cat("\n=== Avg-T + Precip with Merel–Gammans penalty ===\n")
print(summary(mod_mg))

# -----------------------------
# Niceties: compact coefficient tables (optional)
# -----------------------------
grab <- function(m) broom::tidy(m, conf.int = TRUE)
suppressWarnings({
  if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
})
res_tables <- list(
  selected_gs = sel,
  bins_model  = grab(mod_bins),
  avg_model   = grab(mod_avg),
  mg_model    = grab(mod_mg)
)

# Example: write CSVs (optional)
# write.csv(res_tables$bins_model, "out_bins_model.csv", row.names = FALSE)
# write.csv(res_tables$avg_model,  "out_avg_model.csv",  row.names = FALSE)
# write.csv(res_tables$mg_model,   "out_mg_model.csv",   row.names = FALSE)

# -----------------------------
# Quick sanity prints
# -----------------------------
cat("\nSelected GS window:", sel$gsstart, "to", sel$gsstop, "\n")
cat("Obs used:", nrow(df), "  Counties:", dplyr::n_distinct(df$fips), "  Years:", dplyr::n_distinct(df$year), "\n")
