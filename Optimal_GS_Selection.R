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

bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")
yield_csv <- file.path(root, "sunflower_yields_1980_2024.csv")

years_use <- 1981:2022

# Thinned search grid (Apr..Oct)
s1_start_grid <- c(5)
s2_start_grid <- c(7)
s2_end_grid   <- c(9, 10)

gdd_bases_grid <- c(5)            # per-season GDD base (°C)
kdd_thr_grid   <- c(28, 30, 32)   # per-season KDD threshold (°C)
gdd_cap_grid   <- kdd_thr_grid    # <-- GDD caps we will precompute (cap = KDD thr)
use_centers    <- FALSE           # bin centers vs integer bins

# Quiet fixest notes
if ("setFixest_notes" %in% ls(getNamespace("fixest"))) fixest::setFixest_notes(FALSE) else options(fixest_notes = FALSE)

# -----------------------------
# LOAD & STANDARDIZE (then force data.frame)
# -----------------------------
tempbinsdata_dt <- fread(bins_csv)
weatherdata_dt  <- fread(weath_csv)
yld_dt          <- fread(yield_csv)

# Yield column canonicalization
data.table::setnames(
  yld_dt,
  old = c("State ANSI","County ANSI","Year","VALUE","value"),
  new = c("State.ANSI","County.ANSI","year","Value","Value"),
  skip_absent = TRUE
)

stopifnot(all(c("fips","year","month") %in% names(tempbinsdata_dt)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata_dt)))

# Types
tempbinsdata_dt[, `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
weatherdata_dt[,  `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]

# Precip column -> ppt
ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata_dt))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
data.table::setnames(weatherdata_dt, ppt_col[1], "ppt")

# Filter window (keep all months; seasons will pick 5..9/10)
tempbinsdata_dt <- tempbinsdata_dt[year %in% years_use & month >= 1 & month <= 12]
weatherdata_dt  <- weatherdata_dt[ year %in% years_use & month >= 1 & month <= 12]

# Yields: clean & log
yld_dt[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
yld_dt <- yld_dt[is.finite(Value) & Value > 0]
yld_dt[, `:=`(
  year  = as.integer(year),
  fips  = 1000L * as.integer(State.ANSI) + as.integer(County.ANSI),
  yield = log(Value)
)]
yld_dt <- yld_dt[year %in% years_use]

# ---- Force data.frame versions (ALL modeling uses these) ----
tempbinsdata <- as.data.frame(tempbinsdata_dt); rm(tempbinsdata_dt)
weatherdata  <- as.data.frame(weatherdata_dt);  rm(weatherdata_dt)
yld_raw      <- as.data.frame(yld_dt);          rm(yld_dt)

# Ensure ppt is numeric and non-NA
weatherdata$ppt <- as.numeric(weatherdata$ppt)
weatherdata$ppt[is.na(weatherdata$ppt)] <- 0

# Base keys (data.frame)
base_keys_df <- unique(yld_raw[c("fips","year","State.ANSI")])

# -----------------------------
# BIN DISCOVERY & MONTHLY DEGREE-DAYS
# -----------------------------
discover_bins <- function(TB, use_centers = FALSE){
  TBdf <- as.data.frame(TB)
  bins <- grep("^bin_-?\\d+$", names(TBdf), value = TRUE)
  if (!length(bins)) stop("No bin_* columns found. Ensure tempbinsdata is passed.")
  k    <- as.integer(sub("^bin_", "", bins))
  ord  <- order(k); bins <- bins[ord]; k <- k[ord]
  k_rep <- if (use_centers) ifelse(k >= 0, k + 0.5, k - 0.5) else k
  list(TBdf = TBdf, bins = bins, k = k, k_rep = k_rep)
}

monthly_gdd_from_bins <- function(TB, base, cap = 30, use_centers = FALSE){
  info <- discover_bins(TB, use_centers)
  DF   <- info$TBdf[, c("fips","year","month", info$bins), drop = FALSE]
  H    <- as.matrix(DF[, info$bins, drop = FALSE])  # hours
  days <- H / 24
  w    <- pmax(pmin(info$k_rep, cap) - base, 0)
  out  <- DF[, c("fips","year","month")]
  out$gdd_m <- as.numeric(days %*% w)
  out
}

monthly_kdd_from_bins <- function(TB, thr, use_centers = FALSE){
  info <- discover_bins(TB, use_centers)
  DF   <- info$TBdf[, c("fips","year","month", info$bins), drop = FALSE]
  H    <- as.matrix(DF[, info$bins, drop = FALSE])
  days <- H / 24
  w    <- pmax(info$k_rep - thr, 0)
  out  <- DF[, c("fips","year","month")]
  out$kdd_m <- as.numeric(days %*% w)
  out
}

# Precompute monthly KDD (per thr)
message("Precomputing monthly GDD/KDD tables ...")
kdd_map <- setNames(
  lapply(kdd_thr_grid, function(h) monthly_kdd_from_bins(tempbinsdata, thr = h, use_centers = use_centers)),
  paste0("h", kdd_thr_grid)
)

# Precompute monthly GDD for every (base, cap) where cap ∈ kdd_thr_grid
gdd_keys <- expand.grid(base = gdd_bases_grid, cap = gdd_cap_grid, KEEP.OUT.ATTRS = FALSE)
gdd_map  <- list()
for (i in seq_len(nrow(gdd_keys))) {
  b <- gdd_keys$base[i]; c <- gdd_keys$cap[i]
  key <- sprintf("b%d_c%d", b, c)
  gdd_map[[key]] <- monthly_gdd_from_bins(tempbinsdata, base = b, cap = c, use_centers = use_centers)
}

# -----------------------------
# BUILD SEASON FEATURES (cap = KDD thr per season)
# -----------------------------
build_2season_features_bins <- function(s1_start, s2_start, s2_end,
                                        gdd1_base, kdd1_thr,
                                        gdd2_base, kdd2_thr,
                                        gdd_map, kdd_map,
                                        weather_monthly){
  
  s1_months <- s1_start:(s2_start - 1)
  s2_months <- s2_start:s2_end
  
  key1 <- sprintf("b%d_c%d", gdd1_base, kdd1_thr)
  key2 <- sprintf("b%d_c%d", gdd2_base, kdd2_thr)
  
  if (is.null(gdd_map[[key1]]))
    stop("Missing GDD map for season 1: ", key1, ". Ensure grids include base=", gdd1_base, " cap=", kdd1_thr)
  if (is.null(gdd_map[[key2]]))
    stop("Missing GDD map for season 2: ", key2, ". Ensure grids include base=", gdd2_base, " cap=", kdd2_thr)
  
  G1 <- gdd_map[[key1]]
  H1 <- kdd_map[[paste0("h", kdd1_thr)]]
  G2 <- gdd_map[[key2]]
  H2 <- kdd_map[[paste0("h", kdd2_thr)]]
  
  W  <- weather_monthly[, c("fips","year","month","ppt")]
  W$ppt[is.na(W$ppt)] <- 0
  W$ppt <- W$ppt / 100
  
  S1 <- merge(merge(G1, H1, by = c("fips","year","month"), all = TRUE), W, by = c("fips","year","month"), all.x = TRUE)
  S2 <- merge(merge(G2, H2, by = c("fips","year","month"), all = TRUE), W, by = c("fips","year","month"), all.x = TRUE)
  
  A1 <- S1[S1$month %in% s1_months, ]
  A1 <- aggregate(. ~ fips + year, data = A1[, c("fips","year","gdd_m","kdd_m","ppt")], FUN = function(z) sum(z, na.rm = TRUE))
  names(A1) <- c("fips","year","gdd_s1","kdd_s1","ppt_s1")
  
  A2 <- S2[S2$month %in% s2_months, ]
  A2 <- aggregate(. ~ fips + year, data = A2[, c("fips","year","gdd_m","kdd_m","ppt")], FUN = function(z) sum(z, na.rm = TRUE))
  names(A2) <- c("fips","year","gdd_s2","kdd_s2","ppt_s2")
  
  out <- merge(A1, A2, by = c("fips","year"), all = TRUE)
  for (v in c("gdd_s1","kdd_s1","ppt_s1","gdd_s2","kdd_s2","ppt_s2")) {
    if (!v %in% names(out)) out[[v]] <- 0
    out[[v]][is.na(out[[v]])] <- 0
    out[[v]] <- as.numeric(out[[v]])
  }
  out
}

# -----------------------------
# CV SCORER (county FE in CV, train-centered trend; predict fixef=FALSE)
# -----------------------------
cv_score_fe <- function(df, xvars, yvar = "yield",
                        min_test = 10, min_train = 30,
                        add_quad = FALSE){
  DF  <- as.data.frame(df)
  yrs <- sort(unique(DF$year))
  r2s <- numeric(0)
  
  for (yv in yrs){
    tr <- DF[DF$year != yv, , drop = FALSE]
    te <- DF[DF$year == yv, , drop = FALSE]
    if (nrow(te) < min_test || nrow(tr) < min_train) next
    
    # Train-centered trend (no leakage)
    t0    <- mean(tr$year, na.rm = TRUE)
    tr$t  <- tr$year - t0
    te$t  <- te$year - t0
    if (add_quad) { tr$t2 <- tr$t^2; te$t2 <- te$t^2 }
    
    cols <- intersect(xvars, names(tr))
    cols <- c(cols, "t", if (add_quad) "t2")
    vrs  <- sapply(tr[, cols, drop = FALSE], function(x) stats::var(x, na.rm = TRUE))
    keep <- names(vrs)[is.finite(vrs) & vrs > 0]
    if (!length(keep)) next
    
    # County FE only during CV
    fml <- as.formula(paste0(yvar, " ~ ", paste(keep, collapse = " + "), " | fips"))
    fit <- tryCatch(fixest::feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(fit)) next
    
    # Predict without FE so held-out FE are not required
    yhat <- tryCatch(predict(fit, newdata = te, fixef = FALSE), error = function(e) rep(NA_real_, nrow(te)))
    ok   <- is.finite(yhat) & is.finite(te[[yvar]])
    if (!any(ok)) next
    
    mse <- mean((te[[yvar]][ok] - yhat[ok])^2)
    vte <- stats::var(te[[yvar]][ok])
    if (!is.finite(mse) || !is.finite(vte) || vte <= 0) next
    
    r2s <- c(r2s, 1 - mse / vte)
  }
  if (length(r2s) == 0) -Inf else mean(r2s)
}

# -----------------------------
# (A) SINGLE-CANDIDATE SANITY CHECK
# -----------------------------
pick <- function(v, i) v[pmin(i, length(v))]  # safe pick when grid has length 1

sanity <- list(
  s1_start = pick(s1_start_grid, 1),
  s2_start = pick(s2_start_grid, 1),
  s2_end   = pick(s2_end_grid,   1),
  g1 = pick(gdd_bases_grid, 1),
  h1 = pick(kdd_thr_grid,    1),
  g2 = pick(gdd_bases_grid, 2),  # if one base only, falls back to 1
  h2 = pick(kdd_thr_grid,    2)   # second threshold if available; else first
)

# preflight keys
k1 <- sprintf("b%d_c%d", sanity$g1, sanity$h1)
k2 <- sprintf("b%d_c%d", sanity$g2, sanity$h2)
message(sprintf("Sanity candidate: S1 %02d-%02d | S2 %02d-%02d | GDD maps (%s, %s) KDD(h%d,h%d)",
                sanity$s1_start, sanity$s2_start-1, sanity$s2_start, sanity$s2_end, k1, k2, sanity$h1, sanity$h2))
stopifnot(k1 %in% names(gdd_map), k2 %in% names(gdd_map),
          paste0("h", sanity$h1) %in% names(kdd_map),
          paste0("h", sanity$h2) %in% names(kdd_map))

feats_test <- build_2season_features_bins(
  s1_start = sanity$s1_start, s2_start = sanity$s2_start, s2_end = sanity$s2_end,
  gdd1_base = sanity$g1, kdd1_thr = sanity$h1,
  gdd2_base = sanity$g2, kdd2_thr = sanity$h2,
  gdd_map = gdd_map, kdd_map = kdd_map, weather_monthly = weatherdata
)

df_test <- merge(yld_raw[, c("fips","year","yield","State.ANSI")], feats_test, by = c("fips","year"), all = FALSE)
df_test$ppt2_s1 <- df_test$ppt_s1^2
df_test$ppt2_s2 <- df_test$ppt_s2^2

xvars <- c("gdd_s1","kdd_s1","ppt_s1","ppt2_s1","gdd_s2","kdd_s2","ppt_s2","ppt2_s2")

cat("\n[Sanity] rows:", nrow(df_test), " years:", length(unique(df_test$year)), "\n")
cat("[Sanity] any NA RHS? ", any(!complete.cases(df_test[, xvars, drop = FALSE])), "\n")
cat("[Sanity] var(RHS):\n"); print(round(sapply(df_test[, xvars, drop = FALSE], var, na.rm = TRUE), 6))

# Small fold printout
yrs_demo <- head(sort(unique(df_test$year)), 5)
cat("[Sanity] mini LOYO years:", paste(yrs_demo, collapse = ", "), "\n")
for (yv in yrs_demo){
  tr <- df_test[df_test$year != yv, , drop = FALSE]
  te <- df_test[df_test$year == yv, , drop = FALSE]
  t0 <- mean(tr$year); tr$t <- tr$year - t0; te$t <- te$year - t0
  cols <- intersect(xvars, names(tr)); cols <- c(cols, "t")
  vrs  <- sapply(tr[, cols, drop = FALSE], function(x) var(x, na.rm = TRUE))
  keep <- names(vrs)[is.finite(vrs) & vrs > 0]
  cat("  Year", yv, " keep:", if (length(keep)) paste(keep, collapse = ", ") else "<none>", "\n")
}

sanity_r2 <- cv_score_fe(df_test, xvars, "yield")
cat("[Sanity] LOYO R² =", sanity_r2, "\n")
if (!is.finite(sanity_r2) || sanity_r2 <= -1) stop("Sanity CV failed. Check printed diagnostics.")

# -----------------------------
# (B) FULL GRID SEARCH
# -----------------------------
gen_two_season_grid <- function(s1_starts, s2_starts, s2_ends){
  grid <- expand.grid(s1_start = s1_starts, s2_start = s2_starts, s2_end = s2_ends, KEEP.OUT.ATTRS = FALSE)
  subset(grid, s1_start < s2_start & s2_start <= s2_end)
}
parts_grid <- gen_two_season_grid(s1_start_grid, s2_start_grid, s2_end_grid)

top_tbl <- list()
best    <- list(score = -Inf); combo_counter <- 0L

for (i in seq_len(nrow(parts_grid))){
  s1_start <- parts_grid$s1_start[i]; s2_start <- parts_grid$s2_start[i]; s2_end <- parts_grid$s2_end[i]
  for (g1 in gdd_bases_grid) for (g2 in gdd_bases_grid) {
    for (h1 in kdd_thr_grid) for (h2 in kdd_thr_grid) {
      
      combo_counter <- combo_counter + 1L
      if (combo_counter %% 20 == 0) message(sprintf("... tried %d combos so far", combo_counter))
      
      feats <- build_2season_features_bins(
        s1_start, s2_start, s2_end,
        gdd1_base = g1, kdd1_thr = h1,
        gdd2_base = g2, kdd2_thr = h2,
        gdd_map = gdd_map, kdd_map = kdd_map,
        weather_monthly = weatherdata
      )
      
      dfX <- merge(yld_raw[, c("fips","year","yield","State.ANSI")], feats, by = c("fips","year"), all = FALSE)
      for (v in c("gdd_s1","kdd_s1","ppt_s1","gdd_s2","kdd_s2","ppt_s2")) {
        dfX[[v]][is.na(dfX[[v]])] <- 0
        dfX[[v]] <- as.numeric(dfX[[v]])
      }
      dfX$ppt2_s1 <- dfX$ppt_s1^2; dfX$ppt2_s2 <- dfX$ppt_s2^2
      
      sc  <- cv_score_fe(dfX, xvars, "yield")
      rec <- data.frame(score = sc, s1_start, s2_start, s2_end, gdd1_base = g1, kdd1_thr = h1, gdd2_base = g2, kdd2_thr = h2)
      top_tbl[[length(top_tbl) + 1L]] <- rec
      
      if (is.finite(sc) && sc > best$score) best <- c(as.list(rec), list(dfX = dfX, xvars = xvars))
    }
  }
}

top_tbl <- do.call(rbind, top_tbl)
top_tbl <- top_tbl[order(-top_tbl$score), ]

cat(sprintf(
  "\nBest LOYO R^2 = %.3f | S1 %02d-%02d | S2 %02d-%02d | GDD(%g,%g) KDD(%g,%g)\n",
  best$score, best$s1_start, best$s2_start - 1, best$s2_start, best$s2_end,
  best$gdd1_base, best$gdd2_base, best$kdd1_thr, best$kdd2_thr
))

# Final refit on full data (with year FE now)
fml_best <- as.formula(paste0("yield ~ ", paste(xvars, collapse = " + "), " | fips + year"))
mod_best <- feols(fml_best, data = best$dfX, cluster = ~ State.ANSI)
print(summary(mod_best))

# Top-15 plot
topK <- head(top_tbl, 15)
topK$label <- with(topK, sprintf("S1 %02d-%02d | S2 %02d-%02d | GDD(%g,%g) KDD(%g,%g)",
                                 s1_start, s2_start - 1, s2_start, s2_end,
                                 gdd1_base, gdd2_base, kdd1_thr, kdd2_thr))
gg <- ggplot(topK, aes(x = reorder(label, score), y = score)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "LOYO CV R²", title = "Top-15 two-season configurations (bin-derived GDD/KDD, cap=KDD)") +
  theme_minimal(base_size = 12)
print(gg)

