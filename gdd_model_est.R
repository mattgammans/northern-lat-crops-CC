# ============================================================
# One script → CV + FE model for:
#   (1) YIELD  (log(Value) from sugarbeet_yields.csv)
#   (2) SUCROSE (log(Value) from sugarbeet_sucrose.csv)
# Everything else identical (features, grids, CV, clustering)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

# -----------------------------
# PATHS (EDIT THESE)
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")
yield_csv <- file.path(root, "sugarbeet_yields.csv")
sucro_csv <- file.path(root, "sugarbeet_sucrose.csv")

# -----------------------------
# CV KNOBS (same for both outcomes)
# -----------------------------
plant_start    <- 4                 # April (planting window)
plant_end      <- 5                 # May
end_grid       <- c(9, 10)          # Growing end: Sep or Oct
heat_grid      <- c(28, 30)         # HDD threshold (and GDD cap)
gdd_base_grid  <- c(5, 10)          # GDD base to search
years_use      <- 1981:2022
filt_states    <- NULL              # e.g., c(27,38,26); NULL = all states

# -----------------------------
# LOAD WEATHER (once)
# -----------------------------
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)

stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))

tempbinsdata[, `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
weatherdata[,  `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]

# precip column → ppt
ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
setnames(weatherdata, ppt_col[1], "ppt")

# -----------------------------
# HELPERS
# -----------------------------
disc_bins <- function(TB){
  nm <- names(TB)
  bins <- nm[grepl("^bin_-?\\d+$", nm)]
  if (!length(bins)) stop("No bin_* columns.")
  k <- as.integer(sub("^bin_", "", bins))
  ord <- order(k)
  list(bins = bins[ord], k = k[ord])
}

make_monthly_dd <- function(tempbinsdata, gdd_base, heat_thr){
  db  <- disc_bins(tempbinsdata)
  # Use data.table selection but return plain data.frame for base ops
  MM  <- as.data.frame(tempbinsdata[, c("fips","year","month", db$bins), with = FALSE])
  Hdy <- as.matrix(MM[, db$bins, drop = FALSE]) / 24  # hours→days
  
  w_gdd <- pmax(pmin(db$k, heat_thr) - gdd_base, 0)
  w_hdd <- pmax(db$k - heat_thr, 0)
  
  out <- MM[, c("fips","year","month"), drop = FALSE]
  out$gdd_m   <- as.numeric(Hdy %*% w_gdd)
  out$hdd_m   <- as.numeric(Hdy %*% w_hdd)
  out$frost_m <- if (any(db$k < 0)) rowSums(Hdy[, db$k < 0, drop = FALSE], na.rm = TRUE) else 0
  out
}

# LOYO (county FE in training; predict without FE)
loyo_cv_fe <- function(df, xvars, yvar = "yield", min_test = 10, min_train = 30){
  yrs <- sort(unique(df$year))
  r2s <- numeric(0)
  for (yv in yrs){
    tr <- df[df$year != yv, , drop = FALSE]
    te <- df[df$year == yv, , drop = FALSE]
    if (nrow(tr) < min_train || nrow(te) < min_test) next
    t0   <- mean(tr$year, na.rm = TRUE)
    tr$t <- tr$year - t0
    te$t <- te$year - t0
    cols <- c(intersect(xvars, names(tr)), "t")
    vrs  <- sapply(tr[, cols, drop = FALSE], function(x) var(x, na.rm = TRUE))
    keep <- names(vrs)[is.finite(vrs) & vrs > 0]
    if (!length(keep)) next
    fml  <- as.formula(paste0(yvar, " ~ ", paste(keep, collapse = " + "), " | fips"))
    fit  <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(fit)) next
    yhat <- tryCatch(predict(fit, newdata = te, fixef = FALSE), error = function(e) rep(NA_real_, nrow(te)))
    ok   <- is.finite(yhat) & is.finite(te[[yvar]])
    if (!any(ok)) next
    mse  <- mean((te[[yvar]][ok] - yhat[ok])^2)
    vte  <- var(te[[yvar]][ok])
    if (!is.finite(mse) || !is.finite(vte) || vte <= 0) next
    r2s  <- c(r2s, 1 - mse / vte)
  }
  if (length(r2s) == 0) -Inf else mean(r2s)
}

# Feature builder (Apr–May planting; early PPT & frost; growing has GDD/HDD + late PPT + late^2)
build_two_window <- function(end_month, heat_thr, gdd_base,
                             tempbinsdata, weatherdata,
                             plant_start = 4, plant_end = 5){
  DD <- make_monthly_dd(tempbinsdata, gdd_base, heat_thr)
  
  W  <- as.data.frame(weatherdata[, c("fips","year","month","ppt")])
  W$ppt[is.na(W$ppt)] <- 0; W$ppt <- as.numeric(W$ppt)/100
  
  M <- merge(DD, W, by = c("fips","year","month"), all.x = TRUE)
  M <- as.data.frame(M)
  for (v in c("gdd_m","hdd_m","frost_m","ppt")) { M[[v]][is.na(M[[v]])] <- 0 }
  
  # Planting (Apr–May) → ppt_early, frost_early
  P <- M[M$month %in% (plant_start:plant_end), , drop = FALSE]
  early <- aggregate(P[, c("ppt","frost_m")],
                     by = list(fips = P$fips, year = P$year),
                     FUN = function(z) sum(z, na.rm = TRUE))
  names(early) <- c("fips","year","ppt_early","frost_early")
  
  # Growing (June..end) → gdd_grow, hdd_grow, ppt_late (+²)
  grow_months <- (plant_end + 1):end_month
  G <- M[M$month %in% grow_months, , drop = FALSE]
  grow <- aggregate(G[, c("gdd_m","hdd_m","ppt")],
                    by = list(fips = G$fips, year = G$year),
                    FUN = function(z) sum(z, na.rm = TRUE))
  names(grow) <- c("fips","year","gdd_grow","hdd_grow","ppt_late")
  grow$ppt_late2 <- grow$ppt_late^2
  
  X <- merge(early, grow, by = c("fips","year"), all = TRUE)
  for (v in c("ppt_early","frost_early","gdd_grow","hdd_grow","ppt_late","ppt_late2")){
    if (!v %in% names(X)) X[[v]] <- 0
    X[[v]][is.na(X[[v]])] <- 0
    X[[v]] <- as.numeric(X[[v]])
  }
  X
}

# Load an outcome (yield OR sucrose) → yld_raw with yield = log(Value)
load_outcome <- function(csv_path, years_use){
  dat <- fread(csv_path)
  setnames(dat,
           old = c("State ANSI","County ANSI","Year","VALUE","value","SUCROSE","Sucrose"),
           new = c("State.ANSI","County.ANSI","year","Value","Value","Value","Value"),
           skip_absent = TRUE)
  dat[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
  dat <- dat[is.finite(Value) & Value > 0]
  dat[, `:=`(
    year = as.integer(year),
    fips = 1000L * as.integer(State.ANSI) + as.integer(County.ANSI),
    yield = log(Value)
  )]
  dat <- dat[is.finite(fips) & year %in% years_use]
  as.data.frame(dat)
}

# Run the full CV + refit for a given outcome table
run_pipeline <- function(yld_raw, label = "Outcome"){
  recs <- list(); best <- list(score = -Inf); idx <- 0L
  xvars <- c("gdd_grow","hdd_grow","ppt_early","frost_early","ppt_late","ppt_late2")
  
  for (end_m in end_grid) {
    for (thr in heat_grid) {
      for (base in gdd_base_grid) {
        feats <- build_two_window(end_month = end_m, heat_thr = thr, gdd_base = base,
                                  tempbinsdata = tempbinsdata, weatherdata = weatherdata,
                                  plant_start = plant_start, plant_end = plant_end)
        
        df <- merge(yld_raw[, c("fips","year","yield","State.ANSI")], feats,
                    by = c("fips","year"), all = FALSE)
        
        if (!is.null(filt_states)) df <- df[df$State.ANSI %in% filt_states, , drop = FALSE]
        if (nrow(df) < 120) next
        
        sc  <- loyo_cv_fe(df, xvars, "yield")
        idx <- idx + 1L
        recs[[idx]] <- data.frame(end_month = end_m, heat_thr = thr, gdd_base = base, score = sc)
        
        if (is.finite(sc) && sc > best$score)
          best <- list(end = end_m, thr = thr, base = base, score = sc, df = df)
      }
    }
  }
  
  cv_table <- if (length(recs)) {
    out <- do.call(rbind, recs)
    out[order(-out$score), ]
  } else data.frame()
  
  cat("\n================ ", toupper(label),
      " — CV (end × heat_thr × gdd_base) ================\n", sep = "")
  print(cv_table, row.names = FALSE)
  
  if (is.finite(best$score)) {
    fml <- as.formula(paste0("yield ~ ", paste(xvars, collapse = " + "), " | fips + year"))
    mod <- feols(fml, data = best$df, cluster = ~ fips + year)
    
    cat("\n=== Best spec (", label, ") ===\n", sep = "")
    cat(sprintf("Planting: %02d–%02d | Growing: %02d–%02d | GDD base: %g | Heat thr: %g\n",
                plant_start, plant_end, plant_end + 1, best$end, best$base, best$thr))
    cat(sprintf("LOYO CV R^2: %.3f\n\n", best$score))
    print(summary(mod))
    invisible(list(cv_table = cv_table, best = best, model = mod))
  } else {
    cat("\nNo valid spec (", label, ").\n", sep = "")
    invisible(NULL)
  }
}

# =========================
# RUN BOTH OUTCOMES
# =========================
yld_raw <- load_outcome(yield_csv, years_use)
yld_raw <- yld_raw[yld_raw$State.ANSI%in% c(27,38,26),]
res_yield <- run_pipeline(yld_raw, label = "Yield")

yld_raw <- load_outcome(sucro_csv, years_use)
res_sucrose <- run_pipeline(yld_raw, label = "Sucrose")

