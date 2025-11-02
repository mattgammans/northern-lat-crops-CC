# ============================================================
# One script → CV + FE model for:
#   (1) YIELD  (log(Value) from sugarbeet_yields.csv)
#   (2) SUCROSE (log(Value) from sugarbeet_sucrose.csv)
#
# Weather spec:
#  - Planting (Apr–May): ppt_early, frost_early
#  - Growing (Jun–end): temperature BINS (as SHARES of grow-days)
#       keep:   <0 (freeze), 0–5, 5–10, 20–25, 25–30,  >30  (pooled)
#       omit:   10–20 (10–15 & 15–20 are excluded entirely)
#    + ppt_late, ppt_late^2, and grow_days (season-length control)
#
# CV over growing end month: {Sep (9), Oct (10)}
# FE model clustered by fips + year
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
plant_start <- 4          # April
plant_end   <- 5          # May
end_grid    <- c(9, 10)   # growing season end: Sep or Oct
years_use   <- 1981:2022
filt_states <- NULL       # e.g., c(27,38,26); NULL = all states

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

# Build *monthly* day counts per bin (hours → days)
make_monthly_bin_days <- function(tempbinsdata){
  db <- disc_bins(tempbinsdata)
  MM <- as.data.frame(tempbinsdata[, c("fips","year","month", db$bins), with = FALSE])
  H  <- as.matrix(MM[, db$bins, drop = FALSE]) / 24  # hours → days
  out <- MM[, c("fips","year","month"), drop = FALSE]
  attr(out, "H") <- H
  attr(out, "k") <- db$k
  out
}

# Aggregate monthly bin days into groups over a set of months
# Groups kept: freeze(<0), 0-5, 5-10, 20-25, 25-30, >=30; 10-15 and 15-20 are omitted entirely
sum_bins_over_window <- function(monthly_df, months_vec){
  H <- attr(monthly_df, "H"); k <- attr(monthly_df, "k")
  D <- monthly_df[monthly_df$month %in% months_vec, c("fips","year"), drop = FALSE]
  if (nrow(D) == 0) {
    return(data.frame(fips = integer(0), year = integer(0),
                      bin_freeze = numeric(0), bin0_5 = numeric(0), bin5_10 = numeric(0),
                      bin20_25 = numeric(0), bin25_30 = numeric(0), bin_gt30 = numeric(0)))
  }
  idx <- which(monthly_df$month %in% months_vec)
  Hsub <- H[idx, , drop = FALSE]
  
  sum_range <- function(L, U) {
    cols <- which(k >= L & k <= U)
    if (length(cols)) rowSums(Hsub[, cols, drop = FALSE], na.rm = TRUE) else 0
  }
  D$bin_freeze <- if (any(k < 0)) rowSums(Hsub[, k < 0, drop = FALSE], na.rm = TRUE) else 0
  D$bin0_5     <- sum_range(0, 4)
  D$bin5_10    <- sum_range(5, 9)
  # OMIT 10-15 and 15-20 entirely (no columns created)
  D$bin20_25   <- sum_range(20, 24)
  D$bin25_30   <- sum_range(25, 29)
  # POOLED tail >30
  tail_cols <- which(k >= 30)
  D$bin_gt30  <- if (length(tail_cols)) rowSums(Hsub[, tail_cols, drop = FALSE], na.rm = TRUE) else 0
  
  agg <- aggregate(D[, -c(1:2), drop = FALSE], by = list(fips = D$fips, year = D$year),
                   FUN = function(z) sum(z, na.rm = TRUE))
  return(agg)
}

# LOYO (county FE in training; predict without FE)
loyo_cv_fe <- function(df, xvars, yvar = "yield", min_test = 10, min_train = 30){
  yrs <- sort(unique(df$year)); r2s <- numeric(0)
  for (yv in yrs){
    tr <- df[df$year != yv, , drop = FALSE]
    te <- df[df$year == yv, , drop = FALSE]
    if (nrow(tr) < min_train || nrow(te) < min_test) next
    t0   <- mean(tr$year, na.rm = TRUE)
    tr$t <- tr$year - t0; te$t <- te$year - t0
    cols <- c(intersect(xvars, names(tr)), "t")
    vrs  <- sapply(tr[, cols, drop = FALSE], function(x) var(x, na.rm = TRUE))
    keep <- names(vrs)[is.finite(vrs) & vrs > 0]; if (!length(keep)) next
    fml  <- as.formula(paste0(yvar," ~ ", paste(keep, collapse=" + "), " | fips"))
    fit  <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(fit)) next
    yhat <- tryCatch(predict(fit, newdata = te, fixef = FALSE), error = function(e) rep(NA_real_, nrow(te)))
    ok   <- is.finite(yhat) & is.finite(te[[yvar]]); if (!any(ok)) next
    mse  <- mean((te[[yvar]][ok] - yhat[ok])^2); vte <- var(te[[yvar]][ok])
    if (!is.finite(mse) || !is.finite(vte) || vte <= 0) next
    r2s  <- c(r2s, 1 - mse / vte)
  }
  if (length(r2s) == 0) -Inf else mean(r2s)
}

# -----------------------------
# Feature builder (bins as SHARES)
# Planting: ppt_early, frost_early
# Growing : shares of day-bins over June..end + ppt_late + ppt_late^2 + grow_days
# Kept bins: freeze, 0–5, 5–10, 20–25, 25–30, >30  (10–20 omitted)
# -----------------------------
build_bins_two_window <- function(end_month,
                                  tempbinsdata, weatherdata,
                                  plant_start = 4, plant_end = 5){
  # monthly bin-day table
  MB <- make_monthly_bin_days(tempbinsdata)
  
  # monthly ppt (/100)
  W  <- as.data.frame(weatherdata[, c("fips","year","month","ppt")])
  W$ppt[is.na(W$ppt)] <- 0; W$ppt <- as.numeric(W$ppt)/100
  
  # Merge for planting ppt aggregation
  M <- merge(as.data.frame(MB), W, by = c("fips","year","month"), all.x = TRUE)
  for (v in c("ppt")) M[[v]][is.na(M[[v]])] <- 0
  
  # Planting window (Apr–May) → ppt_early, frost_early from freeze bin over planting months
  P_months <- plant_start:plant_end
  P_bins   <- sum_bins_over_window(MB, P_months)  # has freeze bin
  P <- merge(P_bins, aggregate(M[M$month %in% P_months, "ppt", drop = FALSE],
                               by = list(fips = M$fips[M$month %in% P_months],
                                         year = M$year[M$month %in% P_months]),
                               FUN = function(z) sum(z, na.rm = TRUE)),
             by = c("fips","year"), all = TRUE)
  names(P)[names(P) == "ppt"] <- "ppt_early"
  # frost_early = freeze-day count during planting
  P$frost_early <- if ("bin_freeze" %in% names(P)) P$bin_freeze else 0
  early <- P[, c("fips","year","ppt_early","frost_early")]
  for (v in names(early)[-c(1,2)]) { early[[v]][is.na(early[[v]])] <- 0 }
  
  # Growing window (June..end): bin shares + late PPT (+sq)
  grow_months <- (plant_end + 1):end_month
  G_bins <- sum_bins_over_window(MB, grow_months)  # counts kept bins
  
  # Total grow days (for shares and optional control)
  G_bins$grow_days <- with(G_bins,
                           bin_freeze + bin0_5 + bin5_10 + bin20_25 + bin25_30 + bin_gt30
  )
  
  # Convert to shares; keep grow_days as a control
  to_share <- c("bin_freeze","bin0_5","bin5_10","bin20_25","bin25_30","bin_gt30")
  for (v in to_share) {
    G_bins[[paste0(v,"_sh")]] <- ifelse(G_bins$grow_days > 0, G_bins[[v]] / G_bins$grow_days, 0)
  }
  
  G_bins_sh <- G_bins[, c("fips","year","grow_days",
                          paste0(to_share, "_sh"))]
  
  # Late PPT over growing months (+ squared)
  G_ppt <- aggregate(M[M$month %in% grow_months, "ppt", drop = FALSE],
                     by = list(fips = M$fips[M$month %in% grow_months],
                               year = M$year[M$month %in% grow_months]),
                     FUN = function(z) sum(z, na.rm = TRUE))
  names(G_ppt) <- c("fips","year","ppt_late")
  G_ppt$ppt_late2 <- G_ppt$ppt_late^2
  
  # Assemble
  X <- Reduce(function(a,b) merge(a,b, by = c("fips","year"), all = TRUE),
              list(early, G_bins_sh, G_ppt))
  
  for (v in setdiff(names(X), c("fips","year"))) {
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

# Run the full CV + refit for a given outcome table (BIN TEMP version, with shares)
run_pipeline_bins <- function(yld_raw, label = "Outcome"){
  recs <- list(); best <- list(score = -Inf); idx <- 0L
  
  # bin shares used (10–20 omitted; >30 pooled); include grow_days control
  bin_vars <- c("bin_freeze_sh","bin0_5_sh","bin5_10_sh","bin20_25_sh","bin25_30_sh","bin_gt30_sh")
  xvars <- c(bin_vars, "ppt_early","frost_early","ppt_late","ppt_late2", "grow_days")
  
  for (end_m in end_grid) {
    feats <- build_bins_two_window(end_month = end_m,
                                   tempbinsdata = tempbinsdata,
                                   weatherdata = weatherdata,
                                   plant_start = plant_start, plant_end = plant_end)
    
    df <- merge(yld_raw[, c("fips","year","yield","State.ANSI")], feats,
                by = c("fips","year"), all = FALSE)
    if (!is.null(filt_states)) df <- df[df$State.ANSI %in% filt_states, , drop = FALSE]
    if (nrow(df) < 120) next
    
    sc  <- loyo_cv_fe(df, xvars, "yield")
    idx <- idx + 1L
    recs[[idx]] <- data.frame(end_month = end_m, score = sc)
    
    if (is.finite(sc) && sc > best$score)
      best <- list(end = end_m, score = sc, df = df)
  }
  
  cv_table <- if (length(recs)) {
    out <- do.call(rbind, recs)
    out[order(-out$score), ]
  } else data.frame()
  
  cat("\n================ ", toupper(label),
      " — CV (end month; BIN TEMP with shares; omit 10–20; pool >30) ================\n", sep = "")
  print(cv_table, row.names = FALSE)
  
  if (is.finite(best$score)) {
    fml <- as.formula(paste0("yield ~ ",
                             paste(xvars, collapse = " + "),
                             " | fips + year"))
    mod <- feols(fml, data = best$df, cluster = ~ fips + year)
    
    cat("\n=== Best spec (", label, "; BIN TEMP SHARES) ===\n", sep = "")
    cat(sprintf("Planting: %02d–%02d | Growing: %02d–%02d (bins)\n",
                plant_start, plant_end, plant_end + 1, best$end))
    cat("Bins kept (shares): freeze, 0–5, 5–10, 20–25, 25–30, >30; omitted: 10–15 & 15–20\n")
    cat(sprintf("LOYO CV R^2: %.3f\n\n", best$score))
    print(summary(mod))
    invisible(list(cv_table = cv_table, best = best, model = mod))
  } else {
    cat("\nNo valid spec (", label, ").\n", sep = "")
    invisible(NULL)
  }
}

# =========================
# RUN BOTH OUTCOMES (BIN TEMP SHARES)
# =========================
yld_raw <- load_outcome(yield_csv, years_use)
res_yield_bins <- run_pipeline_bins(yld_raw, label = "Yield")

yld_raw <- load_outcome(sucro_csv, years_use)
res_sucrose_bins <- run_pipeline_bins(yld_raw, label = "Sucrose")
