# ============================================================
# Single-season CV (Apr–Oct) — select by MEDIAN SKILL
# No centering, NO winsorization of precipitation
# Baseline: yield ~ t + i(State.ANSI, t) | fips
# Candidate RHS:
#   GDD+HDD+PPT (linear), GDD+HDD+PPT+PPT^2
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

# -----------------------------
# PATHS / SETTINGS
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")

years_use <- c(1981:2019,2021:2022)
months_in_season <- 4:10
gdd_base_grid    <- c(5, 10)
heat_thr_grid    <- c(28, 30, 32)
filt_states_north <- c(8,16,26,27,30,31,38,39,41,53,56)

# Candidate RHS (no centering; NO winsorization; PPT as-is)
rhs_menu <- list(
  GDD_HDD_P2 = "GDD_season + HDD_season + PPT_season + I(PPT_season^2)",
  GDD_HDD_P  = "GDD_season + HDD_season + PPT_season"
)

# -----------------------------
# UTILITIES
# -----------------------------
ensure_state_ansi <- function(DT){
  DT <- as.data.table(DT)
  if (!("State.ANSI" %in% names(DT)) || any(!is.finite(suppressWarnings(as.numeric(DT$State.ANSI))))) {
    DT[, State.ANSI := as.integer(fips %/% 1000L)]
  } else {
    DT[, State.ANSI := as.integer(State.ANSI)]
  }
  DT
}

disc_bins <- function(TB){
  nm <- names(TB)
  bins <- nm[grepl("^bin_-?\\d+$", nm)]
  if (!length(bins)) stop("No bin_* columns found.")
  k <- as.integer(sub("^bin_", "", bins))
  keep <- is.finite(k)
  bins <- bins[keep]; k <- k[keep]
  ord <- order(k)
  list(bins = bins[ord], k = k[ord])
}

make_monthly_dd <- function(tempbinsdata, gdd_base, heat_thr){
  db <- disc_bins(tempbinsdata)
  MM <- as.data.table(tempbinsdata[, c("fips","year","month", db$bins), with = FALSE])
  for (b in db$bins) { MM[, (b) := suppressWarnings(as.numeric(get(b)))]; MM[is.na(get(b)), (b) := 0] }
  Hdy <- as.matrix(MM[, db$bins, with = FALSE]) / 24
  storage.mode(Hdy) <- "double"
  gdd_base <- as.numeric(gdd_base)[1]; heat_thr <- as.numeric(heat_thr)[1]
  w_gdd <- pmax(pmin(db$k, heat_thr) - gdd_base, 0)
  w_hdd <- pmax(db$k - heat_thr, 0)
  if (ncol(Hdy) != length(w_gdd)) stop("Non-conformable bins/weights.")
  out <- MM[, .(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
  out[, GDD_mon := as.numeric(Hdy %*% w_gdd)]
  out[, HDD_mon := as.numeric(Hdy %*% w_hdd)]
  out[]
}

build_season_features <- function(tempbinsdata, weatherdata, gdd_base, heat_thr){
  DD <- make_monthly_dd(tempbinsdata, gdd_base = gdd_base, heat_thr = heat_thr)
  W  <- as.data.table(weatherdata[, .(fips, year, month, ppt)])
  W[is.na(ppt), ppt := 0]
  W[, ppt := as.numeric(ppt) / 100]  # keep your scaling
  M  <- merge(DD, W, by = c("fips","year","month"), all.x = TRUE)
  M  <- M[month %in% months_in_season]
  M[, .(
    GDD_season = sum(GDD_mon, na.rm = TRUE),
    HDD_season = sum(HDD_mon, na.rm = TRUE),
    PPT_season = sum(ppt,     na.rm = TRUE)
  ), by = .(fips, year)]
}

# -----------------------------
# LOYO SKILL for a given RHS (string) on a given df
# - Baseline: t + i(State.ANSI, t) | fips
# - NO winsorization, NO centering
# -----------------------------
loyo_skill_for_rhs <- function(df, rhs_string){
  DT <- as.data.table(df)
  if (!"State.ANSI" %in% names(DT)) DT[, State.ANSI := as.integer(fips %/% 1000L)]
  yrs <- sort(unique(DT$year))
  skills <- c()
  
  for (yv in yrs){
    tr <- copy(DT[year != yv]); te <- copy(DT[year == yv])
    if (nrow(tr) < 30 || nrow(te) < 10) next
    
    # center time on TRAIN
    t0 <- mean(tr$year, na.rm = TRUE)
    tr[, t := year - t0]; te[, t := year - t0]
    
    # Baseline
    base <- tryCatch(feols(yield ~ t + i(State.ANSI, t) | fips, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(base)) next
    pb <- predict(base, newdata = te, fixef = FALSE)
    
    # Candidate
    fml <- as.formula(paste0("yield ~ t + i(State.ANSI, t) + ", rhs_string, " | fips"))
    mod <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(mod)) next
    pm <- predict(mod, newdata = te, fixef = FALSE)
    
    ok <- is.finite(te$yield) & is.finite(pb) & is.finite(pm)
    if (!any(ok)) next
    mse_b <- mean((te$yield[ok] - pb[ok])^2)
    mse_m <- mean((te$yield[ok] - pm[ok])^2)
    if (is.finite(mse_b) && mse_b > 0 && is.finite(mse_m)){
      skills <- c(skills, 1 - mse_m/mse_b)
    }
  }
  list(skills = skills,
       mean_skill = if (length(skills)) mean(skills) else -Inf,
       median_skill = if (length(skills)) median(skills) else -Inf)
}

# -----------------------------
# CV SEARCH (grid × RHS), select by MEDIAN SKILL (tie-break by MEAN SKILL)
# -----------------------------
cv_search_single <- function(yld_raw){
  recs <- list(); best <- NULL; idx <- 0L
  for (b in gdd_base_grid){
    for (H in heat_thr_grid){
      feats <- build_season_features(tempbinsdata, weatherdata, gdd_base = b, heat_thr = H)
      df <- merge(yld_raw[, .(fips, year, yield, State.ANSI)], feats, by = c("fips","year"), all = FALSE)
      df <- ensure_state_ansi(df)[State.ANSI %in% filt_states_north]
      if (nrow(df) < 120) next
      
      for (nm in names(rhs_menu)){
        rhs <- rhs_menu[[nm]]
        sc  <- loyo_skill_for_rhs(df, rhs)
        idx <- idx + 1L
        recs[[idx]] <- data.table(
          gdd_base = b, heat_thr = H, rhs = nm,
          median_skill = sc$median_skill, mean_skill = sc$mean_skill,
          n_folds = length(sc$skills)
        )
        if (is.null(best) ||
            sc$median_skill > best$median_skill + 1e-12 ||
            (abs(sc$median_skill - best$median_skill) <= 1e-12 && sc$mean_skill > best$mean_skill)){
          best <- list(
            gdd_base = b, heat_thr = H, rhs = nm,
            median_skill = sc$median_skill, mean_skill = sc$mean_skill,
            df = df
          )
        }
      }
    }
  }
  cv_table <- if (length(recs)) rbindlist(recs)[order(-median_skill, -mean_skill)] else data.table()
  list(cv_table = cv_table, best = best)
}

# -----------------------------
# Final refit on selected combo (NO winsorization)
# -----------------------------
refit_best_single <- function(df_best, rhs_name){
  D <- copy(as.data.table(df_best))
  if (!"t" %in% names(D)) { t0 <- mean(D$year, na.rm = TRUE); D[, t := year - t0] }
  rhs <- rhs_menu[[rhs_name]]
  fml <- as.formula(paste0("yield ~ ", rhs, " + i(State.ANSI, t) | fips"))
  feols(fml, data = D, cluster = ~ State.ANSI)
}

# -----------------------------
# Load outcome
# -----------------------------
load_outcome <- function(csv_path, years_use){
  dat <- fread(csv_path)
  setnames(dat,
           old = c("State ANSI","County ANSI","Year","VALUE","value","SUCROSE","Sucrose"),
           new = c("State.ANSI","County.ANSI","year","Value","Value","Value","Value"),
           skip_absent = TRUE)
  dat[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
  dat <- dat[is.finite(Value) & Value > 0]
  if (!("fips" %in% names(dat))) {
    dat[, State.ANSI  := as.integer(State.ANSI)]
    dat[, County.ANSI := as.integer(County.ANSI)]
    dat[, fips := 1000L * State.ANSI + County.ANSI]
  } else dat[, fips := as.integer(fips)]
  dat[, `:=`(year = as.integer(year), yield = log(Value))]
  dat <- dat[is.finite(fips) & year %in% years_use]
  ensure_state_ansi(dat)[]
}

# -----------------------------
# LOAD WEATHER & RUN
# -----------------------------
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)

stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))

tempbinsdata[, `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
weatherdata[,  `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
setnames(weatherdata, ppt_col[1], "ppt")

run_pipeline_single <- function(yld_raw, crop_slug, save_dir = root, run_fit = TRUE){
  label <- tools::toTitleCase(crop_slug)
  cat(sprintf("\n\n########## %s (Single Season Apr–Oct) ##########\n", toupper(label)))
  res <- cv_search_single(yld_raw)
  
  cat("\n================ ", toupper(label), " — CV (grid × RHS; north only) ================\n", sep = "")
  if (nrow(res$cv_table)) print(res$cv_table, row.names = FALSE) else cat("No valid CV results.\n")
  if (is.null(res$best)) { cat("\nNo valid spec.\n"); return(invisible(NULL)) }
  
  out_csv <- file.path(save_dir, sprintf("%s_single_season_df.csv", tolower(crop_slug)))
  fwrite(res$best$df, out_csv)
  cat(sprintf("\nSaved single-season features → %s\n", out_csv))
  
  cat("\n=== Best combo ===\n")
  cat(sprintf("GDD base: %g | Heat thr: %g | RHS: %s | median SKILL: %.3f | mean SKILL: %.3f\n",
              res$best$gdd_base, res$best$heat_thr, res$best$rhs,
              res$best$median_skill, res$best$mean_skill))
  
  if (isTRUE(run_fit)){
    mod <- refit_best_single(res$best$df, res$best$rhs)
    cat("\n=== Final model (single season; selected RHS) ===\n")
    print(summary(mod))
    invisible(list(cv_table = res$cv_table, best = res$best, model = mod, df_path = out_csv))
  } else {
    invisible(list(cv_table = res$cv_table, best = res$best, df_path = out_csv))
  }
}

for (crop in c("sunflower","canola","flaxseed")) {
  yfile <- file.path(root, sprintf("%s_yields.csv", crop))
  if (!file.exists(yfile)) { cat(sprintf("\n-- %s yields not found: %s (skipping) --\n", toupper(crop), yfile)); next }
  yld_raw <- load_outcome(yfile, years_use)
  run_pipeline_single(yld_raw, crop_slug = crop, save_dir = root, run_fit = TRUE)
}







### estimate models with best mean error

# ============================================================
# Refit models chosen by MEAN skill (instead of median skill)
# ============================================================

library(data.table)
library(fixest)

root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"

# -----------------------------
# Function to build formula safely (handles single or staged)
# -----------------------------
make_selected_formula <- function(n_stages, hdd, ppt, avail) {
  n_stages <- ifelse(is.na(n_stages) || n_stages < 2, 1, n_stages)
  rhs <- c("GDD_total", "t")
  
  # HDD block
  if (hdd == "shared_all") {
    need <- intersect(c("HDD_s1","HDD_s2", if (n_stages == 3) "HDD_s3"), avail)
    if (length(need)) rhs <- c(rhs, sprintf("I(%s)", paste(need, collapse = " + ")))
  } else {
    if (n_stages == 1) {
      need <- intersect("HDD_total", avail)
    } else {
      need <- intersect(c("HDD_s1","HDD_s2", if (n_stages == 3) "HDD_s3"), avail)
    }
    rhs <- c(rhs, need)
  }
  
  # Precipitation block (linear + quadratic)
  if (ppt == "shared_all") {
    need <- intersect(c("P_s1","P_s2", if (n_stages == 3) "P_s3"), avail)
    if (length(need)) {
      rhs <- c(rhs, sprintf("I(%s)", paste(need, collapse = " + ")))
      rhs <- c(rhs, sprintf("I((%s)^2)", paste(need, collapse = " + ")))
    }
  } else {
    if (n_stages == 1) {
      need <- intersect("P_total", avail)
      rhs <- c(rhs, need, sprintf("I(%s^2)", need))
    } else {
      for (s in c("s1","s2", if (n_stages == 3) "s3")) {
        p <- paste0("P_", s)
        if (p %in% avail) rhs <- c(rhs, p, sprintf("I(%s^2)", p))
      }
    }
  }
  
  as.formula(paste("yield ~ i(State.ANSI, year) +", paste(rhs, collapse = " + "), "| fips"))
}

# -----------------------------
# Function to refit model for given crop
# -----------------------------
refit_best_mean <- function(crop_slug) {
  cv_path <- file.path(root, sprintf("%s_cv_comparison_table.csv", tolower(crop_slug)))
  df_path <- file.path(root, sprintf("%s_best_stages_df.csv", tolower(crop_slug)))
  
  if (!file.exists(cv_path)) stop("CV results not found for ", crop_slug)
  if (!file.exists(df_path)) stop("Best-stages dataframe not found for ", crop_slug)
  
  cv <- fread(cv_path)
  best <- cv[which.max(mean_skill)]  # select by mean_skill
  
  cat("\n=== Selected by MEAN SKILL ===\n")
  print(best)
  
  df_best <- fread(df_path)
  if (!"State.ANSI" %in% names(df_best)) df_best[, State.ANSI := as.integer(fips %/% 1000L)]
  if (!"t" %in% names(df_best)) {
    t0 <- mean(df_best$year, na.rm = TRUE)
    df_best[, t := year - t0]
  }
  
  fml <- make_selected_formula(best$n_stages, best$hdd, best$ppt, names(df_best))
  mod <- feols(fml, data = df_best, cluster = ~ year)
  
  cat("\n=== Final model (selected by mean skill) ===\n")
  print(summary(mod))
  invisible(mod)
}

# -----------------------------
# Run for all three crops
# -----------------------------
for (crop in c("sunflower", "canola", "flaxseed")) {
  cat(sprintf("\n\n########## %s ##########\n", toupper(crop)))
  refit_best_mean(crop)
}

