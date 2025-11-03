# ============================================================
# Multi-season CV (S1, S2A only) with constrained specs
# - Temp: GDD/HDD (bases {5,10}, caps {26..32}) OR 5°C bins (start=10)
# - S2A forbids mixing (both seasons GDD or both bins)
# - PPT: ALWAYS quadratic for all seasons in the split
# - Metrics: median skill, trimmed mean (drop top2 & bottom2), mean
# Baseline: yield ~ t + i(State.ANSI, t) | fips
# ============================================================

suppressPackageStartupMessages({ library(data.table); library(fixest) })

# -----------------------------
# PATHS / SETTINGS
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")
years_use <- c(1981:2022)

# Northern-state filter you used earlier (adjust if you want SD=46 included)
filt_states_north <- c(8,16,26,27,30,31,38,39,41,53,56)

# GDD/HDD grids under the new constraints
gdd_base_grid <- c(5, 10)
hdd_cap_grid  <- 26:32

# -----------------------------
# LOAD WEATHER
# -----------------------------
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)
stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))
tempbinsdata[, `:=`(fips=as.integer(fips), year=as.integer(year), month=as.integer(month))]
weatherdata[,  `:=`(fips=as.integer(fips), year=as.integer(year), month=as.integer(month))]

ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
setnames(weatherdata, ppt_col[1], "ppt")

# ============================================================
# Helpers
# ============================================================
ensure_state_ansi <- function(DT){
  DT <- as.data.table(DT)
  if (!"State.ANSI" %in% names(DT)) DT[, State.ANSI := as.integer(fips %/% 1000L)]
  DT[, State.ANSI := as.integer(State.ANSI)]
  DT
}

disc_bins <- function(TB){
  nm   <- names(TB)
  bins <- nm[grepl("^bin_-?\\d+$", nm)]
  if (!length(bins)) stop("No bin_* columns found in temp bins file.")
  k    <- as.integer(sub("^bin_", "", bins))
  keep <- is.finite(k)
  bins <- bins[keep]; k <- k[keep]
  ord  <- order(k)
  list(bins = bins[ord], k = k[ord])
}

# Monthly GDD/HDD (hours -> days)
make_monthly_gddhdd <- function(tempbins, base, cap){
  stopifnot(all(c("fips","year","month") %in% names(tempbins)))
  DB <- disc_bins(tempbins)
  TB <- as.data.table(tempbins)[, c("fips","year","month", DB$bins), with = FALSE]
  for (b in DB$bins) {
    TB[, (b) := suppressWarnings(as.numeric(get(b)))]
    TB[is.na(get(b)), (b) := 0]
    TB[, (b) := get(b) / 24]  # hours -> days
  }
  Hdy <- as.matrix(TB[, DB$bins, with = FALSE])
  storage.mode(Hdy) <- "double"
  base <- as.numeric(base); cap <- as.numeric(cap)
  w_gdd <- pmax(pmin(DB$k, cap) - base, 0)
  w_hdd <- pmax(DB$k - cap, 0)
  if (ncol(Hdy) != length(w_gdd)) stop("Non-conformable bins/weights.")
  out <- TB[, .(fips = as.integer(fips),
                year = as.integer(year),
                month = as.integer(month))]
  out[, GDD_mon := as.numeric(Hdy %*% w_gdd)]
  out[, HDD_mon := as.numeric(Hdy %*% w_hdd)]
  out[]
}

# Collapse to 5°C groups for selected months; county-year wide
collapse_bins_5C_by_months <- function(tempbins, months_keep, suffix){
  stopifnot(all(c("fips","year","month") %in% names(tempbins)))
  DB <- disc_bins(tempbins)
  TB <- as.data.table(tempbins)[month %in% months_keep, c("fips","year","month", DB$bins), with = FALSE]
  for (b in DB$bins) {
    TB[, (b) := suppressWarnings(as.numeric(get(b)))]
    TB[is.na(get(b)), (b) := 0]
    TB[, (b) := get(b) / 24]  # hours -> days
  }
  g5 <- floor(DB$k / 5) * 5
  groups <- split(DB$bins, g5)
  summed <- TB[, {
    sums <- lapply(groups, function(cols) rowSums(.SD[, ..cols], na.rm = TRUE))
    as.data.table(sums)
  }, by = .(fips, year)]
  wide <- summed[, lapply(.SD, sum, na.rm = TRUE), by = .(fips, year)]
  gnames <- names(groups)
  setnames(wide, old = gnames, new = paste0("bin5C_", gnames, "_", suffix))
  wide[]
}

# ============================================================
# CACHES (shared across crops)
# ============================================================
.dd_cache       <- new.env(parent = emptyenv())  # monthly GDD/HDD per (base, cap)
.bin5_cache     <- new.env(parent = emptyenv())  # 5°C bins per season window
.ppt_cache      <- new.env(parent = emptyenv())  # PPT per season window
.baseline_cache <- new.env(parent = emptyenv())  # LOYO baseline preds per crop/fold

key_dd   <- function(base, cap)   paste0("dd::", base, "_", cap)
key_bins <- function(season)      paste0("bin5::", season)
key_ppt  <- function(season)      paste0("ppt::", season)
fold_key <- function(crop, yrs)   paste0(crop, "::", paste(sort(as.integer(yrs)), collapse=","))

# Precompute monthly DD/HDD for all (base, cap)
precompute_dd <- function(tempbins){
  for (base in gdd_base_grid) {
    for (cap in hdd_cap_grid) {
      k <- key_dd(base, cap)
      if (!exists(k, .dd_cache, inherits = FALSE)) {
        assign(k, make_monthly_gddhdd(tempbins, base, cap), .dd_cache)
      }
    }
  }
}

# Precompute 5°C bins per season window
precompute_bin5 <- function(tempbins){
  windows <- list(SPRING = 4:6, SUMMER = 7:10, ALL = 4:10)
  for (nm in names(windows)){
    k <- key_bins(nm)
    if (!exists(k, .bin5_cache, inherits = FALSE)) {
      assign(k, collapse_bins_5C_by_months(tempbins, windows[[nm]], suffix = tolower(nm)), .bin5_cache)
    }
  }
}

# Precompute PPT per season window (always quadratic later)
precompute_ppt <- function(weather){
  stopifnot(all(c("fips","year","month","ppt") %in% names(weather)))
  W <- as.data.table(weather)
  W[, ppt := as.numeric(ppt)]; W[is.na(ppt), ppt := 0]
  windows <- list(SPRING = 4:6, SUMMER = 7:10, ALL = 4:10)
  for (nm in names(windows)){
    k <- key_ppt(nm)
    if (!exists(k, .ppt_cache, inherits = FALSE)) {
      tmp <- W[month %in% windows[[nm]], .(PPT = sum(ppt, na.rm = TRUE)), by = .(fips, year)]
      setnames(tmp, "PPT", paste0("PPT_", nm))
      assign(k, tmp[], .ppt_cache)
    }
  }
}

# -----------------------------
# Parse a "gdd{base}_{cap}" string
# -----------------------------
parse_gdd_spec <- function(spec){
  m <- regexec("^gdd(\\d+)_(\\d+)$", spec)
  mm <- regmatches(spec, m)[[1]]
  if (length(mm) != 3) stop("Bad GDD spec: ", spec)
  list(base = as.integer(mm[2]), cap = as.integer(mm[3]))
}

# -----------------------------
# Build features for a candidate
# split: "S2A" (SPRING/SUMMER) or "S1" (ALL)
# temp_assign: named list per season: either "bin5C_10" OR "gdd{base}_{cap}"
# Always append PPT_<season> and squared
# -----------------------------
build_candidate_features <- function(split, temp_assign){
  if (split == "S2A") seasons <- c("SPRING","SUMMER")
  if (split == "S1")  seasons <- c("ALL")
  feats <- list()
  
  # Temperature blocks
  for (sn in seasons){
    spec <- temp_assign[[sn]]
    if (is.null(spec)) next
    
    if (spec == "bin5C_10"){
      feats[[length(feats) + 1L]] <- get(key_bins(sn), .bin5_cache)
    } else if (grepl("^gdd\\d+_\\d+$", spec)) {
      pp <- parse_gdd_spec(spec)
      DD <- get(key_dd(pp$base, pp$cap), .dd_cache)
      months <- if (sn == "SPRING") 4:6 else if (sn == "SUMMER") 7:10 else 4:10
      agg <- DD[month %in% months, .(GDD = sum(GDD_mon), HDD = sum(HDD_mon)), by = .(fips, year)]
      setnames(agg, c("GDD","HDD"), c(paste0("GDD_", sn), paste0("HDD_", sn)))
      feats[[length(feats) + 1L]] <- agg
    } else {
      stop("Unknown temp spec: ", spec)
    }
  }
  
  # Always add PPT (linear + quadratic) for all seasons in the split
  for (sn in seasons){
    P <- get(key_ppt(sn), .ppt_cache)
    set(P, j = paste0("PPT_", sn, "_sq"), value = P[[paste0("PPT_", sn)]]^2)
    feats[[length(feats) + 1L]] <- P
  }
  
  if (!length(feats)) return(data.table(fips = integer(), year = integer()))
  out <- Reduce(function(a, b) merge(a, b, by = c("fips","year"), all = TRUE), feats)
  out[is.na(out)] <- 0
  out[]
}

# -----------------------------
# RHS builder from available columns
# -----------------------------
build_rhs_terms <- function(df_names, split, temp_assign){
  rhs <- character(0)
  seasons <- if (split == "S2A") c("SPRING","SUMMER") else "ALL"
  
  # Temperature terms
  for (sn in seasons){
    spec <- temp_assign[[sn]]
    if (spec == "bin5C_10"){
      suf <- tolower(sn)
      bn <- grep(paste0("^bin5C_(-?\\d+)_", suf, "$"), df_names, value = TRUE)
      if (length(bn)) {
        centers <- as.integer(sub(paste0("^bin5C_(-?\\d+)_", suf, "$"), "\\1", bn))
        bn <- bn[centers >= 10]              # start at 10°C by design
        if (length(bn)) rhs <- c(rhs, sprintf("`%s`", bn))
      }
    } else if (grepl("^gdd\\d+_\\d+$", spec)) {
      for (nm in c(paste0("GDD_", sn), paste0("HDD_", sn))){
        if (nm %in% df_names) rhs <- c(rhs, nm)
      }
    }
  }
  
  # Always add PPT blocks
  for (sn in seasons){
    v <- paste0("PPT_", sn)
    if (v %in% df_names) rhs <- c(rhs, v, sprintf("I(%s^2)", v))
  }
  
  if (!length(rhs)) "0" else paste(rhs, collapse = " + ")
}

# -----------------------------
# LOYO SKILL (adds trimmed mean)
# -----------------------------
get_baseline_preds <- function(crop_slug, tr, te){
  k <- fold_key(crop_slug, unique(te$year))
  if (exists(k, .baseline_cache, inherits = FALSE)) return(get(k, .baseline_cache))
  base <- feols(yield ~ t + i(State.ANSI, t) | fips, data = tr, cluster = ~ State.ANSI)
  pb   <- predict(base, newdata = te, fixef = FALSE)
  assign(k, pb, .baseline_cache)
  pb
}

loyo_skill_df <- function(crop_slug, df, rhs_string, min_train = 30, min_test = 10){
  DT <- as.data.table(df)
  yrs <- sort(unique(DT$year))
  skills <- numeric(0)
  
  for (yv in yrs){
    tr <- copy(DT[year != yv])
    te <- copy(DT[year == yv])
    if (nrow(tr) < min_train || nrow(te) < min_test) next
    
    # center t on TRAIN
    t0 <- mean(tr$year, na.rm = TRUE)
    tr[, t := year - t0]
    te[, t := year - t0]
    
    # baseline (cached)
    pb <- get_baseline_preds(crop_slug, tr, te)
    
    # candidate
    fml <- as.formula(paste0("yield ~ t + i(State.ANSI, t) + ", rhs_string, " | fips"))
    fit <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(fit)) next
    
    pm <- tryCatch(predict(fit, newdata = te, fixef = FALSE),
                   error = function(e) rep(NA_real_, nrow(te)))
    
    ok   <- is.finite(te$yield) & is.finite(pb) & is.finite(pm)
    if (!any(ok)) next
    mseB <- mean((te$yield[ok] - pb[ok])^2)
    mseM <- mean((te$yield[ok] - pm[ok])^2)
    if (is.finite(mseB) && mseB > 0 && is.finite(mseM)){
      skills <- c(skills, 1 - mseM/mseB)
    }
  }
  
  # trimmed mean: drop 2 highest and 2 lowest folds if possible
  mean_trim2 <- if (length(skills) > 4) {
    xs <- sort(skills)
    mean(xs[(2+1):(length(xs)-2)])
  } else {
    mean(skills)
  }
  
  list(skills = skills,
       median_skill = if (length(skills)) median(skills) else -Inf,
       mean_trim2   = if (length(skills)) mean_trim2 else -Inf,
       mean_skill   = if (length(skills)) mean(skills) else -Inf)
}

# -----------------------------
# Candidate counting & enumeration (under constraints)
# -----------------------------
season_splits <- list(
  S2A = c("SPRING","SUMMER"),   # Apr–Jun, Jul–Oct
  S1  = c("ALL")                # Apr–Oct
)

count_candidates <- function(){
  s1_gdd  <- length(gdd_base_grid) * length(hdd_cap_grid)        # 14
  s1_bins <- 1
  s2_gdd  <- (length(gdd_base_grid) * length(hdd_cap_grid))^2    # 196
  s2_bins <- 1
  s1_gdd + s1_bins + s2_gdd + s2_bins
}

# -----------------------------
# CV runner
# -----------------------------
run_cv_multiseason <- function(yld_raw, crop_slug, save_dir, years_use,
                               keep_states = filt_states_north,
                               verbose = TRUE){
  
  # Restrict outcome
  Y <- copy(as.data.table(yld_raw))
  Y <- Y[year %in% years_use]
  Y <- ensure_state_ansi(Y)
  if (!is.null(keep_states)) Y <- Y[State.ANSI %in% keep_states]
  Y <- Y[, .(fips, year, yield, State.ANSI)]
  if (nrow(Y) < 120) stop("Too few rows after filtering for ", crop_slug)
  
  total_candidates <- count_candidates()
  if (verbose){
    cat(sprintf("\n######### %s — CV (S1 & S2A only) ##########\n", toupper(crop_slug)))
    cat("Total candidates: ", total_candidates, "\n", sep = "")
  }
  
  recs <- vector("list", total_candidates)
  idx  <- 0L
  
  # -------- S1 (ALL): GDD variants + Bins --------
  # GDD: bases {5,10}, caps {26..32}
  for (b in gdd_base_grid){
    for (H in hdd_cap_grid){
      ta  <- list(ALL = paste0("gdd", b, "_", H))
      fts <- build_candidate_features("S1", ta)
      df  <- merge(Y, fts, by = c("fips","year"), all.x = TRUE)
      df[is.na(df)] <- 0
      
      rhs <- build_rhs_terms(names(df), "S1", ta)
      sc  <- loyo_skill_df(crop_slug, df, rhs)
      
      idx <- idx + 1L
      recs[[idx]] <- data.table(
        split = "S1",
        temp_model = "gdd",
        gdd_base_ALL = b,
        hdd_cap_ALL  = H,
        rhs = rhs,
        median_skill = sc$median_skill,
        mean_trim2   = sc$mean_trim2,
        mean_skill   = sc$mean_skill,
        n_folds      = length(sc$skills)
      )
      if (verbose && (idx %% 50 == 0)) cat(sprintf("  ... %d / %d\n", idx, total_candidates))
    }
  }
  # Bins (ALL)
  {
    ta  <- list(ALL = "bin5C_10")
    fts <- build_candidate_features("S1", ta)
    df  <- merge(Y, fts, by = c("fips","year"), all.x = TRUE)
    df[is.na(df)] <- 0
    
    rhs <- build_rhs_terms(names(df), "S1", ta)
    sc  <- loyo_skill_df(crop_slug, df, rhs)
    
    idx <- idx + 1L
    recs[[idx]] <- data.table(
      split = "S1",
      temp_model = "bins",
      rhs = rhs,
      median_skill = sc$median_skill,
      mean_trim2   = sc$mean_trim2,
      mean_skill   = sc$mean_skill,
      n_folds      = length(sc$skills)
    )
    if (verbose && (idx %% 50 == 0)) cat(sprintf("  ... %d / %d\n", idx, total_candidates))
  }
  
  # -------- S2A (SPRING,SUMMER): GDD×GDD + Bins×Bins --------
  # GDD x GDD with independent thresholds
  for (bS in gdd_base_grid){
    for (HS in hdd_cap_grid){
      for (bU in gdd_base_grid){
        for (HU in hdd_cap_grid){
          ta <- list(SPRING = paste0("gdd", bS, "_", HS),
                     SUMMER = paste0("gdd", bU, "_", HU))
          fts <- build_candidate_features("S2A", ta)
          df  <- merge(Y, fts, by = c("fips","year"), all.x = TRUE)
          df[is.na(df)] <- 0
          
          rhs <- build_rhs_terms(names(df), "S2A", ta)
          sc  <- loyo_skill_df(crop_slug, df, rhs)
          
          idx <- idx + 1L
          recs[[idx]] <- data.table(
            split = "S2A",
            temp_model = "gdd",
            gdd_base_SPRING = bS,
            hdd_cap_SPRING  = HS,
            gdd_base_SUMMER = bU,
            hdd_cap_SUMMER  = HU,
            rhs = rhs,
            median_skill = sc$median_skill,
            mean_trim2   = sc$mean_trim2,
            mean_skill   = sc$mean_skill,
            n_folds      = length(sc$skills)
          )
          if (verbose && (idx %% 50 == 0)) cat(sprintf("  ... %d / %d\n", idx, total_candidates))
        }
      }
    }
  }
  # Bins × Bins
  {
    ta  <- list(SPRING = "bin5C_10", SUMMER = "bin5C_10")
    fts <- build_candidate_features("S2A", ta)
    df  <- merge(Y, fts, by = c("fips","year"), all.x = TRUE)
    df[is.na(df)] <- 0
    
    rhs <- build_rhs_terms(names(df), "S2A", ta)
    sc  <- loyo_skill_df(crop_slug, df, rhs)
    
    idx <- idx + 1L
    recs[[idx]] <- data.table(
      split = "S2A",
      temp_model = "bins",
      rhs = rhs,
      median_skill = sc$median_skill,
      mean_trim2   = sc$mean_trim2,
      mean_skill   = sc$mean_skill,
      n_folds      = length(sc$skills)
    )
    if (verbose && (idx %% 50 == 0)) cat(sprintf("  ... %d / %d\n", idx, total_candidates))
  }
  
  # Compact results
  cv_table <- rbindlist(recs[seq_len(idx)], use.names = TRUE, fill = TRUE)
  setorder(cv_table, -median_skill, -mean_trim2, -mean_skill)
  
  # Save
  out_csv <- file.path(save_dir, sprintf("%s_cv_comparison_table.csv", tolower(crop_slug)))
  fwrite(cv_table, out_csv)
  if (verbose){
    cat("\n=== Best (by median skill) ===\n")
    print(cv_table[1L])
    cat(sprintf("Saved CV table → %s\n", out_csv))
  }
  
  invisible(cv_table)
}

# -----------------------------
# Outcome loader (same as before)
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

# ============================================================
# PRECOMPUTE (ONCE), THEN RUN FOR CROPS
# ============================================================
precompute_dd(tempbinsdata)    # base in {5,10}, cap in {26..32}
precompute_bin5(tempbinsdata)  # SPRING, SUMMER, ALL
precompute_ppt(weatherdata)    # SPRING, SUMMER, ALL (always quadratic later)

cat("\nSanity check: total candidate models per crop = ", count_candidates(), "\n", sep = "")  # 212

for (crop in c("sunflower","canola","flaxseed")) {
  yfile <- file.path(root, sprintf("%s_yields.csv", crop))
  if (!file.exists(yfile)) {
    cat(sprintf("\n-- %s yields not found: %s (skipping) --\n", toupper(crop), yfile))
    next
  }
  yld_raw <- load_outcome(yfile, years_use)
  run_cv_multiseason(yld_raw, crop_slug = crop, save_dir = root,
                     years_use = years_use, keep_states = filt_states_north,
                     verbose = TRUE)
}
