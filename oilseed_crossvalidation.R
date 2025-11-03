# ============================================================
# Multi-season CV with constrained temperature specs
# Splits: S3 (4-6,7,8-10), S2A (4-6,7-10), S1 (4-10)
# Temp per season: none, gdd5_30, gdd10_30, gdd5_28, gdd10_28, bin5C_10
#  - Constraint: in a candidate, all non-NONE seasons share the same temp spec
# PPT per season: none, quad (independent by season)
# Scoring: SKILL vs baseline yield ~ t + i(State.ANSI, t) | fips
# ============================================================

suppressPackageStartupMessages({ library(data.table); library(fixest) })

# -----------------------------
# PATHS / SETTINGS (edit if needed)
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")
years_use <- c(1981:2019, 2021:2022)
filt_states_north <- c(8,16,26,27,30,31,38,39,41,53,56)

# -----------------------------
# LOAD WEATHER
# -----------------------------
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)
stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))
tempbinsdata[, `:=`(fips=as.integer(fips), year=as.integer(year), month=as.integer(month))]
weatherdata[,  `:=`(fips=as.integer(fips), year=as.integer(year), month=as.integer(month))]

# ============================================================
# Fast precompute + caching (weather features across crops)
# + baseline-prediction cache for LOYO
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

# -----------------------------
# Helpers
# -----------------------------
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

# Monthly GDD/HDD from histogram-of-temp columns (hours -> days)
make_monthly_gddhdd <- function(tempbins, base, cap){
  stopifnot(all(c("fips","year","month") %in% names(tempbins)))
  DB <- disc_bins(tempbins)
  TB <- as.data.table(tempbins)[, c("fips","year","month", DB$bins), with = FALSE]
  # numeric + hours -> days
  for (b in DB$bins) {
    TB[, (b) := suppressWarnings(as.numeric(get(b)))]
    TB[is.na(get(b)), (b) := 0]
    TB[, (b) := get(b) / 24]
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

# Collapse to 5°C groups for selected months; returns county-year wide table
collapse_bins_5C_by_months <- function(tempbins, months_keep, suffix){
  stopifnot(all(c("fips","year","month") %in% names(tempbins)))
  DB <- disc_bins(tempbins)
  TB <- as.data.table(tempbins)[month %in% months_keep, c("fips","year","month", DB$bins), with = FALSE]
  # numeric + hours -> days
  for (b in DB$bins) {
    TB[, (b) := suppressWarnings(as.numeric(get(b)))]
    TB[is.na(get(b)), (b) := 0]
    TB[, (b) := get(b) / 24]
  }
  # Map each fine bin to a 5°C group label
  g5 <- floor(DB$k / 5) * 5
  groups <- split(DB$bins, g5)  # list: names are group labels
  # Sum across member columns, then sum over months -> county-year
  summed <- TB[, {
    # per-row rowSums for each 5C group
    sums <- lapply(groups, function(cols) rowSums(.SD[, ..cols], na.rm = TRUE))
    as.data.table(sums)
  }, by = .(fips, year)]
  # Now collapse (if multiple months) to county-year totals
  wide <- summed[, lapply(.SD, sum, na.rm = TRUE), by = .(fips, year)]
  # Rename columns to bin5C_<g>_<suffix>
  gnames <- names(groups)
  setnames(wide, old = gnames, new = paste0("bin5C_", gnames, "_", suffix))
  wide[]
}

# -----------------------------
# Global caches live across crops
# -----------------------------
.dd_cache    <- new.env(parent = emptyenv())  # monthly GDD/HDD for (base, cap)
.bin5_cache  <- new.env(parent = emptyenv())  # 5°C bins per season window
.ppt_cache   <- new.env(parent = emptyenv())  # PPT per season window
.baseline_cache <- new.env(parent = emptyenv())  # LOYO baseline predictions per crop/fold

key_dd   <- function(base, cap)     paste0("dd::", base, "_", cap)
key_bins <- function(season_name)   paste0("bin5::", season_name)
key_ppt  <- function(season_name)   paste0("ppt::", season_name)
fold_key <- function(crop, yrs)     paste0(crop, "::", paste(sort(as.integer(yrs)), collapse=","))

# -----------------------------
# Precompute monthly DD/HDD for all (base, cap)
# -----------------------------
precompute_dd <- function(tempbins){
  for (base in c(5, 10)) {
    for (cap in c(28, 30, 32)) {
      k <- key_dd(base, cap)
      if (!exists(k, .dd_cache, inherits = FALSE)) {
        assign(k, make_monthly_gddhdd(tempbins, base, cap), .dd_cache)
      }
    }
  }
}

# -----------------------------
# Precompute 5°C bins per season window
# -----------------------------
precompute_bin5 <- function(tempbins){
  windows <- list(SPRING = 4:6, JULY = 7, LATE = 8:10, SUMMER = 7:10, ALL = 4:10)
  for (nm in names(windows)){
    k <- key_bins(nm)
    if (!exists(k, .bin5_cache, inherits = FALSE)) {
      assign(k, collapse_bins_5C_by_months(tempbins, windows[[nm]], suffix = tolower(nm)), .bin5_cache)
    }
  }
}

# -----------------------------
# Precompute PPT per season window (expects weather$ppt numeric)
# -----------------------------
precompute_ppt <- function(weather){
  stopifnot(all(c("fips","year","month","ppt") %in% names(weather)))
  W <- as.data.table(weather)
  W[, ppt := as.numeric(ppt)]; W[is.na(ppt), ppt := 0]
  windows <- list(SPRING = 4:6, JULY = 7, LATE = 8:10, SUMMER = 7:10, ALL = 4:10)
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
# Build features for a candidate spec using caches
# split: "S3" = SPRING/JULY/LATE; "S2A" = SPRING/SUMMER; "S1" = ALL
# temp_assign: named list per season -> one of
#   "none", "gdd5_28", "gdd5_30", "gdd10_28", "gdd10_30", "bin5C_10"
# ppt_assign:  named list per season -> one of "none", "quad"
# -----------------------------
build_candidate_features <- function(split, temp_assign, ppt_assign){
  if (split == "S3")  seas <- list(SPRING = 4:6, JULY = 7, LATE = 8:10)
  if (split == "S2A") seas <- list(SPRING = 4:6, SUMMER = 7:10)
  if (split == "S1")  seas <- list(ALL = 4:10)
  feats <- list()
  
  # Temperature blocks from caches
  for (sn in names(seas)){
    spec <- temp_assign[[sn]]
    if (is.null(spec) || is.na(spec) || spec == "none") next
    
    if (startsWith(spec, "gdd")) {
      base <- if (grepl("^gdd5",  spec)) 5 else 10
      cap  <- if (grepl("_30$",  spec)) 30 else 28
      DD   <- get(key_dd(base, cap), .dd_cache)
      # sum within this season
      agg <- DD[month %in% seas[[sn]], .(GDD = sum(GDD_mon), HDD = sum(HDD_mon)), by = .(fips, year)]
      setnames(agg, c("GDD","HDD"), c(paste0("GDD_", sn), paste0("HDD_", sn)))
      feats[[length(feats) + 1]] <- agg
      
    } else if (spec == "bin5C_10") {
      # pre-collapsed 5C groups for this season
      feats[[length(feats) + 1]] <- get(key_bins(sn), .bin5_cache)
    }
  }
  
  # PPT blocks from caches
  for (sn in names(seas)){
    pps <- ppt_assign[[sn]]
    if (is.null(pps) || is.na(pps) || pps == "none") next
    P <- get(key_ppt(sn), .ppt_cache)
    set(P, j = paste0("PPT_", sn, "_sq"), value = P[[paste0("PPT_", sn)]]^2)
    feats[[length(feats) + 1]] <- P
  }
  
  if (!length(feats)) return(data.table(fips = integer(), year = integer()))
  out <- Reduce(function(a, b) merge(a, b, by = c("fips","year"), all = TRUE), feats)
  out[is.na(out)] <- 0
  out[]
}

# -----------------------------
# Baseline prediction cache for LOYO
# -----------------------------
get_baseline_preds <- function(crop_slug, tr, te){
  k <- fold_key(crop_slug, unique(te$year))
  if (exists(k, .baseline_cache, inherits = FALSE)) return(get(k, .baseline_cache))
  base <- feols(yield ~ t + i(State.ANSI, t) | fips, data = tr, cluster = ~ State.ANSI)
  pb   <- predict(base, newdata = te, fixef = FALSE)
  assign(k, pb, .baseline_cache)
  pb
}


# ============================================================
# Candidate enumeration + LOYO CV loops + save outputs
# ============================================================

suppressPackageStartupMessages({ library(data.table); library(fixest) })

# ---------- Splits and menus ----------
season_splits <- list(
  S3  = c("SPRING","JULY","LATE"),  # Apr–Jun, Jul, Aug–Oct
  S2A = c("SPRING","SUMMER"),       # Apr–Jun, Jul–Oct
  S1  = c("ALL")                    # Apr–Oct
)

# Temperature menu (identical type across non‑NONE seasons)
temp_types <- c("none", "gdd5_30","gdd10_30","gdd5_28","gdd10_28","bin5C_10")

# PPT menu (per season)
ppt_menu <- c("none","quad")

# ---------- Utilities ----------
# Power set (list of all subsets) of a character vector 'x'
all_subsets <- function(x, include_empty = TRUE){
  out <- list(character(0))
  if (length(x) == 0L) return(if (include_empty) out else list())
  for (el in x){
    out <- c(out, lapply(out, function(s) c(s, el)))
  }
  if (!include_empty) out <- Filter(length, out)
  out
}

# Build RHS term string from a merged df (so we can see which columns exist)
build_rhs_terms <- function(df_names, split_name, temp_assign, ppt_assign){
  rhs <- character(0)
  
  # Temperature terms
  for (sn in names(temp_assign)){
    spec <- temp_assign[[sn]]
    if (is.null(spec) || spec == "none") next
    if (startsWith(spec, "gdd")) {
      for (nm in c(paste0("GDD_", sn), paste0("HDD_", sn))){
        if (nm %in% df_names) rhs <- c(rhs, nm)
      }
    } else if (spec == "bin5C_10") {
      suf <- tolower(sn)
      # keep 5°C groups with center >= 10 (bin5C_<center>_<suffix>)
      bn <- grep(paste0("^bin5C_(-?\\d+)_", suf, "$"), df_names, value = TRUE)
      if (length(bn)) {
        centers <- as.integer(sub(paste0("^bin5C_(-?\\d+)_", suf, "$"), "\\1", bn))
        keep <- centers >= 10
        bn <- bn[keep]
        # drop zero-variance columns if any
        # (we don't have the data here; defensive check happens at fit time)
        rhs <- c(rhs, sprintf("`%s`", bn))
      }
    }
  }
  
  # PPT terms
  for (sn in names(ppt_assign)){
    ps <- ppt_assign[[sn]]
    if (is.null(ps) || ps == "none") next
    v <- paste0("PPT_", sn)
    if (v %in% df_names){
      rhs <- c(rhs, v, sprintf("I(%s^2)", v))
    }
  }
  
  if (!length(rhs)) return("0") else paste(rhs, collapse = " + ")
}

# LOYO skill vs baseline (uses baseline cache)
loyo_skill_df <- function(crop_slug, df, rhs_string, min_train = 30, min_test = 10){
  DT <- as.data.table(df)
  yrs <- sort(unique(DT$year))
  skills <- numeric(0)
  
  for (yv in yrs){
    tr <- copy(DT[year != yv])
    te <- copy(DT[year == yv])
    if (nrow(tr) < min_train || nrow(te) < min_test) next
    
    # center time on TRAIN
    t0 <- mean(tr$year, na.rm = TRUE)
    tr[, t := year - t0]; te[, t := year - t0]
    
    # --- baseline (cached) ---
    pb <- get_baseline_preds(crop_slug, tr, te)
    
    # --- candidate ---
    fml <- as.formula(paste0("yield ~ t + i(State.ANSI, t) + ", rhs_string, " | fips"))
    fit <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI),
                    error = function(e) NULL)
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
  
  list(skills = skills,
       median_skill = if (length(skills)) median(skills) else -Inf,
       mean_skill   = if (length(skills)) mean(skills)   else -Inf)
}

# Count total candidates (for a split) given the rules
count_candidates_for_split <- function(seasons){
  # Temperature:
  #  - one "all NONE" candidate
  #  - for each non-NONE type (5 types), choose any non-empty subset of seasons for which it's applied
  n_temp <- 1 + 5 * (2^length(seasons) - 1)
  # PPT: each season picks independently in {none, quad}
  n_ppt  <- 2^length(seasons)
  n_temp * n_ppt
}

# ---------- CV runner ----------
run_cv_multiseason <- function(yld_raw, crop_slug, save_dir, years_use,
                               keep_states = filt_states_north,
                               verbose = TRUE){
  
  # restrict outcome
  Y <- copy(as.data.table(yld_raw))
  Y <- Y[year %in% years_use]
  Y <- ensure_state_ansi(Y)
  if (!is.null(keep_states)) Y <- Y[State.ANSI %in% keep_states]
  Y <- Y[, .(fips, year, yield, State.ANSI)]
  if (nrow(Y) < 120) stop("Too few rows after filtering for ", crop_slug)
  
  # pre-count total candidates
  totals <- sapply(season_splits, count_candidates_for_split)
  total_candidates <- sum(totals)
  
  if (verbose){
    cat(sprintf("\n######### %s — Multi-season CV ##########\n", toupper(crop_slug)))
    cat("Total candidates: ", total_candidates, "\n", sep = "")
  }
  
  recs <- vector("list", total_candidates)
  idx  <- 0L
  
  # Enumerate splits
  for (split_name in names(season_splits)){
    seasons <- season_splits[[split_name]]
    
    # --- Temperature assignment enumerations ---
    # 1) All NONE
    temp_assign_list <- list(setNames(as.list(rep("none", length(seasons))), seasons))
    
    # 2) For each non-NONE type, apply to each non-empty subset of seasons
    non_none_types <- setdiff(temp_types, "none")
    subsets <- all_subsets(seasons, include_empty = FALSE)
    for (tt in non_none_types){
      for (S in subsets){
        ta <- setNames(as.list(rep("none", length(seasons))), seasons)
        if (length(S)) for (sn in S) ta[[sn]] <- tt
        temp_assign_list[[length(temp_assign_list) + 1L]] <- ta
      }
    }
    
    # --- PPT assignment enumerations (independent per season) ---
    ppt_assign_list <- list()
    ppt_subsets <- all_subsets(seasons, include_empty = TRUE) # subset that receives 'quad'
    for (P in ppt_subsets){
      pa <- setNames(as.list(rep("none", length(seasons))), seasons)
      if (length(P)) for (sn in P) pa[[sn]] <- "quad"
      ppt_assign_list[[length(ppt_assign_list) + 1L]] <- pa
    }
    
    # --- Cross product of temp and ppt assignments ---
    for (ta in temp_assign_list){
      # Build features once for this temp assignment; PPT will add linear/quad terms only
      # but features for PPT are also pulled from caches, so we must rebuild per PPT pattern.
      for (pa in ppt_assign_list){
        
        # features from caches
        feats <- build_candidate_features(split = split_name,
                                          temp_assign = ta,
                                          ppt_assign  = pa)
        
        # Merge with outcome (all yields kept; features fill missing with 0 in builder)
        if (!nrow(feats)) next
        df <- merge(Y, feats, by = c("fips","year"), all.x = TRUE)
        df[is.na(df)] <- 0
        
        # Build RHS
        rhs <- build_rhs_terms(names(df), split_name, ta, pa)
        
        # Score via LOYO vs baseline
        sc  <- loyo_skill_df(crop_slug, df, rhs)
        
        idx <- idx + 1L
        row <- list(
          split        = split_name,
          rhs          = rhs,
          median_skill = sc$median_skill,
          mean_skill   = sc$mean_skill,
          n_folds      = length(sc$skills)
        )
        # add temp_* and ppt_* columns explicitly
        for (nm in names(ta)) row[[paste0("temp_", nm)]] <- ta[[nm]]
        for (nm in names(pa)) row[[paste0("ppt_",  nm)]] <- pa[[nm]]
        
        recs[[idx]] <- as.data.table(row)
        
        if (verbose && (idx %% 50 == 0)) cat(sprintf("  ... %d / %d\n", idx, total_candidates))
      }
    }
  }
  
  # Compact results
  cv_table <- rbindlist(recs[seq_len(idx)], use.names = TRUE, fill = TRUE)
  setorder(cv_table, -median_skill, -mean_skill)
  
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

# ---------- Outcome loader (same as your earlier helper) ----------
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
# Example: RUN across the three crops
# (make sure caches were precomputed once before this)
# ============================================================

root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")

# Load raw weather once
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)
stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))

# Normalize types
tempbinsdata[, `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
weatherdata[,  `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]

# Ensure ppt column is named 'ppt'
ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
setnames(weatherdata, ppt_col[1], "ppt")

# Precompute caches ONCE
precompute_dd(tempbinsdata)
precompute_bin5(tempbinsdata)
precompute_ppt(weatherdata)

# Years / states
years_use <- c(1981:2019, 2021:2022)
filt_states_north <- c(8,16,26,27,30,31,38,39,41,53,56)  # CO, ID, MI, MN, MT, NE, ND, OH, OR, WA, WY (etc.)

# Run
for (crop in c("sunflower","canola","flaxseed")) {
  yfile <- file.path(root, sprintf("%s_yields.csv", crop))
  if (!file.exists(yfile)) { cat(sprintf("\n-- %s yields not found: %s (skipping) --\n", toupper(crop), yfile)); next }
  yld_raw <- load_outcome(yfile, years_use)
  run_cv_multiseason(yld_raw, crop_slug = crop, save_dir = root, years_use = years_use,
                     keep_states = filt_states_north, verbose = TRUE)
}



