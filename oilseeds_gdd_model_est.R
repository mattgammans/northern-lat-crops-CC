# ============================================================
# CV + FE model with GDD-triggered stage splits (monthly data)
# Crops: sunflower, canola, flaxseed
# - CV chooses: 2 vs 3 stages; split1 (and split2); HDD/PPT sharing
# - Constraints: GDD slope shared; HDD threshold = 30C; PPT quadratic
# - Model: county FE + linear time trend; LOYO by year
# - Northern states filter; robust State.ANSI derivation
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

# -----------------------------
# PATHS
# -----------------------------
root <- "/Users/gammansm/Dropbox/NorthernLatitudeCrops_CCimpacts/Data"
bins_csv  <- file.path(root, "USco_PRISM_bins_monthly_1981_2022_cropland.csv")
weath_csv <- file.path(root, "USco_PRISM_weather_monthly_1981_2022_cropland.csv")

# -----------------------------
# Northern states only (ANSI)
# -----------------------------
filt_states_north <- c(8,16,26,27,30,31,38,39,41,53,56) # add 46 (SD), 55 (WI) if desired

# -----------------------------
# CV SETTINGS
# -----------------------------
years_use     <- 1981:2022
gdd_base_grid <- 10   # search bases (°C) used to compute GDD_mon
hdd_thresh    <- 30           # fixed HDD threshold (°C)
# Stage GDD cutpoints (°C·day) to search:
split1_grid   <- seq(600, 650, by = 50)   # begin flowering
split2_grid   <- seq(850, 900, by = 50)   # end flowering (3-stage only), requires > split1 + 100
hdd_patterns  <- c("shared_all", "separate_all")
ppt_patterns  <- c("shared_all", "separate_all")

# -----------------------------
# LOAD WEATHER
# -----------------------------
tempbinsdata <- fread(bins_csv)
weatherdata  <- fread(weath_csv)

stopifnot(all(c("fips","year","month") %in% names(tempbinsdata)))
stopifnot(all(c("fips","year","month") %in% names(weatherdata)))

tempbinsdata[, `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]
weatherdata[,  `:=`(fips = as.integer(fips), year = as.integer(year), month = as.integer(month))]

# precip column -> ppt
ppt_col <- intersect(c("ppt","ppt_mm","prcp_mm","precip_mm","prcp"), names(weatherdata))
if (!length(ppt_col)) stop("Monthly weather must include a precip column.")
setnames(weatherdata, ppt_col[1], "ppt")

# -----------------------------
# UTILITIES
# -----------------------------
ensure_state_ansi <- function(DT){
  if (!("State.ANSI" %in% names(DT)) || any(!is.finite(suppressWarnings(as.numeric(DT$State.ANSI))))) {
    if (!("fips" %in% names(DT))) stop("ensure_state_ansi(): 'fips' not found.")
    DT[, State.ANSI := as.integer(fips %/% 1000L)]
  } else {
    DT[, State.ANSI := as.integer(State.ANSI)]
  }
  DT
}

disc_bins <- function(TB){
  nm <- names(TB)
  bins <- nm[grepl("^bin_-?\\d+$", nm)]
  if (!length(bins)) stop("No bin_* columns in bins file.")
  k <- as.integer(sub("^bin_", "", bins))
  ord <- order(k)
  list(bins = bins[ord], k = k[ord])
}

# monthly GDD & HDD30 from discrete temp bins
make_monthly_dd <- function(tempbinsdata, gdd_base, heat_thr){
  db  <- disc_bins(tempbinsdata)
  MM  <- as.data.frame(tempbinsdata[, c("fips","year","month", db$bins), with = FALSE])
  Hdy <- as.matrix(MM[, db$bins, drop = FALSE]) / 24  # hours -> days
  w_gdd <- pmax(pmin(db$k, heat_thr) - gdd_base, 0)
  w_hdd <- pmax(db$k - heat_thr, 0)
  out <- MM[, c("fips","year","month"), drop = FALSE]
  out$GDD_mon   <- as.numeric(Hdy %*% w_gdd)
  out$HDD30_mon <- as.numeric(Hdy %*% w_hdd)
  out
}

# fractional share of a month before trigger T* given cum_prev and the month’s GDD
frac_before_in_month <- function(cum_prev, gdd_m, Tstar){
  if (is.na(gdd_m) || gdd_m <= 0) return(1)
  need <- Tstar - cum_prev
  if (need <= 0) return(0)      # already crossed
  if (need >= gdd_m) return(1)  # won't cross this month
  max(0, min(1, need / gdd_m))
}

# Build 2- or 3-stage features via GDD triggers with fractional month apportionment
# Outputs county-year: GDD_total; HDD_s1..s3; P_s1..s3; State.ANSI
build_stage_features <- function(tempbinsdata, weatherdata, gdd_base, split1, split2 = NA_real_, n_stages = 2){
  stopifnot(n_stages %in% c(2,3))
  DD <- make_monthly_dd(tempbinsdata, gdd_base = gdd_base, heat_thr = hdd_thresh)
  W  <- weatherdata[, .(fips, year, month, ppt)]
  W[, ppt := as.numeric(ppt) / 100]
  M  <- merge(DD, W, by = c("fips","year","month"), all.x = TRUE)
  setDT(M); setorder(M, fips, year, month)
  M[is.na(ppt), ppt := 0]
  M[, GDD_cum := cumsum(GDD_mon), by = .(fips, year)]
  M[, `:=`(w_s1=0, w_s2=0, w_s3=0)]
  
  byfun <- function(df){
    df <- df[order(month)]
    cum_prev <- shift(df$GDD_cum, fill = 0)
    gdd_m    <- df$GDD_mon
    if (n_stages == 2){
      p1 <- mapply(frac_before_in_month, cum_prev, gdd_m, MoreArgs = list(Tstar = split1))
      df$w_s1 <- p1
      df$w_s2 <- 1 - p1
      df$w_s3 <- 0
    } else {
      p1 <- mapply(frac_before_in_month, cum_prev, gdd_m, MoreArgs = list(Tstar = split1))
      g_before1 <- p1 * gdd_m
      g_after1  <- (1 - p1) * gdd_m
      cum_prev2 <- cum_prev + g_before1
      p2 <- mapply(frac_before_in_month, cum_prev2, g_after1, MoreArgs = list(Tstar = split2))
      df$w_s1 <- p1
      df$w_s2 <- (1 - p1) * p2
      df$w_s3 <- (1 - p1) * (1 - p2)
    }
    df
  }
  M <- M[, byfun(.SD), by = .(fips, year)]
  
  # stage-weight monthly HDD and precip; GDD_total is shared
  for(s in c("s1","s2","s3")){
    ws <- paste0("w_", s)
    if (s != "s3" || n_stages == 3){
      M[, paste0("HDD_", s) := get(ws) * HDD30_mon]
      M[, paste0("P_",   s) := get(ws) * ppt]
    } else {
      M[, paste0("HDD_", s) := 0]
      M[, paste0("P_",   s) := 0]
    }
  }
  
  feat <- M[, .(
    GDD_total = sum(GDD_mon, na.rm = TRUE),
    HDD_s1 = sum(HDD_s1, na.rm = TRUE),
    HDD_s2 = sum(HDD_s2, na.rm = TRUE),
    HDD_s3 = sum(HDD_s3, na.rm = TRUE),
    P_s1   = sum(P_s1,   na.rm = TRUE),
    P_s2   = sum(P_s2,   na.rm = TRUE),
    P_s3   = sum(P_s3,   na.rm = TRUE)
  ), by = .(fips, year)]
  
  feat <- ensure_state_ansi(feat)
  feat[]
}

# Build model formula given stage count and sharing patterns
make_fml_builder <- function(n_stages, hdd_pattern, ppt_pattern){
  function(avail_names){
    rhs <- c("GDD_total", "t")  # shared GDD + linear trend
    # HDD terms
    if (hdd_pattern == "shared_all"){
      need <- c("HDD_s1","HDD_s2", if (n_stages==3) "HDD_s3")
      need <- need[need %in% avail_names]
      rhs  <- c(rhs, sprintf("I(%s)", paste(need, collapse = " + ")))
    } else {
      add <- c("HDD_s1","HDD_s2", if (n_stages==3) "HDD_s3")
      add <- add[add %in% avail_names]
      rhs <- c(rhs, add)
    }
    # Precip quadratic
    if (ppt_pattern == "shared_all"){
      need <- c("P_s1","P_s2", if (n_stages==3) "P_s3")
      need <- need[need %in% avail_names]
      rhs  <- c(rhs, sprintf("I(%s)", paste(need, collapse = " + ")))
      rhs  <- c(rhs, sprintf("I((%s)^2)", paste(need, collapse = " + ")))
    } else {
      for(s in c("s1","s2", if(n_stages==3) "s3")){
        p <- paste0("P_", s)
        if (p %in% avail_names){
          rhs <- c(rhs, p, sprintf("I(%s^2)", p))
        }
      }
    }
    as.formula(paste("yield ~", paste(rhs, collapse = " + "), "| fips"))
  }
}

# LOYO CV (year-only), county FE; cluster by State.ANSI
loyo_cv_fe <- function(df, fml_builder, min_test = 10, min_train = 30){
  df <- as.data.table(df)
  df <- ensure_state_ansi(df)
  yrs <- sort(unique(df$year))
  r2s <- numeric(0)
  for (yv in yrs){
    tr <- df[year != yv]
    te <- df[year == yv]
    if (nrow(tr) < min_train || nrow(te) < min_test) next
    t0 <- mean(tr$year, na.rm = TRUE)
    tr[, t := year - t0]
    te[, t := year - t0]
    fml <- fml_builder(names(tr))
    fit <- tryCatch(feols(fml, data = tr, cluster = ~ State.ANSI), error = function(e) NULL)
    if (is.null(fit)) next
    yhat <- tryCatch(predict(fit, newdata = te, fixef = FALSE), error = function(e) rep(NA_real_, nrow(te)))
    ok <- is.finite(yhat) & is.finite(te$yield)
    if (!any(ok)) next
    mse <- mean((te$yield[ok] - yhat[ok])^2)
    vte <- var(te$yield[ok])
    if (is.finite(mse) && is.finite(vte) && vte > 0) r2s <- c(r2s, 1 - mse/vte)
  }
  if (length(r2s) == 0) -Inf else mean(r2s)
}

# Load outcome (annual log yield)
load_outcome <- function(csv_path, years_use){
  dat <- fread(csv_path)
  setnames(dat,
           old = c("State ANSI","County ANSI","Year","VALUE","value","SUCROSE","Sucrose"),
           new = c("State.ANSI","County.ANSI","year","Value","Value","Value","Value"),
           skip_absent = TRUE)
  dat[, Value := as.numeric(gsub("[^0-9.\\-]", "", as.character(Value)))]
  dat <- dat[is.finite(Value) & Value > 0]
  # build/clean IDs
  if (!("fips" %in% names(dat))) {
    dat[, State.ANSI  := as.integer(State.ANSI)]
    dat[, County.ANSI := as.integer(County.ANSI)]
    dat[, fips := 1000L * State.ANSI + County.ANSI]
  } else dat[, fips := as.integer(fips)]
  dat[, `:=`(year = as.integer(year), yield = log(Value))]
  dat <- dat[is.finite(fips) & year %in% years_use]
  dat <- ensure_state_ansi(dat)
  as.data.frame(dat)
}

# -----------------------------
# CV SEARCH WRAPPER
# -----------------------------
cv_search <- function(yld_raw, tempbinsdata, weatherdata){
  recs <- list(); best <- list(score = -Inf); idx <- 0L
  
  # candidate split sets
  cand_specs <- list()
  for (s1 in split1_grid) cand_specs <- append(cand_specs, list(list(n_stages=2, split1=s1, split2=NA)))
  for (s1 in split1_grid) for (s2 in split2_grid) if (s2 > s1 + 100)
    cand_specs <- append(cand_specs, list(list(n_stages=3, split1=s1, split2=s2)))
  
  for (base in gdd_base_grid){
    for (cs in cand_specs){
      feats <- build_stage_features(tempbinsdata, weatherdata,
                                    gdd_base = base, split1 = cs$split1,
                                    split2 = cs$split2, n_stages = cs$n_stages)
      df <- merge(yld_raw[, c("fips","year","yield","State.ANSI")],
                  feats[,   c("fips","year","GDD_total","HDD_s1","HDD_s2","HDD_s3","P_s1","P_s2","P_s3","State.ANSI")],
                  by = c("fips","year"), all = FALSE, suffixes = c(".y",".x"))
      df <- as.data.table(df)
      # collapse any .x/.y State.ANSI
      if ("State.ANSI.x" %in% names(df) || "State.ANSI.y" %in% names(df)){
        if (!"State.ANSI" %in% names(df))
          df[, State.ANSI := fifelse(is.finite(as.numeric(`State.ANSI.x`)), as.integer(`State.ANSI.x`), as.integer(`State.ANSI.y`))]
        df[, c("State.ANSI.x","State.ANSI.y") := NULL]
      }
      df <- ensure_state_ansi(df)
      # filter north
      df <- df[State.ANSI %in% filt_states_north]
      if (nrow(df) < 120) next
      
      for (hpat in hdd_patterns){
        for (ppat in ppt_patterns){
          fml_builder <- make_fml_builder(cs$n_stages, hpat, ppat)
          sc <- loyo_cv_fe(df, fml_builder)
          idx <- idx + 1L
          recs[[idx]] <- data.frame(
            gdd_base = base, n_stages = cs$n_stages,
            split1 = cs$split1, split2 = ifelse(is.na(cs$split2), NA, cs$split2),
            hdd = hpat, ppt = ppat, score = sc
          )
          if (is.finite(sc) && sc > best$score){
            best <- list(score = sc, df = df, gdd_base = base,
                         n_stages = cs$n_stages, split1 = cs$split1, split2 = cs$split2,
                         hdd = hpat, ppt = ppat)
          }
        }
      }
    }
  }
  
  cv_table <- if (length(recs)) { out <- rbindlist(recs); out[order(-score)] } else data.frame()
  list(cv_table = cv_table, best = best)
}

# -----------------------------
# Final refit on best spec
# -----------------------------
refit_best <- function(best){
  fml_builder <- make_fml_builder(best$n_stages, best$hdd, best$ppt)
  t0 <- mean(best$df$year, na.rm = TRUE)
  best$df[, t := year - t0]
  fml <- fml_builder(names(best$df))
  feols(fml, data = best$df, cluster = ~ State.ANSI)
}

# -----------------------------
# Runner for a crop
# -----------------------------
run_pipeline <- function(yld_raw, label = "Outcome"){
  cat(sprintf("\n\n########## %s ##########\n", toupper(label)))
  res <- cv_search(yld_raw, tempbinsdata, weatherdata)
  cat("\n================ ", toupper(label), " — CV (stages × splits × sharing; northern states only) ================\n", sep = "")
  if (nrow(res$cv_table)) print(res$cv_table, row.names = FALSE) else cat("No valid CV results.\n")
  
  if (is.finite(res$best$score)){
    mod <- refit_best(res$best)
    b <- res$best
    cat("\n=== Best spec ===\n")
    cat(sprintf("GDD base: %g | Stages: %d | split1: %g | split2: %s | HDD: %s | PPT: %s\n",
                b$gdd_base, b$n_stages, b$split1,
                ifelse(is.na(b$split2), "NA", as.character(b$split2)),
                b$hdd, b$ppt))
    cat(sprintf("LOYO CV R^2: %.3f\n\n", b$score))
    print(summary(mod))
    invisible(list(cv_table = res$cv_table, best = b, model = mod))
  } else {
    cat("\nNo valid spec.\n")
    invisible(NULL)
  }
}

# =========================
# RUN FOR THE THREE OILSEEDS
# =========================
for (crop in c("sunflower","canola","flaxseed")) {
  yfile <- file.path(root, sprintf("%s_yields.csv", crop))
  if (!file.exists(yfile)) {
    cat(sprintf("\n-- %s yields not found: %s (skipping) --\n", toupper(crop), yfile))
    next
  }
  yld_raw <- load_outcome(yfile, years_use)
  run_pipeline(yld_raw, label = tools::toTitleCase(crop))
}
