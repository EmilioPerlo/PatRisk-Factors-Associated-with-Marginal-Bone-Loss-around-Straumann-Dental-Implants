# =============================================================================
# MARGINAL BONE LOSS ANALYSIS — v7 COMPLETE
# "Radiographic Evaluation of Patient, Implant, Surgical, and Prosthetic
#  Risk Factors Associated with Marginal Bone Loss around Straumann Implants"
#
# This is the fully self-contained script merging v6 (all original sections)
# with v7 (all peer-review implementations). Nothing needs to be pasted in.
#
# CHANGELOG v6 -> v7
#  FIX  1: emmeans() — added nuisance= argument to resolve reference grid
#           overflow (>2M rows) that caused all emmeans calls to fail in v6.
#  IMPL 1: Multiple Imputation (mice) for abutment_height_bin in Model B.
#  IMPL 2: MNAR tipping-point sensitivity analysis.
#  IMPL 3: Tipping-point summary table and figure.
#  IMPL 4: Formal included-vs-excluded comparison tables (Model B + abt_height).
#  IMPL 5: GEE AR(1) + 3-way LMM / GEE(exch) / GEE(AR1) comparison.
#  IMPL 6: Collapsed ea_distal_bin3 (3 levels) — removes near-empty >=30° cell.
#  IMPL 7: Dual definition of abutment_angulation_present (A = any value,
#           B = >= clinical threshold), controlled via CONFIG.
#  IMPL 8: Reference category supplement table saved as CSV.
#
# INHERITED FROM v6 (all bugs already fixed):
#  BUG 1: dir.create before sink().
#  BUG 2: on.exit(sink()).
#  BUG 3: GEE id = as.integer(patient_id).
#  BUG 4: VIF excludes high-level factors.
#  BUG 5: FDR applied to all terms of multi-level factors.
#  BUG 6: implant_site_fdi excluded from Model A fixed effects.
#  BUG 7: EPV check documented; implant_family sensitivity.
# =============================================================================


# ===========================================================================
# 0.  CONFIG
# ===========================================================================

PATH_PATIENTS   <- "Patient_final.xlsx"
PATH_IMPLANTS   <- "Implant_final.xlsx"
PATH_PROSTHETIC <- "Prosthetic_final.xlsx"
PATH_MBL        <- "MBL_final.xlsx"

BASELINE_YEAR         <- 1
BASELINE_LABEL        <- "year 1 (carico protesico)"
MCID                  <- 0.5
OUTLIER_ABS_MM        <- 10
MAX_MISSING_FOR_MODEL <- 0.30
SYNTH_PATIENT_CUTOFF  <- 9999
set.seed(20260607)

# ── v7 additions ─────────────────────────────────────────────────────────────
MI_M          <- 20      # number of imputed datasets (20-50 recommended)
MI_MAXIT      <- 30      # mice iteration cycles
MNAR_DELTAS   <- c(0, 0.10, 0.20, 0.30, 0.50)  # MNAR shift values (mm)

# abutment_angulation_present definition (IMPL 7):
#   "A" = any non-missing measurement (v6 default, n~80)
#   "B" = angulation >= ABT_ANG_THRESHOLD degrees (clinical threshold)
ABT_ANG_DEFINITION <- "A"
ABT_ANG_THRESHOLD  <- 15


# ===========================================================================
# 1.  PACKAGES
# ===========================================================================

pkgs <- c("readxl","tidyverse","lme4","lmerTest","emmeans",
          "tableone","car","performance","broom.mixed","geepack",
          "patchwork","pwr","mice","miceadds")

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))


# ===========================================================================
# 2.  OUTPUT DIRECTORIES + SINK
# ===========================================================================

dir.create("outputs",    showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/MI", showWarnings = FALSE, recursive = TRUE)

sink("outputs/console_output_v7.txt", split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("=============================================================\n")
cat("  MBL ANALYSIS v7 — ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("=============================================================\n\n")


# ===========================================================================
# 3.  HELPER FUNCTIONS
# ===========================================================================

patient_num <- function(pid) as.integer(stringr::str_extract(pid, "(?<=^p)\\d+"))

report_fit_n <- function(fit, label) {
  g <- summary(fit)$ngrps
  cat(sprintf("  [%s] analysed obs = %d | %s\n", label, nobs(fit),
              paste(names(g), unname(g), sep = " = ", collapse = " | ")))
}

RE_LADDER_DEFAULT <- c(
  "(1 + year | patient_id) + (1 + year | patient_id:implant_id)",
  "(1 | patient_id) + (1 + year | patient_id:implant_id)",
  "(1 + year | patient_id) + (1 | patient_id:implant_id)",
  "(1 | patient_id) + (1 | patient_id:implant_id)",
  "(1 | patient_id)"
)

fit_re_ladder <- function(fixed_rhs, data, reml = TRUE,
                          ladder = RE_LADDER_DEFAULT, label = "model") {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
  for (re in ladder) {
    fm  <- as.formula(paste(fixed_rhs, "+", re))
    fit <- tryCatch(suppressWarnings(lmer(fm, data = data, REML = reml,
                                          control = ctrl)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    singular  <- lme4::isSingular(fit, tol = 1e-4)
    converged <- length(fit@optinfo$conv$lme4$messages) == 0
    if (!singular && converged) {
      attr(fit, "re_used") <- re
      cat(sprintf("  [%s] random structure used: %s\n", label, re))
      return(fit)
    }
  }
  fm  <- as.formula(paste(fixed_rhs, "+", tail(ladder, 1)))
  fit <- suppressWarnings(lmer(fm, data = data, REML = reml, control = ctrl))
  attr(fit, "re_used") <- paste0(tail(ladder, 1),
                                 "  (FORCED — richer structures singular/non-converged)")
  cat(sprintf("  [%s] WARNING: %s\n", label, attr(fit, "re_used")))
  fit
}

annotate_mcid <- function(tidy_df, mcid = MCID) {
  tidy_df %>%
    mutate(
      sig          = p.value < 0.05,
      exceeds_mcid = abs(estimate) >= mcid,
      verdict = case_when(
        !sig & !exceeds_mcid ~ "neither detected nor clinically meaningful",
        sig  & !exceeds_mcid ~ "statistically detected but < MCID (clinically negligible)",
        !sig &  exceeds_mcid ~ "point estimate >= MCID but not significant (check CI width)",
        sig  &  exceeds_mcid ~ "significant AND >= MCID"
      )
    )
}

usable_terms <- function(terms, data, max_miss = MAX_MISSING_FOR_MODEL) {
  keep <- character(0); dropped <- character(0)
  for (t in terms) {
    if (!t %in% names(data)) {
      dropped <- c(dropped, sprintf("%s (not in data)", t)); next
    }
    x    <- data[[t]]
    miss <- mean(is.na(x))
    if (miss > max_miss) {
      dropped <- c(dropped, sprintf("%s (%.0f%% missing)", t, 100*miss)); next
    }
    ok <- if (is.factor(x)) nlevels(droplevels(x[!is.na(x)])) >= 2
          else length(unique(x[!is.na(x)])) >= 2
    if (isTRUE(ok)) keep <- c(keep, t)
    else dropped <- c(dropped, sprintf("%s (single level)", t))
  }
  if (length(dropped))
    cat("  dropped:", paste(dropped, collapse = "; "), "\n")
  keep
}

count_fixed_params <- function(terms, data) {
  total <- 1L
  for (t in terms) {
    x <- data[[t]]
    if (is.factor(x)) total <- total + nlevels(droplevels(x[!is.na(x)])) - 1L
    else total <- total + 1L
  }
  total
}

# ── v7: Rubin's rules pooling for lmer fitted on MI datasets ──────────────
pool_lmer_mi <- function(tidy_list) {
  terms_all <- unique(unlist(lapply(tidy_list, function(x) x$term)))
  m <- length(tidy_list)
  map_dfr(terms_all, function(trm) {
    ests <- map_dbl(tidy_list, function(df) {
      v <- df$estimate[df$term == trm]; if (length(v) == 0) NA_real_ else v[1] })
    ses  <- map_dbl(tidy_list, function(df) {
      v <- df$std.error[df$term == trm]; if (length(v) == 0) NA_real_ else v[1] })
    valid <- !is.na(ests) & !is.na(ses)
    if (sum(valid) < 2)
      return(tibble(term = trm, estimate = NA_real_, std.error = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_,
                    p.value = NA_real_, mi_m_valid = sum(valid)))
    ests <- ests[valid]; ses <- ses[valid]; m_v <- sum(valid)
    Q_bar <- mean(ests); U_bar <- mean(ses^2); B <- var(ests)
    T_var <- U_bar + (1 + 1/m_v) * B; T_se <- sqrt(T_var)
    lambda <- (1 + 1/m_v) * B / T_var
    df_r   <- (m_v - 1) / lambda^2
    t_stat <- Q_bar / T_se
    p_val  <- 2 * pt(-abs(t_stat), df = df_r)
    tibble(term = trm, estimate = Q_bar, std.error = T_se,
           conf.low  = Q_bar - qt(0.975, df_r) * T_se,
           conf.high = Q_bar + qt(0.975, df_r) * T_se,
           p.value = p_val, mi_m_valid = m_v)
  })
}

fit_lmer_mi <- function(formula_str, imp_data_list, label = "MI model") {
  cat(sprintf("  Fitting %s on %d imputed datasets...\n", label, length(imp_data_list)))
  tidy_list <- map(seq_along(imp_data_list), function(i) {
    fit <- tryCatch(
      fit_re_ladder(formula_str, data = imp_data_list[[i]], reml = TRUE,
                    label = sprintf("%s [imp %d]", label, i)),
      error = function(e) { cat(sprintf("    imp %d failed: %s\n", i, e$message)); NULL })
    if (is.null(fit)) return(NULL)
    broom.mixed::tidy(fit, effects = "fixed", conf.int = FALSE) %>%
      filter(term != "(Intercept)")
  })
  tidy_list <- Filter(Negate(is.null), tidy_list)
  if (length(tidy_list) == 0) { cat("  ERROR: all imputed models failed.\n"); return(NULL) }
  cat(sprintf("  Pooling %d/%d successful models...\n", length(tidy_list), length(imp_data_list)))
  pool_lmer_mi(tidy_list)
}


# ===========================================================================
# 4.  IMPORT
# ===========================================================================

raw_patients   <- read_excel(PATH_PATIENTS,   sheet = "patients")
raw_implants   <- read_excel(PATH_IMPLANTS,   sheet = "implants")
raw_prosthetic <- read_excel(PATH_PROSTHETIC, sheet = "prosthetic")
raw_mbl_wide   <- read_excel(PATH_MBL,        sheet = "mbl_wide_yr")

raw_mbl <- raw_mbl_wide %>%
  pivot_longer(
    cols          = matches("^mbl_y\\d+_(mesial|distal)_mm$"),
    names_to      = c("year","surface"),
    names_pattern = "mbl_y(\\d+)_(mesial|distal)_mm",
    values_to     = "mbl_mm"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(mbl_mm))

raw_mbl <- raw_mbl %>% filter(patient_id %in% raw_patients$patient_id)

cat("=== IMPORT ===\n")
cat("Patients:", nrow(raw_patients),
    "| Implants:", nrow(raw_implants),
    "| MBL rows (long):", nrow(raw_mbl), "\n")

synth_ids <- raw_implants %>%
  mutate(pn = patient_num(patient_id)) %>%
  filter(pn >= SYNTH_PATIENT_CUTOFF) %>% pull(implant_id)

if (length(synth_ids) > 0) {
  cat("\n*** DATA-PROVENANCE FLAG ***\n")
  cat(sprintf("  Synthetic implants (id >= %d): %d / %d (%.1f%%)\n",
              SYNTH_PATIENT_CUTOFF, length(synth_ids), nrow(raw_implants),
              100 * length(synth_ids) / nrow(raw_implants)))
} else {
  cat("\n*** DATA-PROVENANCE FLAG: All data REAL ***\n\n")
}


# ===========================================================================
# 5.  CLEAN: PATIENTS
# ===========================================================================

bin01 <- function(x) factor(x, levels = c(0,1), labels = c("No","Yes"))

df_patients <- raw_patients %>%
  filter(!is.na(patient_id)) %>%
  mutate(
    sex             = factor(sex, levels = c("Male","Female")),
    bisphosphonates = bin01(bisphosphonates),
    heart_disease   = bin01(heart_disease),
    hypertension    = bin01(hypertension),
    diabetes        = bin01(diabetes),
    periodontitis   = bin01(periodontitis),
    smoking         = bin01(smoking)
  ) %>%
  select(-any_of("notes"))


# ===========================================================================
# 6.  CLEAN: IMPLANTS
# ===========================================================================

IMPLANT_TYPE_LEVELS <- c("BL","BLT","BLX","NN","RN","TE","TL","TLC","TLX","WN")

df_implants <- raw_implants %>%
  filter(!is.na(implant_id)) %>%
  mutate(
    data_provenance  = factor(if_else(implant_id %in% synth_ids,"generated","real"),
                              levels = c("real","generated")),
    jaw              = factor(jaw,    levels = c("Maxilla","Mandible")),
    region           = factor(region, levels = c("Anterior","Posterior")),
    implant_level    = factor(implant_level,
                              levels = c("Bone Level","Tissue Level")),
    loading_protocol = factor(loading_protocol,
                              levels = c("Delayed","Immediate")),
    prosthesis_type  = factor(prosthesis_type,
                              levels = c("Crown","Bridge","Full-Arch")),
    surgical_protocol = factor(surgical_protocol,
                               levels = c("Type 4","Type 1","Type 2","Type 3")),
    implant_type = factor(
      if_else(implant_type %in% IMPLANT_TYPE_LEVELS, implant_type, NA_character_),
      levels = IMPLANT_TYPE_LEVELS
    ),
    implant_family = case_when(
      as.character(implant_type) %in% c("BL","BLT","BLX","NN")           ~ "Bone-Level",
      as.character(implant_type) %in% c("TL","TLC","TLX","WN","RN","TE") ~ "Tissue-Level",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Bone-Level","Tissue-Level")),
    implant_site_fdi = if ("implant_site_fdi" %in% names(pick(everything())))
                         factor(as.character(implant_site_fdi))
                       else NA_character_,
    gbr      = bin01(as.integer(gbr)),
    scais    = bin01(as.integer(scais)),
    flapless = bin01(as.integer(flapless)),
    length_mm            = as.numeric(length_mm),
    diameter_mm          = as.numeric(diameter_mm),
    age_at_placement     = as.numeric(age_at_placement),
    insertion_torque_ncm = as.numeric(insertion_torque_ncm),
    age_group = factor(
      case_when(
        age_at_placement <= 49 ~ "<=49",
        age_at_placement <= 70 ~ "50-70",
        age_at_placement >= 71 ~ ">=71",
        TRUE                   ~ NA_character_
      ), levels = c("<=49","50-70",">=71")),
    torque_group = factor(
      case_when(
        insertion_torque_ncm <= 20 ~ "<=20",
        insertion_torque_ncm <= 35 ~ "21-35",
        insertion_torque_ncm >= 36 ~ ">=36",
        TRUE                       ~ NA_character_
      ), levels = c("<=20","21-35",">=36"))
  ) %>%
  select(-any_of(c("implant_type_raw","age_note","insertion_note")))


# ===========================================================================
# 7.  CLEAN: PROSTHETIC
# ===========================================================================

df_prosthetic <- raw_prosthetic %>%
  filter(!is.na(implant_id)) %>%
  group_by(implant_id) %>%
  summarise(
    across(c(crown_implant_ratio,
             ea_mesial_deg, ea_distal_deg,
             outer_contour_mesial_mm, outer_contour_distal_mm,
             abutment_height_mm, abutment_angulation_deg,
             interprox_contact_mesial_mm, interprox_contact_distal_mm,
             implant_to_tooth_mesial_mm, implant_to_tooth_distal_mm,
             interimplant_dist_mesial_mm, interimplant_dist_distal_mm),
           ~ mean(.x, na.rm = TRUE)),
    ep_mesial = first(na.omit(ep_mesial)),
    ep_distal = first(na.omit(ep_distal)),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), ~ if_else(is.nan(.x), NA_real_, .x)),
    ep_mesial = factor(ep_mesial, levels = c("Straight","Concave","Convex")),
    ep_distal = factor(ep_distal, levels = c("Straight","Concave","Convex")),
    crown_implant_ratio_bin = factor(
      case_when(is.na(crown_implant_ratio) ~ NA_character_,
                crown_implant_ratio < 1    ~ "<1", TRUE ~ ">=1"),
      levels = c("<1",">=1")),
    ## 4-level mesial ea (kept — n=17 in >=30 is marginal but acceptable)
    ea_mesial_bin = factor(
      case_when(is.na(ea_mesial_deg) ~ NA_character_,
                ea_mesial_deg < 10   ~ "<10",
                ea_mesial_deg < 20   ~ "[10,20)",
                ea_mesial_deg < 30   ~ "[20,30)", TRUE ~ ">=30"),
      levels = c("<10","[10,20)","[20,30)",">=30")),
    ## Original 4-level distal ea — kept for Table 3 / descriptives only
    ea_distal_bin = factor(
      case_when(is.na(ea_distal_deg) ~ NA_character_,
                ea_distal_deg < 10   ~ "<10",
                ea_distal_deg < 20   ~ "[10,20)",
                ea_distal_deg < 30   ~ "[20,30)", TRUE ~ ">=30"),
      levels = c("<10","[10,20)","[20,30)",">=30")),
    ## IMPL 6: Collapsed 3-level distal ea — used in Model B v7
    ea_distal_bin3 = factor(
      case_when(is.na(ea_distal_deg) ~ NA_character_,
                ea_distal_deg < 10   ~ "<10",
                ea_distal_deg < 20   ~ "[10,20)", TRUE ~ ">=20"),
      levels = c("<10","[10,20)",">=20")),
    outer_contour_mesial_bin = factor(
      case_when(is.na(outer_contour_mesial_mm) ~ NA_character_,
                outer_contour_mesial_mm < 2    ~ "<2",
                outer_contour_mesial_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    outer_contour_distal_bin = factor(
      case_when(is.na(outer_contour_distal_mm) ~ NA_character_,
                outer_contour_distal_mm < 2    ~ "<2",
                outer_contour_distal_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    abutment_height_bin = factor(
      case_when(is.na(abutment_height_mm) ~ NA_character_,
                abutment_height_mm < 2    ~ "<2", TRUE ~ ">=2"),
      levels = c("<2",">=2")),
    ## IMPL 7: two definitions of angulation presence
    abt_ang_A = factor(
      if_else(!is.na(abutment_angulation_deg), "present", "absent"),
      levels = c("absent","present")),
    abt_ang_B = factor(
      case_when(is.na(abutment_angulation_deg)               ~ "absent",
                abutment_angulation_deg >= ABT_ANG_THRESHOLD  ~ "present",
                TRUE                                           ~ "absent"),
      levels = c("absent","present")),
    ## operative definition selected by CONFIG
    abutment_angulation_present = if (ABT_ANG_DEFINITION == "A") abt_ang_A else abt_ang_B,
    interprox_contact_mesial_bin = factor(
      case_when(is.na(interprox_contact_mesial_mm) ~ NA_character_,
                interprox_contact_mesial_mm < 2    ~ "<2",
                interprox_contact_mesial_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    interprox_contact_distal_bin = factor(
      case_when(is.na(interprox_contact_distal_mm) ~ NA_character_,
                interprox_contact_distal_mm < 2    ~ "<2",
                interprox_contact_distal_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    implant_to_tooth_mesial_bin = factor(
      case_when(is.na(implant_to_tooth_mesial_mm) ~ NA_character_,
                implant_to_tooth_mesial_mm < 2    ~ "<2",
                implant_to_tooth_mesial_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    implant_to_tooth_distal_bin = factor(
      case_when(is.na(implant_to_tooth_distal_mm) ~ NA_character_,
                implant_to_tooth_distal_mm < 2    ~ "<2",
                implant_to_tooth_distal_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    interimplant_dist_mesial_bin = factor(
      case_when(is.na(interimplant_dist_mesial_mm) ~ NA_character_,
                interimplant_dist_mesial_mm < 2    ~ "<2",
                interimplant_dist_mesial_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4")),
    interimplant_dist_distal_bin = factor(
      case_when(is.na(interimplant_dist_distal_mm) ~ NA_character_,
                interimplant_dist_distal_mm < 2    ~ "<2",
                interimplant_dist_distal_mm < 4    ~ "[2,4)", TRUE ~ ">=4"),
      levels = c("<2","[2,4)",">=4"))
  )

cat("=== PROSTHETIC: high-missingness variable counts ===\n")
cat("abutment_height_bin:\n");         print(table(df_prosthetic$abutment_height_bin, useNA="always"))
cat("\nabutment_angulation_present (Definition ", ABT_ANG_DEFINITION, "):\n", sep="")
print(table(df_prosthetic$abutment_angulation_present, useNA="always"))
cat(sprintf("\nDefinition A (any non-missing): present=%d  |  Definition B (>=%d°): present=%d\n\n",
            sum(df_prosthetic$abt_ang_A == "present", na.rm=TRUE),
            ABT_ANG_THRESHOLD,
            sum(df_prosthetic$abt_ang_B == "present", na.rm=TRUE)))


# ===========================================================================
# 8.  INTER-READER RELIABILITY — ICC(2,1)
# ===========================================================================

cat("=== RELIABILITY: ICC(2,1) ===\n")
rel_df <- raw_mbl %>%
  filter(!is.na(mbl_mm), abs(mbl_mm) <= OUTLIER_ABS_MM) %>%
  mutate(target    = interaction(implant_id, year, surface, drop = TRUE),
         replicate = factor(replicate))

icc21 <- tryCatch({
  m_rel <- lmer(mbl_mm ~ 1 + (1 | target) + (1 | replicate),
                data = rel_df, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa"))
  vc  <- as.data.frame(VarCorr(m_rel))
  v_t <- vc$vcov[vc$grp == "target"]
  v_r <- vc$vcov[vc$grp == "replicate"]
  v_e <- vc$vcov[vc$grp == "Residual"]
  v_t / (v_t + v_r + v_e)
}, error = function(e) { cat("  ICC model failed:", e$message, "\n"); NA_real_ })
cat(sprintf("  ICC(2,1) = %.3f  (threshold >= 0.75)\n\n", icc21))


# ===========================================================================
# 9.  PROCESS MBL — change score vs baseline
# ===========================================================================

df_cbl <- raw_mbl %>%
  filter(!is.na(mbl_mm)) %>%
  mutate(outlier = abs(mbl_mm) > OUTLIER_ABS_MM,
         mbl_mm  = if_else(outlier, NA_real_, mbl_mm)) %>%
  group_by(patient_id, implant_id, year) %>%
  summarise(cbl_mesial = mean(mbl_mm[surface == "mesial"], na.rm = TRUE),
            cbl_distal = mean(mbl_mm[surface == "distal"], na.rm = TRUE),
            n_outliers = sum(outlier, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(cbl_mesial,cbl_distal), ~ if_else(is.nan(.x), NA_real_, .x)),
         cbl_mean = rowMeans(cbind(cbl_mesial,cbl_distal), na.rm = TRUE),
         cbl_mean = if_else(is.nan(cbl_mean), NA_real_, cbl_mean))

cbl_baseline <- df_cbl %>%
  filter(year == BASELINE_YEAR) %>%
  transmute(implant_id, cbl_base = cbl_mean)

df_cbl <- df_cbl %>%
  left_join(cbl_baseline, by = "implant_id") %>%
  mutate(mbl_change = cbl_mean - cbl_base)

n_no_base <- df_cbl %>% filter(is.na(cbl_base)) %>% distinct(implant_id) %>% nrow()
if (n_no_base > 0)
  cat(sprintf("NOTE: %d implant(s) lack baseline reading -> excluded from change-score models.\n\n",
              n_no_base))

cat("CBL records:", nrow(df_cbl),
    "| Implants:", n_distinct(df_cbl$implant_id),
    "| Outliers removed:", sum(df_cbl$n_outliers > 0, na.rm = TRUE), "\n\n")


# ===========================================================================
# 10. ASSEMBLE ANALYSIS FRAMES
# ===========================================================================

n_imp_pt <- df_implants %>% count(patient_id, name = "n_implants")

df_A <- df_cbl %>%
  left_join(df_implants, by = c("implant_id","patient_id")) %>%
  left_join(df_patients, by = "patient_id") %>%
  left_join(n_imp_pt,   by = "patient_id") %>%
  mutate(multi_implant = factor(if_else(n_implants > 1,"Multiple","Single"),
                                levels = c("Single","Multiple")),
         patient_id = factor(patient_id),
         implant_id = factor(implant_id))

df_B <- df_A %>% inner_join(df_prosthetic, by = "implant_id")

cat("=== ANALYSIS FRAMES ===\n")
cat("  df_A:", nrow(df_A), "rows |", n_distinct(df_A$patient_id), "patients |",
    n_distinct(df_A$implant_id), "implants\n")
cat("  df_B:", nrow(df_B), "rows |", n_distinct(df_B$patient_id), "patients |",
    n_distinct(df_B$implant_id), "implants\n\n")


# ===========================================================================
# 11. PRE-SPECIFIED COVARIATE TIERS
# ===========================================================================

primary_tier <- c(
  "periodontitis","smoking","diabetes","bisphosphonates",
  "loading_protocol","surgical_protocol","implant_level"
)

secondary_tier <- c(
  "age_group","sex","hypertension","heart_disease",
  "jaw","region","prosthesis_type","implant_type",
  "length_mm","diameter_mm","gbr","scais","flapless",
  "torque_group","multi_implant"
)

covars_A <- c(primary_tier, secondary_tier)

interactions_A <- c(
  "year:implant_level",
  "implant_level:diameter_mm",
  "year:length_mm",
  "year:prosthesis_type"
)


# ===========================================================================
# 12. DESCRIPTIVES
# ===========================================================================

df_pt     <- df_A %>% distinct(patient_id, .keep_all = TRUE)
df_base_A <- df_A %>% filter(year == BASELINE_YEAR) %>% distinct(implant_id, .keep_all = TRUE)
df_base_B <- df_B %>% filter(year == BASELINE_YEAR) %>% distinct(implant_id, .keep_all = TRUE)

## Table 1: Patient characteristics
tab1 <- CreateTableOne(
  vars      = c("age_at_placement","age_group","sex",
                "bisphosphonates","heart_disease","hypertension",
                "diabetes","periodontitis","smoking","n_implants"),
  factorVars = c("age_group","sex","bisphosphonates","heart_disease",
                 "hypertension","diabetes","periodontitis","smoking"),
  data = df_pt
)
cat("===== TABLE 1: Patient characteristics =====\n")
print(tab1, showAllLevels = TRUE)

## Table 2: Implant characteristics
tab2 <- CreateTableOne(
  vars       = c("jaw","region","implant_level","loading_protocol",
                 "prosthesis_type","surgical_protocol","implant_type",
                 "implant_site_fdi","length_mm","diameter_mm",
                 "gbr","scais","flapless","age_group","torque_group"),
  factorVars = c("jaw","region","implant_level","loading_protocol","prosthesis_type",
                 "surgical_protocol","implant_type","implant_site_fdi",
                 "gbr","scais","flapless","age_group","torque_group"),
  data = df_base_A
)
cat("\n===== TABLE 2: Implant characteristics at baseline =====\n")
print(tab2, showAllLevels = TRUE)

## Table 3: Prosthetic characteristics
prosth_bin_vars <- c(
  "crown_implant_ratio_bin","ea_mesial_bin","ea_distal_bin",
  "ep_mesial","ep_distal","outer_contour_mesial_bin","outer_contour_distal_bin",
  "abutment_height_bin","abutment_angulation_present",
  "interprox_contact_mesial_bin","interprox_contact_distal_bin",
  "implant_to_tooth_mesial_bin","implant_to_tooth_distal_bin",
  "interimplant_dist_mesial_bin","interimplant_dist_distal_bin"
)
tab3 <- CreateTableOne(vars = prosth_bin_vars, factorVars = prosth_bin_vars, data = df_base_B)
cat("\n===== TABLE 3: Prosthetic characteristics at baseline =====\n")
cat("NOTE: abutment_height_bin: ~229 implants with value\n")
cat("      abutment_angulation_present: definition", ABT_ANG_DEFINITION, "(see CONFIG)\n\n")
print(tab3, showAllLevels = TRUE)

## Table 4: MBL change by year
tab4 <- df_A %>%
  group_by(year) %>%
  summarise(n_implants      = n_distinct(implant_id),
            mean_cbl        = round(mean(cbl_mean,   na.rm = TRUE), 3),
            mean_mbl_change = round(mean(mbl_change, na.rm = TRUE), 3),
            sd_mbl_change   = round(sd(mbl_change,   na.rm = TRUE), 3),
            .groups = "drop")
cat("\n===== TABLE 4: MBL change by year =====\n")
print(tab4, n = 25)


# ===========================================================================
# 13. MODEL A — PRIMARY
# ===========================================================================

cat("\n# ========== MODEL A — PRIMARY (outcome = mbl_change) ==========\n")

covars_A_ok <- usable_terms(covars_A, df_A)

n_params <- count_fixed_params(c("year", covars_A_ok), df_A)
n_units  <- n_distinct(df_A$implant_id)
epv      <- round(n_units / n_params, 1)
cat(sprintf("\n  EPV check: %d fixed params / %d implants -> EPV = %.1f\n",
            n_params, n_units, epv))
if (epv < 10)
  cat("  WARNING: EPV < 10 — use Model A (implant_family) from section 13b as primary.\n\n")

interactions_A_ok <- interactions_A[
  sapply(strsplit(interactions_A, ":"),
         function(t) all(t %in% c("year", covars_A_ok)))
]

fixed_rhs_A <- paste(
  "mbl_change ~ year +",
  paste(covars_A_ok, collapse = " + "),
  if (length(interactions_A_ok)) paste("+", paste(interactions_A_ok, collapse = " + ")) else ""
)

model_A <- fit_re_ladder(fixed_rhs_A, data = df_A, reml = TRUE, label = "Model A")
report_fit_n(model_A, "Model A")

tidy_A <- broom.mixed::tidy(model_A, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>% annotate_mcid()
cat("\n--- Model A fixed effects (95% CI, MCID-annotated) ---\n")
tidy_A %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, sig, exceeds_mcid, verdict) %>%
  print(n = Inf)

cat("\n--- Model A performance ---\n")
print(performance::r2(model_A))
print(performance::icc(model_A))

## VIF
vif_exclude <- c("implant_site_fdi","implant_type")
vif_vars    <- setdiff(covars_A_ok, vif_exclude)
vif_fit <- lm(
  as.formula(paste("mbl_change ~ year +", paste(vif_vars, collapse = " + "))),
  data = df_A %>% drop_na(any_of(c("mbl_change","year",vif_vars)))
)
cat(sprintf("\n--- VIF (excluded: %s; GVIF^(1/(2*Df)) > 2.5 = concern) ---\n",
            paste(vif_exclude, collapse=", ")))
tryCatch(print(car::vif(vif_fit)), error = function(e) cat("  VIF failed:", e$message, "\n"))

## Sensitivity: absolute CBL
cat("\n# ----- Sensitivity: absolute CBL adjusted for baseline -----\n")
fixed_A_abs <- paste(
  "cbl_mean ~ year + cbl_base +",
  paste(covars_A_ok, collapse = " + "),
  if (length(interactions_A_ok)) paste("+", paste(interactions_A_ok, collapse = " + ")) else ""
)
model_A_abs <- fit_re_ladder(fixed_A_abs, data = df_A, reml = TRUE, label = "Model A (abs)")
report_fit_n(model_A_abs, "Model A abs")
broom.mixed::tidy(model_A_abs, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>% annotate_mcid() %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, verdict) %>% print(n = Inf)


# ===========================================================================
# 13b. SENSITIVITY — implant_family (parsimonious primary model)
# ===========================================================================

cat("\n# ===== SENSITIVITY 13b: Model A with implant_family (2 levels) =====\n")
covars_A_fam <- covars_A_ok
if ("implant_type" %in% covars_A_fam)
  covars_A_fam <- c(setdiff(covars_A_fam, "implant_type"), "implant_family")
covars_A_fam <- usable_terms(covars_A_fam, df_A)

n_params_fam <- count_fixed_params(c("year", covars_A_fam), df_A)
cat(sprintf("  EPV (family model): %d params / %d implants -> EPV = %.1f\n",
            n_params_fam, n_units, round(n_units/n_params_fam,1)))

ints_fam <- interactions_A_ok[
  sapply(strsplit(interactions_A_ok, ":"),
         function(t) all(t %in% c("year", covars_A_fam)))
]
fixed_rhs_fam <- paste(
  "mbl_change ~ year +", paste(covars_A_fam, collapse = " + "),
  if (length(ints_fam)) paste("+", paste(ints_fam, collapse = " + ")) else ""
)
model_A_fam <- fit_re_ladder(fixed_rhs_fam, data = df_A, reml = TRUE,
                             label = "Model A (family)")
report_fit_n(model_A_fam, "Model A fam")

tidy_fam <- broom.mixed::tidy(model_A_fam, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>% annotate_mcid()
cat("\n--- Model A (family) fixed effects ---\n")
tidy_fam %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, sig, exceeds_mcid, verdict) %>%
  print(n = Inf)

cat("\n  Comparison PRIMARY: Model A (implant_type) vs Model A (implant_family):\n")
prim_pat <- paste(primary_tier, collapse="|")
tidy_A %>% filter(str_detect(term, prim_pat)) %>%
  select(term, est_type = estimate, p_type = p.value) %>%
  left_join(tidy_fam %>% filter(str_detect(term, prim_pat)) %>%
              select(term, est_fam = estimate, p_fam = p.value), by = "term") %>%
  mutate(delta = round(est_fam - est_type, 4)) %>% print()


# ===========================================================================
# 13c. EXPLORATORY — implant_site_fdi as random effect
# ===========================================================================

cat("\n# ===== EXPLORATORY 13c: variance attributed to implant_site_fdi =====\n")
if ("implant_site_fdi" %in% names(df_A) &&
    nlevels(droplevels(df_A$implant_site_fdi[!is.na(df_A$implant_site_fdi)])) >= 2) {
  site_null <- tryCatch(
    lmer(mbl_change ~ year + (1|patient_id) + (1|implant_site_fdi),
         data = df_A %>% filter(!is.na(implant_site_fdi)), REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("  site_FDI model failed:", e$message, "\n"); NULL }
  )
  if (!is.null(site_null)) {
    vc_site <- as.data.frame(VarCorr(site_null))
    cat("  Variance components (cross-classified):\n")
    print(vc_site %>% select(grp, vcov) %>%
            mutate(pct = round(100 * vcov / sum(vcov), 1)))
  }
} else {
  cat("  SKIP: implant_site_fdi not available or single level.\n")
}


# ===========================================================================
# 14. UNIVARIABLE SCREEN
# ===========================================================================

cat("\n# ===== UNIVARIABLE SCREEN (screening only — not for inference) =====\n")

run_univar <- function(pred, data, outcome = "mbl_change") {
  tryCatch({
    fm  <- as.formula(paste(outcome, "~", pred, "+ (1 | patient_id/implant_id)"))
    fit <- lmer(fm, data = data, REML = FALSE,
                control = lmerControl(optimizer = "bobyqa"))
    broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>% mutate(predictor = pred)
  }, error = function(e)
    tibble(predictor = pred, term = pred, estimate = NA_real_,
           conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
}

all_univar_preds <- unique(c("year", covars_A_ok, "implant_site_fdi", "implant_family"))
all_univar_preds <- all_univar_preds[all_univar_preds %in% names(df_A)]

univar_A <- map_dfr(all_univar_preds, run_univar, data = df_A)
univar_A %>% arrange(p.value) %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(predictor, term, estimate, conf.low, conf.high, p.value) %>%
  print(n = Inf)
cat("Reminder: inference based on multivariable Model A, not this screen.\n")


# ===========================================================================
# 15. MODEL B — PROSTHETIC (exploratory)
# ===========================================================================

cat("\n# ========== MODEL B — PROSTHETIC (exploratory) ==========\n")

core_for_B <- c("periodontitis","smoking","jaw","implant_level",
                "gbr","loading_protocol","prosthesis_type")

## IMPL 6: ea_distal_bin3 replaces ea_distal_bin in Model B
prosth_main <- c(
  "crown_implant_ratio_bin",
  "ea_mesial_bin","ea_distal_bin3",        # <-- collapsed distal (IMPL 6)
  "ep_mesial","ep_distal",
  "outer_contour_mesial_bin","outer_contour_distal_bin",
  "interprox_contact_mesial_bin","interprox_contact_distal_bin",
  "implant_to_tooth_mesial_bin","implant_to_tooth_distal_bin",
  "interimplant_dist_mesial_bin","interimplant_dist_distal_bin"
)

## Missingness report
all_B_check <- c(core_for_B, prosth_main, "abutment_height_bin","abutment_angulation_present")
miss_frac <- df_B %>%
  summarise(across(all_of(intersect(all_B_check, names(df_B))), ~mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "frac_missing")
cat("--- Missingness in df_B ---\n")
print(miss_frac %>% arrange(desc(frac_missing)) %>% mutate(frac_missing = round(frac_missing,3)))

## Spearman mesial/distal correlations
cat("\n--- Spearman mesial-distal correlations (collinearity check) ---\n")
pairs_md <- list(
  c("interprox_contact_mesial_bin","interprox_contact_distal_bin"),
  c("implant_to_tooth_mesial_bin", "implant_to_tooth_distal_bin"),
  c("interimplant_dist_mesial_bin","interimplant_dist_distal_bin"),
  c("outer_contour_mesial_bin",    "outer_contour_distal_bin"),
  c("ea_mesial_bin",               "ea_distal_bin3")
)
for (pr in pairs_md) {
  if (all(pr %in% names(df_B))) {
    x <- as.numeric(df_B[[pr[1]]]); y <- as.numeric(df_B[[pr[2]]])
    r <- cor(x, y, use = "pairwise.complete.obs", method = "spearman")
    cat(sprintf("  %s ~ %s: rho = %.3f\n", pr[1], pr[2], r))
  }
}

## Main Model B
terms_B <- usable_terms(c(core_for_B, prosth_main), df_B)

df_B_cc <- df_B %>%
  drop_na(any_of(c("mbl_change","year", terms_B))) %>%
  droplevels()
terms_B <- usable_terms(terms_B, df_B_cc)

cat(sprintf("\nModel B: %d variables | n_cc = %d obs / %d implants\n",
            length(terms_B), nrow(df_B_cc), n_distinct(df_B_cc$implant_id)))

fixed_rhs_B <- paste("mbl_change ~ year +", paste(terms_B, collapse = " + "))
model_B <- fit_re_ladder(fixed_rhs_B, data = df_B_cc, reml = TRUE, label = "Model B")
report_fit_n(model_B, "Model B")

tidy_B <- broom.mixed::tidy(model_B, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>% annotate_mcid()
cat("\n--- Model B fixed effects ---\n")
tidy_B %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, sig, exceeds_mcid, verdict) %>%
  print(n = Inf)

## VIF Model B
tryCatch({
  vif_B_vars <- setdiff(terms_B, c("implant_type","implant_site_fdi"))
  vif_B_fit  <- lm(as.formula(paste("mbl_change ~ year +",
                                     paste(vif_B_vars, collapse = " + "))), data = df_B_cc)
  cat("\n--- VIF Model B (GVIF^(1/(2*Df)); >2.5 = collinearity) ---\n")
  print(car::vif(vif_B_fit))
}, error = function(e) cat("  VIF Model B failed:", e$message, "\n"))


# ===========================================================================
# 15b. SUB-MODEL B_abt — abutment_height_bin, COMPLETE-CASE (n~229)
# ===========================================================================

cat("\n# ----- Sub-model B_abt: abutment_height_bin (complete-case, ~229 implants) -----\n")
cat("NOTE: estimates NOT comparable with full Model B. See Section 24 for MI version.\n\n")

core_abt    <- usable_terms(c(core_for_B,"abutment_height_bin"),
                             df_B %>% filter(!is.na(abutment_height_bin)))
df_abt      <- df_B %>% filter(!is.na(abutment_height_bin)) %>%
  drop_na(any_of(c("mbl_change","year",core_abt))) %>% droplevels()

cat(sprintf("  Subset abutment: %d implants | %d obs\n",
            n_distinct(df_abt$implant_id), nrow(df_abt)))

if (n_distinct(df_abt$implant_id) >= 20) {
  model_B_abt <- fit_re_ladder(
    paste("mbl_change ~ year +", paste(core_abt, collapse=" + ")),
    data = df_abt, reml = TRUE, label = "B_abt"
  )
  report_fit_n(model_B_abt, "B_abt")
  broom.mixed::tidy(model_B_abt, effects="fixed", conf.int=TRUE) %>%
    filter(term != "(Intercept)") %>% annotate_mcid() %>%
    mutate(across(where(is.numeric), ~round(.x,4))) %>%
    select(term, estimate, conf.low, conf.high, p.value, sig, exceeds_mcid, verdict) %>%
    print(n = Inf)
}


# ===========================================================================
# 15c. SUB-MODEL B_ang — abutment_angulation_present
# ===========================================================================

cat("\n# ----- Sub-model B_ang: abutment_angulation_present -----\n")
cat(sprintf("*** EXPLORATORY ONLY — Definition %s used (see CONFIG). ***\n\n",
            ABT_ANG_DEFINITION))

df_ang2 <- df_B %>%
  drop_na(any_of(c("mbl_change","year","abutment_angulation_present"))) %>%
  droplevels()

n_present <- sum(df_ang2$abutment_angulation_present[df_ang2$year==BASELINE_YEAR] == "present",
                 na.rm=TRUE)
cat(sprintf("  Implants 'present': %d / %d at baseline\n",
            n_present, n_distinct(df_ang2$implant_id)))

if (nlevels(droplevels(df_ang2$abutment_angulation_present)) >= 2) {
  model_B_ang <- tryCatch(
    lmer(mbl_change ~ year + abutment_angulation_present +
           (1 | patient_id/implant_id),
         data = df_ang2, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("  B_ang failed:", e$message, "\n"); NULL }
  )
  if (!is.null(model_B_ang)) {
    broom.mixed::tidy(model_B_ang, effects="fixed", conf.int=TRUE) %>%
      filter(term != "(Intercept)") %>% annotate_mcid() %>%
      mutate(across(where(is.numeric), ~round(.x,4))) %>%
      select(term, estimate, conf.low, conf.high, p.value, sig, exceeds_mcid, verdict) %>%
      print(n = Inf)
  }
} else {
  cat("  SKIP: abutment_angulation_present is single level in this subset.\n")
}


# ===========================================================================
# 16. HYPOTHESIS RESULTS
# ===========================================================================

cat("\n# ===== PRE-SPECIFIED HYPOTHESES (Model A) =====\n")

tidy_A %>% filter(str_detect(term, paste(primary_tier, collapse="|"))) %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, verdict) %>%
  { cat("--- PRIMARY tier ---\n"); print(., n=Inf); . }

tidy_A %>% filter(str_detect(term, paste(secondary_tier, collapse="|"))) %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, verdict) %>%
  { cat(sprintf("\n--- SECONDARY tier (%d terms) ---\n", nrow(.))); print(., n=Inf); . }

## BH-FDR on prosthetic variables
cat("\n--- EXPLORATORY: BH-FDR on prosthetic variables (q = 0.05) ---\n")
expl_preds <- intersect(prosth_main, names(df_B))

expl_list <- map(expl_preds, function(pr) {
  tryCatch({
    fit <- lmer(as.formula(paste("mbl_change ~", pr, "+ year + (1 | patient_id/implant_id)")),
                data = df_B, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    broom.mixed::tidy(fit, effects="fixed", conf.int=TRUE) %>%
      filter(str_detect(term, fixed(pr))) %>% mutate(predictor = pr)
  }, error = function(e) NULL)
})

expl <- bind_rows(expl_list)
if (nrow(expl) > 0) {
  expl <- expl %>%
    mutate(p_adj_BH = p.adjust(p.value, "BH"), sig_FDR = p_adj_BH < 0.05) %>%
    annotate_mcid() %>% arrange(p_adj_BH)
  expl %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
    select(predictor, term, estimate, conf.low, conf.high,
           p.value, p_adj_BH, sig_FDR, verdict) %>% print(n = Inf)
}


# ===========================================================================
# 17. EMMEANS — FIXED: nuisance= argument prevents reference grid overflow
# ===========================================================================

cat("\n# ===== ADJUSTED GROUP CONTRASTS (emmeans, Model A) =====\n")
cat("  NOTE v7 FIX: nuisance= argument added to all calls.\n\n")
emm_options(rg.limit = 2000000)
emm_at_years <- list(year = c(1,5,10))

## Build nuisance list: all covars_A_ok except the variable of interest
nuisance_all <- setdiff(covars_A_ok, "year")

emm_groups <- intersect(
  c("implant_level","gbr","surgical_protocol","loading_protocol",
    "prosthesis_type","implant_family","age_group","torque_group"),
  c(covars_A_ok, covars_A_fam)
)

for (grp in emm_groups) {
  cat("\n---", grp, "(adjusted) ---\n")
  ## choose model: fam model for implant_family, A for everything else
  mod_use <- if (grp == "implant_family") model_A_fam else model_A
  nuis    <- intersect(nuisance_all, setdiff(labels(terms(mod_use)), grp))

  res <- tryCatch({
    if (grp %in% c("implant_level","prosthesis_type","implant_family"))
      emmeans(mod_use, as.formula(paste("~", grp, "| year")),
              at = emm_at_years, nuisance = nuis)
    else
      emmeans(mod_use, as.formula(paste("~", grp)), nuisance = nuis)
  }, error = function(e) { cat("  emmeans failed:", e$message, "\n"); NULL })

  if (!is.null(res)) { print(res); print(pairs(res, adjust = "bonferroni")) }
}


# ===========================================================================
# 18. GEE FALLBACK — exchangeable (inherited from v6, BUG 3 fixed)
# ===========================================================================

cat("\n# ===== GEE FALLBACK (exchangeable, robust SE) =====\n")
gee_df <- df_A %>%
  drop_na(any_of(c("mbl_change","year", covars_A_ok))) %>%
  arrange(patient_id, implant_id, year) %>%
  mutate(patient_id_int = as.integer(patient_id))

gee_formula <- as.formula(paste("mbl_change ~ year +", paste(covars_A_ok, collapse = " + ")))

gee_fit <- tryCatch(
  geeglm(gee_formula, id = patient_id_int, data = gee_df,
         family = gaussian, corstr = "exchangeable"),
  error = function(e) { cat("  GEE exchangeable failed:", e$message, "\n"); NULL }
)
if (!is.null(gee_fit)) {
  broom::tidy(gee_fit, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>% annotate_mcid() %>%
    mutate(across(where(is.numeric), ~round(.x,4))) %>%
    select(term, estimate, conf.low, conf.high, p.value, verdict) %>% print(n = Inf)
  cat("Compare GEE vs LMM signs/CIs as robustness check.\n")
}


# ===========================================================================
# 19. SENSITIVITY: REAL-DATA ONLY
# ===========================================================================

cat("\n# ===== SENSITIVITY: Model A on REAL implants only =====\n")
df_A_real  <- df_A %>% filter(data_provenance == "real") %>% droplevels()
real_terms <- usable_terms(covars_A_ok, df_A_real)
real_ints  <- interactions_A_ok[
  sapply(strsplit(interactions_A_ok,":"), function(t) all(t %in% c("year",real_terms)))
]
model_A_real <- fit_re_ladder(
  paste("mbl_change ~ year +", paste(real_terms, collapse=" + "),
        if(length(real_ints)) paste("+", paste(real_ints, collapse=" + ")) else ""),
  data = df_A_real, reml = TRUE, label = "Model A (real only)"
)
report_fit_n(model_A_real, "Model A real")
broom.mixed::tidy(model_A_real, effects="fixed", conf.int=TRUE) %>%
  filter(term != "(Intercept)") %>% annotate_mcid() %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, conf.low, conf.high, p.value, verdict) %>% print(n=Inf)


# ===========================================================================
# 20. SENSITIVITY: EXCLUDE PATIENTS WITH >= 6 IMPLANTS
# ===========================================================================

cat("\n# ===== SENSITIVITY: Exclude patients with >= 6 implants =====\n")
ipp_sens <- df_A %>% distinct(patient_id, implant_id) %>% count(patient_id, name="n_imp")
df_A_low <- df_A %>% left_join(ipp_sens, by="patient_id") %>%
  filter(n_imp < 6) %>% droplevels()
cat(sprintf("  Reduced dataset: %d patients | %d implants | %d obs\n",
            n_distinct(df_A_low$patient_id), n_distinct(df_A_low$implant_id), nrow(df_A_low)))

lm_terms <- usable_terms(covars_A_ok, df_A_low)
lm_ints  <- interactions_A_ok[
  sapply(strsplit(interactions_A_ok,":"), function(t) all(t %in% c("year",lm_terms)))
]
model_A_low <- fit_re_ladder(
  paste("mbl_change ~ year +", paste(lm_terms, collapse=" + "),
        if(length(lm_ints)) paste("+", paste(lm_ints, collapse=" + ")) else ""),
  data = df_A_low, reml = TRUE, label = "Model A (no>=6)"
)
report_fit_n(model_A_low, "Model A no>=6")

tidy_lm <- broom.mixed::tidy(model_A_low, effects="fixed", conf.int=TRUE) %>%
  filter(term != "(Intercept)") %>% select(term, est_lm=estimate, p_lm=p.value)

compare_multi <- tidy_A %>%
  select(term, estimate, p.value) %>%
  left_join(tidy_lm, by="term") %>%
  mutate(delta     = round(est_lm - estimate, 4),
         pct_delta = round(100 * abs(delta) / abs(estimate), 1),
         sign_flip = !is.na(est_lm) & sign(estimate) != sign(est_lm),
         robust    = !sign_flip & abs(delta) < 0.02)

cat("\n  PRIMARY terms — full vs no>=6imp:\n")
compare_multi %>% filter(str_detect(term, paste(primary_tier, collapse="|"))) %>%
  mutate(across(where(is.numeric), ~round(.x,4))) %>%
  select(term, estimate, est_lm, delta, pct_delta, sign_flip, robust) %>% print(n=Inf)

unstable <- compare_multi %>% filter(!is.na(pct_delta), pct_delta > 20) %>% arrange(desc(pct_delta))
if (nrow(unstable) == 0) {
  cat("  No unstable terms (delta < 20%).\n")
} else {
  cat("\n  Terms with delta > 20%:\n")
  print(unstable %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
          select(term, estimate, est_lm, delta, pct_delta, sign_flip))
}


# ===========================================================================
# 21. DIAGNOSTICS
# ===========================================================================

cat("\n# ===== DIAGNOSTICS (Model A) =====\n")
diag_df <- tibble(fitted = fitted(model_A),
                  resid  = resid(model_A),
                  std    = resid(model_A) / sd(resid(model_A)))
cat(sprintf("  Residual SD: %.3f | scaled range: [%.1f, %.1f]\n",
            sd(diag_df$resid), min(diag_df$std), max(diag_df$std)))
cat("  Heavy tails / fan shape -> prefer GEE fallback for primary inference.\n")

p_diag <- patchwork::wrap_plots(
  ggplot(diag_df, aes(fitted, std)) +
    geom_point(alpha=.15, size=.8) +
    geom_hline(yintercept=0, linetype=2, colour="red") +
    geom_smooth(se=FALSE, method="loess", formula=y~x) +
    labs(title="Residuals vs fitted", x="Fitted", y="Std resid") + theme_minimal(),
  ggplot(diag_df, aes(sample=std)) +
    stat_qq(alpha=.2, size=.8) + stat_qq_line(colour="red") +
    labs(title="Q-Q") + theme_minimal()
)


# ===========================================================================
# 22. SAVE OUTPUTS (v6 figures + tables)
# ===========================================================================

cat("\n=== SAVING OUTPUTS ===\n")

## Figure 1: MBL over time
p_mbl <- ggplot(tab4 %>% filter(year<=15), aes(year, mean_mbl_change)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray50") +
  geom_line(color="#2c3e50", linewidth=1) + geom_point(color="#e74c3c", size=3) +
  geom_errorbar(aes(ymin=mean_mbl_change-(sd_mbl_change/sqrt(n_implants)),
                    ymax=mean_mbl_change+(sd_mbl_change/sqrt(n_implants))),
                width=0.2, color="#34495e") +
  scale_x_continuous(breaks=seq(1,15,1)) +
  labs(title="Figure 1: Mean MBL Change Over Time", subtitle="Error bars = SE",
       x="Follow-up Year", y="MBL Change (mm)") + theme_minimal()
ggsave("outputs/Fig1_MBL_over_time.png", p_mbl, width=8, height=5, dpi=300)

## Figure 2: Bone vs Tissue Level (using nuisance= fix)
emm_data <- tryCatch(
  as.data.frame(emmeans(model_A, ~ implant_level | year,
                        at       = list(year=1:10),
                        nuisance = intersect(nuisance_all,
                                             setdiff(covars_A_ok, "implant_level")))),
  error = function(e) { cat("  emmeans Fig2 failed:", e$message, "\n"); NULL }
)
if (!is.null(emm_data)) {
  ci_lo <- if("lower.CL" %in% names(emm_data)) "lower.CL" else "asymp.LCL"
  ci_hi <- if("upper.CL" %in% names(emm_data)) "upper.CL" else "asymp.UCL"
  names(emm_data)[names(emm_data)==ci_lo] <- "lower.CL"
  names(emm_data)[names(emm_data)==ci_hi] <- "upper.CL"
  p2 <- ggplot(emm_data, aes(year, emmean, color=implant_level, fill=implant_level)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50") +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), alpha=0.2, color=NA) +
    geom_line(linewidth=1) + geom_point(size=3) +
    scale_x_continuous(breaks=1:10) +
    labs(title="Figure 2: MBL by Implant Level (adjusted)",
         subtitle="Model A, nuisance-marginalised (95% CI)",
         x="Year", y="MBL Change (mm)",
         color="Implant Level", fill="Implant Level") +
    theme_minimal() + theme(legend.position="bottom")
  ggsave("outputs/Fig2_ImplantLevel_vs_Time.png", p2, width=8, height=5, dpi=300)
}

## Figure 3: Diagnostics
ggsave("outputs/Fig3_diagnostics_ModelA.png", p_diag, width=10, height=5, dpi=300)

## Tables CSV
write.csv(print(tab1, showAllLevels=TRUE, printToggle=FALSE),
          "outputs/Table1_Patient_Characteristics.csv")
write.csv(print(tab2, showAllLevels=TRUE, printToggle=FALSE),
          "outputs/Table2_Implant_Characteristics.csv")
write.csv(print(tab3, showAllLevels=TRUE, printToggle=FALSE),
          "outputs/Table3_Prosthetic_Characteristics.csv")
write.csv(tab4,      "outputs/Table4_MBL_change_by_year.csv",    row.names=FALSE)
write.csv(tidy_A,    "outputs/Table5_ModelA_Primary.csv",         row.names=FALSE)
write.csv(tidy_fam,  "outputs/Table5b_ModelA_family.csv",         row.names=FALSE)
write.csv(tidy_B,    "outputs/Table6_ModelB_Prosthetic.csv",      row.names=FALSE)
write.csv(univar_A,  "outputs/Table7_Univariable_Screening.csv",  row.names=FALSE)
if (exists("expl") && nrow(expl) > 0)
  write.csv(expl,    "outputs/Table8_Exploratory_FDR.csv",        row.names=FALSE)
write.csv(compare_multi %>% mutate(across(where(is.numeric), ~round(.x,4))),
          "outputs/Table9_Sensitivity_LowMulti.csv", row.names=FALSE)

cat("\n✓ v6 outputs saved.\n")


# ===========================================================================
# 23. POWER / MDE TABLE
# ===========================================================================

z_a <- qnorm(0.975); z_b <- qnorm(0.80)

mde_table <- tidy_A %>%
  filter(term %in% c(
    "year:implant_levelTissue Level",
    "periodontitisYes","smokingYes","bisphosphonatesYes","diabetesYes",
    "gbrYes","jawMandible",
    "prosthesis_typeBridge","prosthesis_typeFull-Arch",
    "age_group50-70","age_group>=71",
    "torque_group21-35","torque_group>=36"
  )) %>%
  mutate(
    n_ref = case_when(
      str_detect(term,"periodontitis|smoking|diabetes|bisphosphonates") ~ 96L,
      str_detect(term,"prosthesis")                                     ~ 150L,
      TRUE                                                              ~ 229L
    ),
    MDE_mm     = round(std.error * (z_a + z_b), 3),
    pwr_obs    = round(100 * pnorm(abs(estimate)/std.error - z_a), 0),
    below_MCID = MDE_mm < 0.5
  ) %>%
  select(term, n_ref, estimate, MDE_mm, pwr_obs, below_MCID)

cat("\n# ===== POWER / MDE TABLE =====\n")
print(mde_table)


# ===========================================================================
# 24. IMPL 4 — FORMAL INCLUDED vs EXCLUDED COMPARISONS
# ===========================================================================

cat("\n# ===== IMPL 4: Included vs Excluded in Model B =====\n")

compare_vars <- c("jaw","implant_level","prosthesis_type","loading_protocol",
                  "gbr","surgical_protocol","periodontitis","smoking",
                  "age_group","sex","length_mm","diameter_mm")

## 24a. All 501 implants: in full Model B vs not
df_B_check <- df_base_A %>%
  mutate(in_B = factor(implant_id %in% df_B_cc$implant_id,
                       levels=c(FALSE,TRUE),
                       labels=c("Excluded from Model B","Included in Model B")))

tab_inclB <- CreateTableOne(vars=compare_vars, strata="in_B",
                            data=df_B_check, test=TRUE, addOverall=TRUE)
cat("\n--- Included vs Excluded in Model B (n=501 at baseline) ---\n")
cat("NOTE: p<0.05 = potential selection bias in Model B.\n\n")
print(tab_inclB, showAllLevels=TRUE)
write.csv(print(tab_inclB, showAllLevels=TRUE, printToggle=FALSE),
          "outputs/Table_InclExcl_ModelB.csv")

## 24b. Among all 501 implants: abutment height observed vs missing
df_B_base <- df_B %>% filter(year==BASELINE_YEAR) %>%
  distinct(implant_id, .keep_all=TRUE)

df_abt_check <- df_B_base %>%
  mutate(in_abt = factor(!is.na(abutment_height_bin),
                         levels=c(FALSE,TRUE),
                         labels=c("abt_height MISSING","abt_height OBSERVED")))

tab_inclAbt <- CreateTableOne(vars=compare_vars, strata="in_abt",
                              data=df_abt_check, test=TRUE, addOverall=TRUE)
cat("\n--- Abutment Height: Observed vs Missing ---\n")
cat("NOTE: p<0.05 = missing abutment data is not MCAR — MI result is more reliable.\n\n")
print(tab_inclAbt, showAllLevels=TRUE)
write.csv(print(tab_inclAbt, showAllLevels=TRUE, printToggle=FALSE),
          "outputs/Table_InclExcl_AbutmentHeight.csv")


# ===========================================================================
# 25. IMPL 1 — MULTIPLE IMPUTATION FOR abutment_height_bin
# ===========================================================================

cat("\n# ===== IMPL 1: Multiple Imputation — abutment_height_bin =====\n")
cat(sprintf("  MI_M = %d | MI_MAXIT = %d\n\n", MI_M, MI_MAXIT))

## 25a. Implant-level imputation dataset
imp_base_vars <- c(
  "implant_id","patient_id","mbl_change",
  "periodontitis","smoking","jaw","implant_level",
  "gbr","loading_protocol","prosthesis_type",
  "abutment_height_mm",           # TARGET
  "crown_implant_ratio","ea_mesial_deg","ea_distal_deg",
  "outer_contour_mesial_mm","interprox_contact_mesial_mm","implant_to_tooth_mesial_mm"
)

df_imp_base <- df_B %>%
  filter(year == BASELINE_YEAR) %>%
  distinct(implant_id, .keep_all = TRUE) %>%
  select(any_of(imp_base_vars))

cat("  Missingness in imputation input:\n")
df_imp_base %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to="var", values_to="frac_missing") %>%
  filter(frac_missing > 0) %>% arrange(desc(frac_missing)) %>%
  mutate(frac_missing=round(frac_missing,3)) %>% print()

## 25b. Run mice
df_mice <- df_imp_base %>%
  select(-implant_id, -patient_id, -mbl_change) %>%
  mutate(across(where(is.factor), droplevels))

pred_matrix <- quickpred(df_mice, mincor = 0.1)
meth <- make.method(df_mice)
meth["abutment_height_mm"] <- "pmm"

cat("\n  Running mice...\n")
imp_obj <- mice(df_mice, m=MI_M, maxit=MI_MAXIT, method=meth,
                pred=pred_matrix, seed=20260607, printFlag=FALSE)
cat("  Imputation complete.\n")

## Convergence plot
png("outputs/MI/MI_convergence_abt_height.png", width=1400, height=700, res=150)
plot(imp_obj, c("abutment_height_mm"))
dev.off()
cat("  -> outputs/MI/MI_convergence_abt_height.png\n")

## Strip plot
png("outputs/MI/MI_stripplot_abt_height.png", width=1000, height=600, res=150)
stripplot(imp_obj, abutment_height_mm ~ .imp, pch=20, cex=1.2,
          main="Observed (blue) vs Imputed (red): abutment_height_mm")
dev.off()
cat("  -> outputs/MI/MI_stripplot_abt_height.png\n")

## 25c. Build imputed longitudinal datasets
missing_abt_ids <- df_imp_base %>% filter(is.na(abutment_height_mm)) %>% pull(implant_id)

imp_data_list <- lapply(seq_len(MI_M), function(i) {
  df_complete_base <- complete(imp_obj, action=i) %>%
    mutate(implant_id = df_imp_base$implant_id,
           abutment_height_bin_mi = factor(
             if_else(abutment_height_mm < 2, "<2", ">=2"), levels=c("<2",">=2"))) %>%
    select(implant_id, abutment_height_bin_mi)

  df_B %>%
    left_join(df_complete_base, by="implant_id") %>%
    mutate(abutment_height_bin = abutment_height_bin_mi) %>%
    select(-abutment_height_bin_mi) %>%
    drop_na(any_of(c("mbl_change","year",
                      usable_terms(c(core_for_B,"abutment_height_bin"), df_B))))
})

cat(sprintf("\n  Implants per imputed dataset: min=%d max=%d (expected ~501)\n",
            min(sapply(imp_data_list, function(d) n_distinct(d$implant_id))),
            max(sapply(imp_data_list, function(d) n_distinct(d$implant_id)))))

## 25d. Fit and pool
core_abt_mi <- usable_terms(c(core_for_B,"abutment_height_bin"), imp_data_list[[1]])
fixed_rhs_abt_mi <- paste("mbl_change ~ year +", paste(core_abt_mi, collapse=" + "))

cat("\n  --- Fitting lmer on each imputed dataset (B_abt_MI) ---\n")
pooled_abt_mi <- fit_lmer_mi(fixed_rhs_abt_mi, imp_data_list, label="B_abt_MI")

if (!is.null(pooled_abt_mi)) {
  pooled_abt_mi <- pooled_abt_mi %>% annotate_mcid()
  cat("\n--- Model B_abt_MI: Pooled estimates (Rubin's rules, 95% CI) ---\n")
  pooled_abt_mi %>% mutate(across(where(is.numeric), ~round(.x,4))) %>%
    select(term, estimate, conf.low, conf.high, p.value, mi_m_valid,
           sig, exceeds_mcid, verdict) %>% print(n=Inf)
  write.csv(pooled_abt_mi, "outputs/MI/Table_MI_ModelB_abt.csv", row.names=FALSE)
  cat("  -> outputs/MI/Table_MI_ModelB_abt.csv\n")
}


# ===========================================================================
# 26. IMPL 2 & 3 — MNAR TIPPING-POINT SENSITIVITY
# ===========================================================================

cat("\n# ===== IMPL 2&3: MNAR Tipping-Point Sensitivity =====\n")
cat(sprintf("  Deltas tested: %s mm\n\n", paste(MNAR_DELTAS, collapse=", ")))
cat(sprintf("  Implants with originally-missing abutment_height_mm: %d\n", length(missing_abt_ids)))

tipping_results <- map_dfr(MNAR_DELTAS, function(delta) {
  cat(sprintf("  delta = %+.2f mm...\n", delta))
  imp_shifted <- lapply(seq_len(MI_M), function(i) {
    df_complete_base <- complete(imp_obj, action=i) %>%
      mutate(
        implant_id = df_imp_base$implant_id,
        abutment_height_mm_shifted = if_else(
          implant_id %in% missing_abt_ids, abutment_height_mm + delta, abutment_height_mm),
        abutment_height_bin_mi = factor(
          if_else(abutment_height_mm_shifted < 2, "<2", ">=2"), levels=c("<2",">=2"))
      ) %>% select(implant_id, abutment_height_bin_mi)
    df_B %>%
      left_join(df_complete_base, by="implant_id") %>%
      mutate(abutment_height_bin = abutment_height_bin_mi) %>%
      select(-abutment_height_bin_mi) %>%
      drop_na(any_of(c("mbl_change","year", core_abt_mi)))
  })
  pooled_shifted <- fit_lmer_mi(fixed_rhs_abt_mi, imp_shifted,
                                label=sprintf("MNAR d=%+.2f", delta))
  if (is.null(pooled_shifted)) return(NULL)
  pooled_shifted %>% filter(str_detect(term, "abutment_height_bin")) %>%
    mutate(delta=delta)
})

cat("\n--- Tipping-Point Table ---\n")
if (!is.null(tipping_results) && nrow(tipping_results) > 0) {
  tipping_results <- tipping_results %>%
    mutate(conclusion_sig = p.value < 0.05, across(where(is.numeric), ~round(.x,4)))
  print(tipping_results %>% select(term, delta, estimate, conf.low, conf.high,
                                    p.value, conclusion_sig))

  base_est <- tipping_results %>% filter(delta==0) %>% pull(estimate)
  base_sig <- tipping_results %>% filter(delta==0) %>% pull(p.value) < 0.05
  if (length(base_est) > 0) {
    tipping_results <- tipping_results %>%
      mutate(sign_flip  = !is.na(estimate) & sign(estimate) != sign(base_est[1]),
             p_crossed  = (p.value < 0.05) != base_sig[1],
             is_tipping = sign_flip | p_crossed)
    tip_delta <- tipping_results %>% filter(is_tipping, delta!=0) %>%
      arrange(abs(delta)) %>% slice(1) %>% pull(delta)
    if (length(tip_delta)==0)
      cat(sprintf("\n  ROBUST: no tipping point found within [%s] mm range.\n",
                  paste(range(MNAR_DELTAS), collapse=" to ")))
    else
      cat(sprintf("\n  TIPPING POINT at delta = %+.2f mm.\n", tip_delta))
  }

  write.csv(tipping_results, "outputs/MI/Table_MNAR_TippingPoint.csv", row.names=FALSE)

  if (nrow(tipping_results) > 1) {
    p_tip <- ggplot(tipping_results,
                    aes(x=delta, y=estimate, ymin=conf.low, ymax=conf.high)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
      geom_ribbon(alpha=0.2, fill="#2C5F8A") +
      geom_line(colour="#2C5F8A", linewidth=1) +
      geom_point(aes(shape=conclusion_sig), size=3, colour="#2C5F8A") +
      scale_shape_manual(values=c("TRUE"=16,"FALSE"=1),
                         labels=c("TRUE"="p<0.05","FALSE"="p>=0.05"),
                         name="Significance") +
      labs(title="MNAR Tipping-Point: abutment_height_bin >=2mm",
           subtitle="delta = assumed shift (mm) for originally-missing values. delta=0 is MAR.",
           x="MNAR shift delta (mm)", y="Estimated effect on MBL change (mm)") +
      theme_minimal(base_size=12)
    ggsave("outputs/MI/Fig_MNAR_TippingPoint.png", p_tip, width=8, height=5, dpi=300)
    cat("  -> outputs/MI/Fig_MNAR_TippingPoint.png\n")
  }
}


# ===========================================================================
# 27. IMPL 5 — GEE AR(1) + 3-WAY COMPARISON
# ===========================================================================

cat("\n# ===== IMPL 5: GEE AR(1) + 3-way comparison =====\n")

gee_ar1 <- tryCatch(
  geeglm(gee_formula, id=patient_id_int, data=gee_df,
         family=gaussian, corstr="ar1"),
  error = function(e) { cat("  GEE AR(1) failed:", e$message, "\n"); NULL }
)

## Re-fit exchangeable for clean comparison
gee_exch <- tryCatch(
  geeglm(gee_formula, id=patient_id_int, data=gee_df,
         family=gaussian, corstr="exchangeable"),
  error = function(e) { cat("  GEE exch re-fit failed:", e$message, "\n"); NULL }
)

build_gee_tidy <- function(gfit, label) {
  if (is.null(gfit)) return(NULL)
  broom::tidy(gfit, conf.int=TRUE) %>% filter(term != "(Intercept)") %>%
    select(term, !!paste0("est_",label) := estimate, !!paste0("p_",label) := p.value)
}

tidy_lmm_cmp  <- tidy_fam %>% select(term, est_LMM=estimate, p_LMM=p.value)
tidy_exch_cmp <- build_gee_tidy(gee_exch, "GEE_exch")
tidy_ar1_cmp  <- build_gee_tidy(gee_ar1,  "GEE_ar1")

compare_3way <- tidy_lmm_cmp
if (!is.null(tidy_exch_cmp)) compare_3way <- left_join(compare_3way, tidy_exch_cmp, by="term")
if (!is.null(tidy_ar1_cmp))  compare_3way <- left_join(compare_3way, tidy_ar1_cmp,  by="term")

est_cols <- intersect(c("est_LMM","est_GEE_exch","est_GEE_ar1"), names(compare_3way))
if (length(est_cols) >= 2) {
  compare_3way <- compare_3way %>% rowwise() %>%
    mutate(
      est_range  = diff(range(c_across(all_of(est_cols)), na.rm=TRUE)),
      est_maxabs = max(abs(c_across(all_of(est_cols))), na.rm=TRUE),
      pct_spread = round(100 * est_range / (est_maxabs + 1e-9), 1),
      sign_agree = n_distinct(sign(c_across(all_of(est_cols))), na.rm=TRUE) == 1,
      flag       = pct_spread > 20 | !sign_agree
    ) %>% ungroup()
}

cat("\n--- 3-Way Comparison: LMM (family model) vs GEE(exch) vs GEE(AR1) ---\n")
compare_3way %>% mutate(across(where(is.numeric), ~round(.x,4))) %>% print(n=Inf)

flagged <- if ("flag" %in% names(compare_3way)) compare_3way %>% filter(flag) else tibble()
if (nrow(flagged) > 0) {
  cat(sprintf("\n  %d terms flagged (>20%% spread or sign disagreement):\n", nrow(flagged)))
  print(flagged %>% select(term, all_of(est_cols), pct_spread, sign_agree))
} else {
  cat("\n  No terms flagged — all three methods broadly agree.\n")
}

if (!is.null(gee_ar1))  cat(sprintf("  GEE AR(1) alpha: %.4f\n",   gee_ar1$geese$alpha))
if (!is.null(gee_exch)) cat(sprintf("  GEE Exch  alpha: %.4f\n",   gee_exch$geese$alpha))

write.csv(compare_3way %>% mutate(across(where(is.numeric), ~round(.x,4))),
          "outputs/Table_3way_LMM_GEE.csv", row.names=FALSE)
cat("  -> outputs/Table_3way_LMM_GEE.csv\n")


# ===========================================================================
# 28. IMPL 8 — REFERENCE CATEGORY SUPPLEMENT TABLE
# ===========================================================================

cat("\n# ===== IMPL 8: Reference Category Supplement Table =====\n")

ref_table <- tibble::tribble(
  ~Variable,                    ~Reference_Level,           ~Comparators,                                        ~Clinical_Rationale,
  "sex",                        "Male",                     "Female",                                            "Biological baseline; hormonal effects modelled as deviation from male",
  "age_group",                  "<=49",                     "50-70; >=71",                                       "Youngest patients as biologically ideal baseline",
  "bisphosphonates",            "No",                       "Yes",                                               "Drug-naive patients as baseline; bisphosphonates alter bone turnover",
  "heart_disease",              "No",                       "Yes",                                               "Absence of comorbidity as baseline",
  "hypertension",               "No",                       "Yes",                                               "Absence of comorbidity as baseline",
  "diabetes",                   "No",                       "Yes",                                               "Absence of comorbidity as baseline",
  "periodontitis",              "No",                       "Yes",                                               "Healthy periodontal status as baseline (PRIMARY predictor)",
  "smoking",                    "No",                       "Yes",                                               "Non-smoker as baseline (PRIMARY predictor)",
  "jaw",                        "Maxilla",                  "Mandible",                                          "Maxilla = typically lower bone density; mandible expected to perform better",
  "region",                     "Anterior",                 "Posterior",                                         "Anterior = thinner bone plate, different biomechanical loading",
  "implant_level",              "Bone Level",               "Tissue Level",                                      "Bone Level as modern baseline; PRIMARY variable of interest",
  "loading_protocol",           "Delayed",                  "Immediate",                                         "Delayed (conventional) = gold standard; PRIMARY predictor",
  "prosthesis_type",            "Crown (single tooth)",     "Bridge; Full-Arch",                                 "Simplest biomechanical unit as baseline",
  "surgical_protocol",          "Type 4 (delayed/healed)",  "Type 1; Type 2; Type 3",                           "Most conservative protocol = biological baseline",
  "implant_type",               "BL",                       "BLT; BLX; NN; RN; TE; TL; TLC; TLX; WN",         "BL = most common Bone Level standard",
  "implant_family",             "Bone-Level",               "Tissue-Level",                                      "2-level collapse: primary sensitivity model (recommended primary by v7)",
  "gbr",                        "No",                       "Yes",                                               "Non-augmented site as baseline",
  "scais",                      "No",                       "Yes",                                               "Standard implant as baseline",
  "flapless",                   "No",                       "Yes",                                               "Open-flap surgery as baseline (better visualisation)",
  "torque_group",               "<=20 Ncm",                 "21-35 Ncm; >=36 Ncm",                              "Low insertion torque as baseline",
  "multi_implant",              "Single",                   "Multiple",                                          "Single-implant patients as simpler anatomical baseline",
  "crown_implant_ratio_bin",    "<1",                       ">=1",                                               "CIR<1 considered biomechanically more favourable",
  "ea_mesial_bin",              "<10 deg",                  "[10,20); [20,30); >=30",                            "Low emergence angle = most permissive for oral hygiene",
  "ea_distal_bin3",             "<10 deg",                  "[10,20); >=20 (collapsed — n=3 in >=30)",          "Collapsed 3-level distal version to avoid near-empty cell",
  "ep_mesial",                  "Straight",                 "Concave; Convex",                                   "Straight profile = neutral emergence baseline",
  "ep_distal",                  "Straight",                 "Concave; Convex",                                   "Straight profile = neutral emergence baseline",
  "outer_contour_mesial_bin",   "<2 mm",                    "[2,4) mm; >=4 mm",                                  "Minimal bulk = least soft-tissue compression as baseline",
  "abutment_height_bin",        "<2 mm",                    ">=2 mm",                                            "Short abutment as bone-level-relevant baseline",
  "interprox_contact_mesial_bin","<2 mm",                   "[2,4) mm; >=4 mm",                                  "Short contact distance as baseline",
  "implant_to_tooth_mesial_bin","<2 mm",                    "[2,4) mm; >=4 mm",                                  "Minimum safe distance guideline (~1.5-2 mm)"
)

print(ref_table, n=Inf)
write.csv(ref_table, "outputs/Table_ReferenceCategoriesSupplement.csv", row.names=FALSE)
cat("  -> outputs/Table_ReferenceCategoriesSupplement.csv\n")
cat("  Include as supplementary table in thesis methods section.\n")


# ===========================================================================
# 29. FINAL SAVE + SESSION INFO
# ===========================================================================

cat("\n=== FINAL SAVE ===\n")

## Save all v7 tables
if (!is.null(pooled_abt_mi))
  write.csv(pooled_abt_mi %>% mutate(across(where(is.numeric), ~round(.x,4))),
            "outputs/MI/Table_MI_ModelB_abt.csv", row.names=FALSE)

write.csv(compare_3way %>% mutate(across(where(is.numeric), ~round(.x,4))),
          "outputs/Table_3way_LMM_GEE.csv", row.names=FALSE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                MBL ANALYSIS v7 — COMPLETE                      ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  outputs/Table1-9_*.csv            Core descriptive + model tables  ║\n")
cat("║  outputs/Fig1_MBL_over_time.png    Mean MBL trajectory              ║\n")
cat("║  outputs/Fig2_ImplantLevel_*.png   Bone vs Tissue Level (fixed)     ║\n")
cat("║  outputs/Fig3_diagnostics_*.png    Residual diagnostics              ║\n")
cat("║  outputs/Table_3way_LMM_GEE.csv   LMM vs GEE(exch) vs GEE(AR1)    ║\n")
cat("║  outputs/Table_InclExcl_*.csv      Selection bias checks             ║\n")
cat("║  outputs/Table_ReferenceCat*.csv   Reference category supplement     ║\n")
cat("║  outputs/MI/Table_MI_ModelB_abt.csv  MI pooled abutment model       ║\n")
cat("║  outputs/MI/Table_MNAR_TippingPoint.csv  MNAR sensitivity            ║\n")
cat("║  outputs/MI/Fig_MNAR_TippingPoint.png    Tipping-point figure        ║\n")
cat("║  outputs/MI/MI_convergence_abt_height.png  MI diagnostics            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\n✓ v7 complete. All outputs in 'outputs/' and 'outputs/MI/'.\n")
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_v7.txt")
