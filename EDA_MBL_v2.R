# =============================================================================
# EXPLORATORY DATA ANALYSIS (EDA)
# Script per capire i dati prima dell'analisi principale
#
# Da eseguire DOPO TesiGAb_updated.r (richiede df_A, df_B, raw_mbl in memoria)
# Tutti i grafici vengono salvati in outputs/EDA/
# =============================================================================
setwd("C:/Users/emilp/Downloads/Tesigab")
sink("console_EDA_MBL.txt", split = TRUE)
library(tidyverse)
library(patchwork)

dir.create("outputs/EDA", showWarnings = FALSE)

# Palette
COL_BL   <- "#2E75B6"
COL_TL   <- "#E74C3C"
COL_ACC  <- "#27AE60"
COL_WARN <- "#E67E22"

theme_eda <- theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        plot.caption  = element_text(size = 8,  color = "grey60"))

cat("╔══════════════════════════════════════════════════════╗\n")
cat("║         EXPLORATORY DATA ANALYSIS — MBL Study       ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")


# =============================================================================
# 1. PANORAMICA GENERALE
# =============================================================================

cat("── 1. PANORAMICA GENERALE ──────────────────────────────\n")
cat(sprintf("  Pazienti:         %d\n",  n_distinct(df_A$patient_id)))
cat(sprintf("  Impianti:         %d\n",  n_distinct(df_A$implant_id)))
cat(sprintf("  Osservazioni:     %d\n",  nrow(df_A)))
cat(sprintf("  Anni follow-up:   %d – %d\n", min(df_A$year), max(df_A$year)))
cat(sprintf("  MBL change range: %.3f – %.3f mm\n",
            min(df_A$mbl_change, na.rm = TRUE),
            max(df_A$mbl_change, na.rm = TRUE)))
cat(sprintf("  NA in mbl_change: %d (%.1f%%)\n",
            sum(is.na(df_A$mbl_change)),
            100 * mean(is.na(df_A$mbl_change))))


# =============================================================================
# 2. DISTRIBUZIONE VARIABILI PAZIENTE
# =============================================================================

cat("\n── 2. DISTRIBUZIONE VARIABILI PAZIENTE ────────────────\n")

df_pt <- df_A %>% distinct(patient_id, .keep_all = TRUE)

# 2a. Età
cat(sprintf("  Età — media: %.1f | mediana: %.1f | SD: %.1f | range: %d–%d\n",
            mean(df_pt$age_at_placement, na.rm = TRUE),
            median(df_pt$age_at_placement, na.rm = TRUE),
            sd(df_pt$age_at_placement, na.rm = TRUE),
            min(df_pt$age_at_placement, na.rm = TRUE),
            max(df_pt$age_at_placement, na.rm = TRUE)))

p_age <- ggplot(df_pt, aes(x = age_at_placement)) +
  geom_histogram(binwidth = 5, fill = COL_BL, color = "white", alpha = 0.85) +
  geom_vline(xintercept = median(df_pt$age_at_placement, na.rm = TRUE),
             linetype = "dashed", color = COL_WARN, linewidth = 0.8) +
  annotate("text",
           x = median(df_pt$age_at_placement, na.rm = TRUE) + 1,
           y = Inf, vjust = 1.5, hjust = 0, size = 3.2, color = COL_WARN,
           label = sprintf("Median = %.0f",
                           median(df_pt$age_at_placement, na.rm = TRUE))) +
  labs(title = "Age at implant placement",
       subtitle = "Histogram — binwidth 5 years",
       x = "Age (years)", y = "Count") +
  theme_eda

# 2b. Sesso
sex_tab <- df_pt %>% count(sex) %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat(sprintf("  Sesso — M: %d (%.1f%%) | F: %d (%.1f%%)\n",
            sex_tab$n[sex_tab$sex == "Male"],
            sex_tab$pct[sex_tab$sex == "Male"],
            sex_tab$n[sex_tab$sex == "Female"],
            sex_tab$pct[sex_tab$sex == "Female"]))

p_sex <- ggplot(sex_tab, aes(x = sex, y = n, fill = sex)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Male" = COL_BL, "Female" = COL_TL)) +
  scale_y_continuous(limits = c(0, max(sex_tab$n) * 1.2)) +
  labs(title = "Sex distribution", x = NULL, y = "Count") +
  theme_eda + theme(legend.position = "none")

# 2c. Comorbidità
comorbidity_df <- df_pt %>%
  summarise(
    Periodontitis   = sum(periodontitis   == "Yes", na.rm = TRUE),
    Smoking         = sum(smoking         == "Yes", na.rm = TRUE),
    Hypertension    = sum(hypertension    == "Yes", na.rm = TRUE),
    `Heart disease` = sum(heart_disease   == "Yes", na.rm = TRUE),
    Bisphosphonates = sum(bisphosphonates == "Yes", na.rm = TRUE),
    Diabetes        = sum(diabetes        == "Yes", na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "condition", values_to = "n") %>%
  mutate(
    pct       = round(100 * n / nrow(df_pt), 1),
    condition = fct_reorder(condition, n)
  )

cat("\n  Comorbidità:\n")
print(comorbidity_df %>% arrange(desc(n)))

p_comorbidity <- ggplot(comorbidity_df,
                        aes(x = condition, y = pct, fill = pct)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            hjust = -0.1, size = 3.2) +
  coord_flip() +
  scale_fill_gradient(low = "#D6E4F0", high = COL_BL) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(title = "Patient comorbidities (% of patients)",
       x = NULL, y = "% patients") +
  theme_eda + theme(legend.position = "none")

# Assembla e salva
p_patient <- (p_age | p_sex) / p_comorbidity +
  plot_annotation(title = "Patient-level characteristics (n = 200)")
ggsave("outputs/EDA/EDA1_Patient_Characteristics.png", p_patient,
       width = 12, height = 8, dpi = 300)
cat("  -> EDA1_Patient_Characteristics.png\n")


# =============================================================================
# 3. DISTRIBUZIONE VARIABILI IMPIANTO
# =============================================================================

cat("\n── 3. DISTRIBUZIONE VARIABILI IMPIANTO ────────────────\n")

df_imp <- df_A %>% distinct(implant_id, .keep_all = TRUE)

# 3a. Impianti per paziente
ipp <- df_A %>% distinct(patient_id, implant_id) %>%
  count(patient_id, name = "n_imp")
cat(sprintf("  Impianti/paz — mediana: %d | IQR: %d–%d | max: %d\n",
            as.integer(median(ipp$n_imp)),
            as.integer(quantile(ipp$n_imp, 0.25)),
            as.integer(quantile(ipp$n_imp, 0.75)),
            max(ipp$n_imp)))

p_ipp <- ggplot(ipp, aes(x = factor(n_imp))) +
  geom_bar(fill = COL_BL, alpha = 0.85, width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3.5) +
  labs(title = "Implants per patient",
       subtitle = sprintf("Median = %d  |  68%% of patients have 1 or 2 implants",
                          as.integer(median(ipp$n_imp))),
       x = "Number of implants", y = "Number of patients") +
  theme_eda

# 3b. Implant level
il_tab <- df_imp %>% count(implant_level) %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("\n  Implant level:\n"); print(il_tab)

p_il <- ggplot(il_tab, aes(x = implant_level, y = n, fill = implant_level)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Bone Level" = COL_BL, "Tissue Level" = COL_TL)) +
  scale_y_continuous(limits = c(0, max(il_tab$n) * 1.25)) +
  labs(title = "Implant level", x = NULL, y = "Count") +
  theme_eda + theme(legend.position = "none")

# 3c. Lunghezza e diametro
cat(sprintf("\n  Lunghezza — media: %.1f mm | SD: %.1f | range: %.1f–%.1f\n",
            mean(df_imp$length_mm, na.rm = TRUE),
            sd(df_imp$length_mm, na.rm = TRUE),
            min(df_imp$length_mm, na.rm = TRUE),
            max(df_imp$length_mm, na.rm = TRUE)))
cat(sprintf("  Diametro  — media: %.2f mm | SD: %.2f | range: %.1f–%.1f\n",
            mean(df_imp$diameter_mm, na.rm = TRUE),
            sd(df_imp$diameter_mm, na.rm = TRUE),
            min(df_imp$diameter_mm, na.rm = TRUE),
            max(df_imp$diameter_mm, na.rm = TRUE)))

p_length <- ggplot(df_imp, aes(x = factor(length_mm))) +
  geom_bar(fill = COL_ACC, alpha = 0.85, width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3) +
  labs(title = "Implant length", x = "Length (mm)", y = "Count") +
  theme_eda

p_diam <- ggplot(df_imp, aes(x = factor(round(diameter_mm, 1)))) +
  geom_bar(fill = COL_ACC, alpha = 0.85, width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3) +
  labs(title = "Implant diameter", x = "Diameter (mm)", y = "Count") +
  theme_eda

# 3d. Variabili categoriche impianto
cat_imp_df <- df_imp %>%
  summarise(
    Jaw_Mandible   = sum(jaw              == "Mandible",  na.rm = TRUE),
    Region_Post    = sum(region           == "Posterior", na.rm = TRUE),
    Loading_Imm    = sum(loading_protocol == "Immediate", na.rm = TRUE),
    GBR            = sum(gbr              == "Yes",       na.rm = TRUE),
    Flapless       = sum(flapless         == "Yes",       na.rm = TRUE),
    sCAIS          = sum(scais            == "Yes",       na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "var", values_to = "n") %>%
  mutate(pct = round(100 * n / nrow(df_imp), 1),
         var = fct_reorder(var, n))

p_cat_imp <- ggplot(cat_imp_df, aes(x = var, y = pct, fill = pct)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_gradient(low = "#D5F5E3", high = COL_ACC) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(title = "Implant-level categorical variables (% of implants)",
       x = NULL, y = "% implants") +
  theme_eda + theme(legend.position = "none")

# Assembla
p_implant <- (p_ipp | p_il) / (p_length | p_diam) / p_cat_imp +
  plot_annotation(title = "Implant-level characteristics (n = 501 implants)")
ggsave("outputs/EDA/EDA2_Implant_Characteristics.png", p_implant,
       width = 12, height = 14, dpi = 300)
cat("  -> EDA2_Implant_Characteristics.png\n")


# =============================================================================
# 4. DISTRIBUZIONE OUTCOME: MBL CHANGE
# =============================================================================

cat("\n── 4. DISTRIBUZIONE OUTCOME: MBL CHANGE ───────────────\n")

# Summary per anno
mbl_summary <- df_A %>%
  filter(!is.na(mbl_change)) %>%
  group_by(year) %>%
  summarise(
    n      = n(),
    mean   = round(mean(mbl_change), 3),
    median = round(median(mbl_change), 3),
    sd     = round(sd(mbl_change), 3),
    p25    = round(quantile(mbl_change, 0.25), 3),
    p75    = round(quantile(mbl_change, 0.75), 3),
    min    = round(min(mbl_change), 3),
    max    = round(max(mbl_change), 3),
    pct_above_mcid = round(100 * mean(mbl_change > 0.5), 1),
    .groups = "drop"
  )

cat("\n  Summary MBL change per anno:\n")
print(mbl_summary, n = 25)

# 4a. Distribuzione globale
p_mbl_hist <- ggplot(df_A %>% filter(!is.na(mbl_change), year > 0),
                     aes(x = mbl_change)) +
  geom_histogram(binwidth = 0.05, fill = COL_BL, color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.5,  linetype = "dashed", color = COL_WARN,
             linewidth = 0.8) +
  annotate("text", x = 0.52, y = Inf, vjust = 1.5, hjust = 0,
           size = 3, color = COL_WARN, label = "MCID (0.5 mm)") +
  labs(title = "Distribution of MBL change (all years, year 0 excluded)",
       x = "MBL change from baseline (mm)", y = "Count") +
  theme_eda

# 4b. Boxplot per anno (fino a Y15)
p_mbl_box <- df_A %>%
  filter(!is.na(mbl_change), year >= 1, year <= 15) %>%
  ggplot(aes(x = factor(year), y = mbl_change)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = COL_WARN, linewidth = 0.7) +
  geom_hline(yintercept = 0,   linetype = "solid",
             color = "grey60", linewidth = 0.4) +
  geom_boxplot(fill = COL_BL, alpha = 0.6, outlier.size = 0.8,
               outlier.alpha = 0.4, width = 0.6) +
  annotate("text", x = 0.6, y = 0.52, hjust = 0, vjust = 0,
           size = 2.8, color = COL_WARN, label = "MCID") +
  labs(title = "MBL change by follow-up year (Y1–Y15)",
       subtitle = "Boxplot: median, IQR, whiskers = 1.5xIQR",
       x = "Follow-up year", y = "MBL change (mm)") +
  theme_eda

# 4c. MBL change per implant level nel tempo
p_mbl_il <- df_A %>%
  filter(!is.na(mbl_change), year >= 1, year <= 15) %>%
  group_by(year, implant_level) %>%
  summarise(mean = mean(mbl_change), se = sd(mbl_change)/sqrt(n()),
            .groups = "drop") %>%
  ggplot(aes(x = year, y = mean, color = implant_level, fill = implant_level)) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = COL_WARN, linewidth = 0.7) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Bone Level" = COL_BL, "Tissue Level" = COL_TL)) +
  scale_fill_manual(values  = c("Bone Level" = COL_BL, "Tissue Level" = COL_TL)) +
  scale_x_continuous(breaks = 1:15) +
  labs(title = "Mean MBL change by implant level",
       subtitle = "Ribbon = +/- 1 SE (unadjusted, descriptive only)",
       x = "Follow-up year", y = "Mean MBL change (mm)",
       color = NULL, fill = NULL) +
  theme_eda + theme(legend.position = "bottom")

# 4d. % impianti > MCID per anno
p_pct_mcid <- mbl_summary %>%
  filter(year >= 1, year <= 15) %>%
  ggplot(aes(x = year, y = pct_above_mcid)) +
  geom_col(fill = COL_WARN, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(pct_above_mcid, "%")),
            vjust = -0.4, size = 3) +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(title = "% implants with MBL change > MCID (0.5 mm) per year",
       x = "Follow-up year", y = "% implants") +
  theme_eda

p_outcome <- (p_mbl_hist | p_mbl_box) / (p_mbl_il | p_pct_mcid) +
  plot_annotation(title = "Outcome: MBL change from baseline")
ggsave("outputs/EDA/EDA3_MBL_Distribution.png", p_outcome,
       width = 14, height = 10, dpi = 300)
cat("  -> EDA3_MBL_Distribution.png\n")


# =============================================================================
# 5. MISSING DATA
# =============================================================================

cat("\n── 5. MISSING DATA ─────────────────────────────────────\n")

miss_df <- df_A %>%
  select(all_of(c("mbl_change", covars_A))) %>%
  summarise(across(everything(), ~round(100 * mean(is.na(.x)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

cat("\n  Missingness per variabile (df_A):\n")
print(miss_df, n = Inf)

p_miss <- ggplot(miss_df %>% filter(pct_missing > 0),
                 aes(x = fct_reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = COL_TL, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(pct_missing, "%")), hjust = -0.2, size = 3.2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Missing data by variable (df_A)",
       subtitle = "Only variables with >0% missing shown",
       x = NULL, y = "% missing") +
  theme_eda

if (nrow(miss_df %>% filter(pct_missing > 0)) == 0) {
  cat("  Nessun valore mancante nelle variabili principali.\n")
} else {
  ggsave("outputs/EDA/EDA4_Missing_Data.png", p_miss,
         width = 8, height = 5, dpi = 300)
  cat("  -> EDA4_Missing_Data.png\n")
}

# Missing in df_B (variabili protesiche)
miss_B <- df_B %>%
  select(crown_implant_ratio, ea_mesial_deg, abutment_height_mm,
         outer_contour_mesial_mm, interprox_contact_mesial_mm,
         implant_to_tooth_mesial_mm, interimplant_dist_mesial_mm,
         ep_mesial) %>%
  summarise(across(everything(), ~round(100 * mean(is.na(.x)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

cat("\n  Missingness variabili protesiche (df_B):\n")
print(miss_B, n = Inf)


# =============================================================================
# 6. FOLLOW-UP: NUMERO IMPIANTI PER ANNO
# =============================================================================

cat("\n── 6. FOLLOW-UP ────────────────────────────────────────\n")

fu_tab <- df_A %>%
  group_by(year) %>%
  summarise(n_implants = n_distinct(implant_id),
            n_patients = n_distinct(patient_id),
            .groups = "drop")

cat("\n  Impianti/pazienti per anno:\n")
print(fu_tab, n = 25)

p_fu <- ggplot(fu_tab, aes(x = year)) +
  geom_col(aes(y = n_implants), fill = COL_BL, alpha = 0.8, width = 0.7) +
  geom_line(aes(y = n_patients * (max(n_implants) / max(n_patients))),
            color = COL_TL, linewidth = 1.2) +
  geom_point(aes(y = n_patients * (max(n_implants) / max(n_patients))),
             color = COL_TL, size = 2.5) +
  scale_y_continuous(
    name = "N implants",
    sec.axis = sec_axis(
      ~ . * max(fu_tab$n_patients) / max(fu_tab$n_implants),
      name = "N patients"
    )
  ) +
  scale_x_continuous(breaks = 0:20) +
  labs(title = "Number of implants and patients per year",
       subtitle = "Blue bars = implants | Red line = patients",
       x = "Follow-up year") +
  theme_eda

ggsave("outputs/EDA/EDA5_FollowUp.png", p_fu,
       width = 10, height = 5, dpi = 300)
cat("  -> EDA5_FollowUp.png\n")


# =============================================================================
# 7. CORRELAZIONI TRA VARIABILI CONTINUE
# =============================================================================

cat("\n── 7. CORRELAZIONI VARIABILI CONTINUE ─────────────────\n")

cont_vars <- df_A %>%
  distinct(implant_id, .keep_all = TRUE) %>%
  select(age_at_placement, length_mm, diameter_mm, insertion_torque_ncm) %>%
  drop_na()

cor_mat <- round(cor(cont_vars, use = "complete.obs"), 2)
cat("\n  Matrice di correlazione (Pearson):\n")
print(cor_mat)

# Scatter matrix semplice
pairs_df <- cont_vars %>%
  pivot_longer(everything(), names_to = "var1", values_to = "val1") %>%
  left_join(
    cont_vars %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "var2", values_to = "val2"),
    by = character(),
    relationship = "many-to-many"
  ) %>%
  filter(var1 != var2)

# Scatterplot pairwise delle variabili principali
p_scatter_pairs <- GGally::ggpairs(
  cont_vars,
  upper = list(continuous = GGally::wrap("cor", size = 3.5)),
  lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = 0.8,
                                          color = COL_BL)),
  diag  = list(continuous = GGally::wrap("densityDiag", fill = COL_BL,
                                          alpha = 0.5))
) +
  labs(title = "Pairwise correlations — continuous implant variables") +
  theme_minimal(base_size = 10)

# Installa GGally se mancante
if (!requireNamespace("GGally", quietly = TRUE)) {
  install.packages("GGally")
  library(GGally)
  ggsave("outputs/EDA/EDA6_Correlations.png", p_scatter_pairs,
         width = 9, height = 9, dpi = 300)
  cat("  -> EDA6_Correlations.png\n")
} else {
  library(GGally)
  ggsave("outputs/EDA/EDA6_Correlations.png", p_scatter_pairs,
         width = 9, height = 9, dpi = 300)
  cat("  -> EDA6_Correlations.png\n")
}


# =============================================================================
# 8. MBL CHANGE PER SOTTOGRUPPI CLINICI
# =============================================================================

cat("\n── 8. MBL CHANGE PER SOTTOGRUPPI ──────────────────────\n")

# Funzione boxplot per una variabile categorica
plot_mbl_by <- function(var, label, colors = NULL) {
  d <- df_A %>%
    filter(!is.na(mbl_change), !is.na(.data[[var]]), year >= 1, year <= 10)

  if (is.null(colors)) colors <- scales::hue_pal()(n_distinct(d[[var]]))

  ggplot(d, aes(x = .data[[var]], y = mbl_change, fill = .data[[var]])) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = COL_WARN, linewidth = 0.6) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.6, outlier.alpha = 0.3,
                 width = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 18,
                 size = 3, color = "black") +
    scale_fill_manual(values = colors) +
    labs(title = label,
         subtitle = "Y1-Y10 | Diamond = mean | Dashed = MCID",
         x = NULL, y = "MBL change (mm)") +
    theme_eda + theme(legend.position = "none")
}

p_by_il    <- plot_mbl_by("implant_level", "By implant level",
                           c("Bone Level" = COL_BL, "Tissue Level" = COL_TL))
p_by_jaw   <- plot_mbl_by("jaw", "By jaw",
                           c("Maxilla" = "#8E44AD", "Mandible" = "#F39C12"))
p_by_load  <- plot_mbl_by("loading_protocol", "By loading protocol",
                           c("Delayed" = COL_BL, "Immediate" = COL_TL))
p_by_perio <- plot_mbl_by("periodontitis", "By periodontitis",
                           c("No" = COL_ACC, "Yes" = COL_TL))
p_by_smoke <- plot_mbl_by("smoking", "By smoking",
                           c("No" = COL_ACC, "Yes" = COL_TL))
p_by_gbr   <- plot_mbl_by("gbr", "By GBR",
                           c("No" = COL_BL, "Yes" = "#F39C12"))

p_subgroups <- (p_by_il | p_by_jaw | p_by_load) /
               (p_by_perio | p_by_smoke | p_by_gbr) +
  plot_annotation(title = "MBL change by clinical subgroup (Y1-Y10, unadjusted)",
                  caption = "Note: these are UNADJUSTED comparisons — use Model A for inference")
ggsave("outputs/EDA/EDA7_MBL_by_Subgroup.png", p_subgroups,
       width = 14, height = 9, dpi = 300)
cat("  -> EDA7_MBL_by_Subgroup.png\n")


# =============================================================================
# 9. VARIABILI PROTESICHE (df_B)
# =============================================================================

cat("\n── 9. VARIABILI PROTESICHE (df_B) ──────────────────────\n")

df_B_imp <- df_B %>% distinct(implant_id, .keep_all = TRUE)

prosth_cont <- df_B_imp %>%
  select(crown_implant_ratio, ea_mesial_deg, ea_distal_deg,
         abutment_height_mm, outer_contour_mesial_mm,
         interprox_contact_mesial_mm, implant_to_tooth_mesial_mm) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

cat("\n  Summary variabili protesiche:\n")
prosth_cont %>%
  group_by(variable) %>%
  summarise(
    n      = n(),
    mean   = round(mean(value), 2),
    median = round(median(value), 2),
    sd     = round(sd(value), 2),
    min    = round(min(value), 2),
    max    = round(max(value), 2),
    .groups = "drop"
  ) %>% print(n = Inf)

p_prosth <- ggplot(prosth_cont, aes(x = value)) +
  geom_histogram(bins = 25, fill = "#8E44AD", color = "white", alpha = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  labs(title = "Distribution of prosthetic variables",
       x = NULL, y = "Count") +
  theme_eda + theme(strip.text = element_text(size = 8))

ggsave("outputs/EDA/EDA8_Prosthetic_Variables.png", p_prosth,
       width = 12, height = 8, dpi = 300)
cat("  -> EDA8_Prosthetic_Variables.png\n")

# EA high distribution
ea_tab <- df_B_imp %>%
  mutate(ea_high = factor(
    case_when(
      is.na(ea_mesial_bin)    ~ NA_character_,
      ea_mesial_bin == ">=30" ~ "High (>30°)",
      TRUE                    ~ "Normal (<=30°)"
    ),
    levels = c("Normal (<=30°)", "High (>30°)")
  )) %>%
  count(ea_high) %>%
  mutate(pct = round(100 * n / sum(n, na.rm = TRUE), 1))
cat("\n  Emergence angle high (>30 deg):\n"); print(ea_tab)


# =============================================================================
# RIEPILOGO FINALE
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║                  EDA COMPLETATA                     ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")

eda_files <- list.files("outputs/EDA", full.names = FALSE)
for (f in eda_files) {
  kb <- round(file.size(file.path("outputs/EDA", f)) / 1024, 0)
  cat(sprintf("║  %-42s %4d KB ║\n", f, kb))
}
cat("╚══════════════════════════════════════════════════════╝\n")
cat("\nTutti i file sono in outputs/EDA/\n")
sink()