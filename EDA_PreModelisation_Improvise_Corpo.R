# ==============================================================================
# EDA PRÉ-MODÉLISATION : DONNÉES CORPORATE
# Objectif : Orienter le choix entre modèle Ridge-pénalisé, hiérarchique complet
#            (nlme/Stan), ou Processus Gaussien pour la courbe de spread.
#
# Structure :
#   Section 1  — Sparsité : combien d'obligations par émetteur / date ?
#   Section 2  — Homogénéité des groupes : rating × secteur sont-ils cohérents ?
#   Section 3  — Structure empirique du spread : quelle forme fonctionnelle ?
#   Section 4  — Variance inter- vs intra-émetteur : justification du shrinkage
#   Section 5  — Stabilité temporelle : les spreads bougent-ils beaucoup ?
#   Section 6  — Corrélations inter-paramètres NSS : faut-il fixer τ ?
#   Section 7  — Synthèse : tableau de décision pour le choix de méthode
# ==============================================================================

rm(list = ls())

# ── PACKAGES ──────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(scales)
library(ggrepel)
source("utils.R")   # nss_func(), calculate_duration(), filter_ecb_criteria()

set.seed(42)

# ── PALETTE ───────────────────────────────────────────────────────────────────
rating_colors <- c(
  "AA"  = "#1B4F8A",
  "A"   = "#2E86C1",
  "BBB" = "#F39C12",
  "HY"  = "#C0392B",
  "NR"  = "#7F8C8D"
)

sector_colors <- c(
  "Communications"          = "#8E44AD",
  "Consumer Discretionary"  = "#E74C3C",
  "Consumer Staples"        = "#E67E22",
  "Energy"                  = "#F1C40F",
  "Financials"              = "#2ECC71",
  "Health Care"             = "#1ABC9C",
  "Industrials"             = "#3498DB",
  "Materials"               = "#95A5A6",
  "Real Estate"             = "#D35400",
  "Technology"              = "#9B59B6",
  "Utilities"               = "#16A085"
)

# ==============================================================================
# 0. CHARGEMENT ET PRÉ-TRAITEMENT
# ==============================================================================

corpo_raw <- read_csv("r_data/NewData_Corpo_clean.csv", show_col_types = FALSE) |>
  mutate(across(starts_with("YLD_"), ~ suppressWarnings(as.numeric(gsub(",", ".", .x)))))

yld_cols <- grep("^YLD_", names(corpo_raw), value = TRUE)

# Format long  — 1 ligne = 1 obligation × 1 date
corpo_long <- corpo_raw |>
  pivot_longer(cols = all_of(yld_cols), names_to = "Date_Label", values_to = "Yield") |>
  mutate(
    Date          = dmy(gsub("^YLD_", "", Date_Label)),
    Maturity      = as.numeric(Residual_Maturity),
    Amt_Out       = suppressWarnings(as.numeric(`Amt Out`)),
    Cpn           = suppressWarnings(as.numeric(Cpn)),
    Yield         = as.numeric(Yield)
  ) |>
  filter(!is.na(Yield), !is.na(Maturity), Maturity > 0.25, Maturity <= 30) |>
  rename(
    Rating    = `BBG Composite`,
    Sector_L1 = `BICS Level 1`,
    Sector_L2 = `BICS Level 2`,
    Country   = `Cntry of Risk`
  ) |>
  mutate(
    Rating_Bucket = case_when(
      Rating %in% c("AAA", "AA+", "AA", "AA-")       ~ "AA",
      Rating %in% c("A+", "A", "A-")                  ~ "A",
      Rating %in% c("BBB+", "BBB", "BBB-")             ~ "BBB",
      Rating %in% c("BB+", "BB", "BB-",
                    "B+", "B", "B-", "CCC", "CC", "C") ~ "HY",
      TRUE                                              ~ "NR"
    )
  )

# Courbe swap risk-free (snapshot statique)
swap_raw <- read_csv("r_data/StatApp_Data1_Swap.csv", show_col_types = FALSE)

clean_swap <- function(df) {
  cvt <- function(t) {
    if (grepl("M", t)) return(as.numeric(gsub("M","",t))/12)
    if (grepl("Y", t)) return(as.numeric(gsub("Y","",t)))
    NA
  }
  out <- df |>
    select(TenorStr = 1, SwapRate = 3) |>
    mutate(Maturity = sapply(TenorStr, cvt),
           Yield    = suppressWarnings(as.numeric(as.character(SwapRate)))) |>
    filter(!is.na(Maturity), !is.na(Yield))
  if (mean(out$Yield, na.rm=TRUE) < 0.5) out$Yield <- out$Yield * 100
  out
}

df_swap <- clean_swap(swap_raw)

# Calibration Svensson sur le swap
p_swap_init <- c(b0=4, b1=-0.5, b2=10, b3=-5, tau1=2.5, tau2=0.7)
swap_opt <- optim(
  par = p_swap_init,
  fn  = nss_objective,
  maturities = df_swap$Maturity,
  yields     = df_swap$Yield,
  weights    = rep(1, nrow(df_swap)),
  method = "L-BFGS-B",
  lower  = c(0, -15, -15, -15, 0.1, 0.1),
  upper  = c(20, 15,  15,  15, 15,  30),
  control = list(maxit = 5000, factr = 1e6)
)
p_swap <- swap_opt$par

# Calcul du spread observé pour chaque obligation × date
corpo_long <- corpo_long |>
  mutate(
    Yield_RF  = nss_func(Maturity, p_swap),
    Spread    = Yield - Yield_RF
  )

message("Dataset : ", nrow(corpo_long), " obs | ",
        length(unique(corpo_long$ISIN)), " ISINs | ",
        length(unique(corpo_long$Ticker)), " tickers | ",
        length(unique(corpo_long$Date)), " dates")

# ==============================================================================
# SECTION 1 — SPARSITÉ : DISTRIBUTION DU NOMBRE D'OBLIGATIONS PAR ÉMETTEUR
# ==============================================================================
# Question : À quel point le problème de sparsité est-il répandu ?
# Décision : Si la majorité des émetteurs a < 6 bonds, le shrinkage est
#            indispensable. Si beaucoup en ont > 10, la Ridge seule suffit.

# ── 1A. Distribution par ticker (toutes dates confondues, snapshot dernière date)
last_date <- max(corpo_long$Date)

sparsity_per_ticker <- corpo_long |>
  filter(Date == last_date) |>
  group_by(Ticker, Rating_Bucket, Sector_L1) |>
  summarise(N_bonds = n_distinct(ISIN), .groups = "drop") |>
  arrange(desc(N_bonds))

p1a <- ggplot(sparsity_per_ticker, aes(x = N_bonds, fill = Rating_Bucket)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 6,  linetype = "dashed", color = "#E74C3C", linewidth = 1) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#E67E22", linewidth = 0.8) +
  annotate("text", x = 6.3,  y = Inf, vjust = 1.5, hjust = 0,
           label = "≥ 6 : NSS identifiable", color = "#E74C3C", size = 3.2, fontface = "italic") +
  annotate("text", x = 10.3, y = Inf, vjust = 3.5, hjust = 0,
           label = "≥ 10 : NSS robuste", color = "#E67E22", size = 3.2, fontface = "italic") +
  scale_fill_manual(values = rating_colors, name = "Rating") +
  scale_x_continuous(breaks = 0:max(sparsity_per_ticker$N_bonds)) +
  labs(
    title    = "Section 1A — Distribution du nombre d'obligations par émetteur",
    subtitle = paste0("Snapshot au ", format(last_date, "%d/%m/%Y"),
                      " | N émetteurs = ", nrow(sparsity_per_ticker)),
    x = "Nombre d'obligations", y = "Nombre d'émetteurs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
print(p1a)

# ── 1B. Tableau récapitulatif : catégorisation des émetteurs
sparsity_summary <- sparsity_per_ticker |>
  mutate(
    Sparsity_Class = case_when(
      N_bonds < 3  ~ "Très sparse (< 3)",
      N_bonds < 6  ~ "Sparse (3–5)",
      N_bonds < 10 ~ "Modéré (6–9)",
      TRUE         ~ "Riche (≥ 10)"
    )
  ) |>
  count(Sparsity_Class, name = "N_emetteurs") |>
  mutate(Pct = round(100 * N_emetteurs / sum(N_emetteurs), 1)) |>
  arrange(match(Sparsity_Class, c("Très sparse (< 3)", "Sparse (3–5)",
                                  "Modéré (6–9)", "Riche (≥ 10)")))

message("\n── Section 1 : Sparsité des émetteurs ──")
print(sparsity_summary)

# ── 1C. Top 10 vs Bottom 10 émetteurs (pour savoir qui "domine" le panel)
message("\nTop 10 émetteurs (plus de bonds) :")
print(head(sparsity_per_ticker, 10))
message("Bottom 10 émetteurs :")
print(tail(sparsity_per_ticker, 10))

# ==============================================================================
# SECTION 2 — HOMOGÉNÉITÉ DES GROUPES rating × secteur
# ==============================================================================
# Question : Les groupes rating × secteur sont-ils suffisamment cohérents
#            pour fournir un a priori informatif ?
# Décision : Si spread moyen par groupe varie beaucoup ENTRE groupes mais
#            peu DANS un groupe → groupes pertinents → shrinkage utile.

# ── 2A. Grille rating × secteur : effectifs et spread moyen
group_stats <- corpo_long |>
  filter(Date == last_date, !is.na(Spread)) |>
  group_by(Rating_Bucket, Sector_L1) |>
  summarise(
    N_bonds      = n(),
    N_tickers    = n_distinct(Ticker),
    Spread_Mean  = round(mean(Spread, na.rm = TRUE), 3),
    Spread_SD    = round(sd(Spread,   na.rm = TRUE), 3),
    CV           = round(Spread_SD / abs(Spread_Mean), 2),  # Coefficient de variation
    .groups = "drop"
  ) |>
  arrange(Rating_Bucket, Sector_L1)

message("\n── Section 2 : Homogénéité des groupes rating × secteur ──")
print(as.data.frame(group_stats))

# ── 2B. Heatmap effectifs
p2b <- group_stats |>
  ggplot(aes(x = Rating_Bucket, y = Sector_L1, fill = N_bonds)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(N_bonds, "\n(", N_tickers, " tic.)")),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#BDC3C7", high = "#1B4F8A",
                      name = "N obligations", na.value = "grey95") +
  labs(
    title    = "Section 2B — Effectifs par groupe rating × secteur",
    subtitle = "Format cellule : N obligations (N tickers)",
    x = "Rating Bucket", y = "Secteur (BICS L1)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))
print(p2b)

# ── 2C. Boxplot spread par groupe (vérifie l'homogénéité intra-groupe)
p2c <- corpo_long |>
  filter(Date == last_date, !is.na(Spread),
         abs(Spread) < 10) |>  # Retire outliers extrêmes pour la visu
  mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA","A","BBB","HY","NR"))) |>
  ggplot(aes(x = Sector_L1, y = Spread, fill = Rating_Bucket)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1, outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = rating_colors, name = "Rating") +
  labs(
    title    = "Section 2C — Distribution du spread par groupe",
    subtitle = "Homogénéité intra-groupe = base solide pour l'a priori",
    x = NULL, y = "Spread vs Swap (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8)
  )
print(p2c)

# ==============================================================================
# SECTION 3 — STRUCTURE EMPIRIQUE DU SPREAD : QUELLE FORME FONCTIONNELLE ?
# ==============================================================================
# Question : Le spread est-il plat, monotone, en cloche ? NSS est-il sur-dimensionné ?
# Décision : Si le spread est globalement monotone/plat, NS3 (4 params) suffit.
#            Si on voit des bosses, NSS6 se justifie.

# ── 3A. Spread vs maturité par rating bucket (courbe LOESS empirique)
p3a <- corpo_long |>
  filter(Date == last_date, !is.na(Spread), abs(Spread) < 8) |>
  mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA","A","BBB","HY"))) |>
  ggplot(aes(x = Maturity, y = Spread, color = Rating_Bucket)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, span = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = rating_colors, name = "Rating") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(
    title    = "Section 3A — Structure empirique du spread par rating",
    subtitle = "Courbe LOESS — indique la forme fonctionnelle naturelle",
    x = "Maturité (Années)", y = "Spread vs Swap (%)"
  ) +
  facet_wrap(~Rating_Bucket, scales = "free_y") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
print(p3a)

# ── 3B. Spread par secteur (tous ratings confondus, maturité bucketée)
p3b <- corpo_long |>
  filter(Date == last_date, !is.na(Spread), abs(Spread) < 8) |>
  mutate(Mat_Bucket = cut(Maturity,
                          breaks = c(0, 2, 5, 10, 30),
                          labels = c("0–2Y", "2–5Y", "5–10Y", "10–30Y"))) |>
  group_by(Sector_L1, Mat_Bucket) |>
  summarise(Spread_Mean = mean(Spread, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = Mat_Bucket, y = Spread_Mean, group = Sector_L1, color = Sector_L1)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = sector_colors, name = "Secteur") +
  labs(
    title    = "Section 3B — Profil de spread moyen par secteur et maturité",
    subtitle = "Révèle si la forme de la term structure varie selon le secteur",
    x = "Bucket de maturité", y = "Spread moyen (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
print(p3b)

# ── 3C. Test : NS3 (3 params) vs NSS6 — pour les émetteurs riches
# On estime les deux sur les tickers avec >= 8 bonds et on compare les RMSE

rich_tickers <- sparsity_per_ticker |>
  filter(N_bonds >= 8) |>
  pull(Ticker)

fit_ns3 <- function(df) {
  # NS3 = Nelson-Siegel pur (sans 4ème terme) : 4 paramètres (b0,b1,b2,tau1)
  cost <- function(params) {
    if (params[4] <= 0) return(1e9)
    tau   <- params[4]
    x1    <- (1 - exp(-df$Maturity/tau)) / (df$Maturity/tau)
    x2    <- x1 - exp(-df$Maturity/tau)
    fitted <- params[1] + params[2]*x1 + params[3]*x2
    sqrt(mean((df$Spread - fitted)^2))
  }
  opt <- optim(c(mean(df$Spread), -0.1, 0, 2), cost, method = "Nelder-Mead",
               control = list(maxit = 3000))
  list(rmse = opt$value, params = opt$par, n = nrow(df))
}

fit_nss6 <- function(df) {
  cost <- function(params) {
    if (params[5] <= 0 || params[6] <= 0) return(1e9)
    fitted <- nss_func(df$Maturity, params)
    sqrt(mean((df$Spread - fitted)^2))
  }
  opt <- optim(c(mean(df$Spread), -0.1, 0, 0, 2, 8), cost,
               method = "Nelder-Mead", control = list(maxit = 5000))
  list(rmse = opt$value, params = opt$par, n = nrow(df))
}

ns3_vs_nss6 <- corpo_long |>
  filter(Date == last_date, Ticker %in% rich_tickers, !is.na(Spread)) |>
  group_by(Ticker, Rating_Bucket) |>
  do({
    df_sub <- .
    r3  <- tryCatch(fit_ns3(df_sub),  error = function(e) list(rmse = NA, n = nrow(df_sub)))
    r6  <- tryCatch(fit_nss6(df_sub), error = function(e) list(rmse = NA, n = nrow(df_sub)))
    data.frame(
      N_bonds  = nrow(df_sub),
      RMSE_NS3 = r3$rmse,
      RMSE_NSS = r6$rmse,
      Delta_RMSE = r3$rmse - r6$rmse  # > 0 : NSS6 gagne
    )
  }) |>
  ungroup() |>
  filter(!is.na(RMSE_NS3), !is.na(RMSE_NSS))

message("\n── Section 3C : NS3 vs NSS6 sur émetteurs riches (RMSE sur spread) ──")
print(ns3_vs_nss6)

p3c <- ns3_vs_nss6 |>
  pivot_longer(cols = c(RMSE_NS3, RMSE_NSS), names_to = "Model", values_to = "RMSE") |>
  ggplot(aes(x = reorder(Ticker, -RMSE), y = RMSE, fill = Model)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("RMSE_NS3" = "#2E86C1", "RMSE_NSS" = "#E74C3C"),
                    labels = c("Nelson-Siegel (4 params)", "Svensson NSS (6 params)"),
                    name = "Modèle") +
  labs(
    title    = "Section 3C — NS3 vs NSS6 sur les émetteurs riches en données",
    subtitle = "Si la différence est faible, NS3 suffit → moins de paramètres à shrinkage",
    x = "Émetteur", y = "RMSE sur le spread (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )
print(p3c)

# ==============================================================================
# SECTION 4 — VARIANCE INTER- vs INTRA-ÉMETTEUR
# ==============================================================================
# Question : Quelle part de la variance du spread est expliquée par l'identité
#            de l'émetteur vs le bruit de marché intra-émetteur ?
# Décision : Si variance inter >> variance intra → groupes pertinents,
#            le shrinkage vers μ_group est très utile.
#            Si variance intra >> variance inter → le bruit domine,
#            l'a priori aidera peu et il faut un lambda de pénalité fort.

# ANOVA à un facteur sur le spread (ticker comme facteur)
df_anova <- corpo_long |>
  filter(Date == last_date, !is.na(Spread), abs(Spread) < 10) |>
  filter(Ticker %in% sparsity_per_ticker$Ticker[sparsity_per_ticker$N_bonds >= 3])

aov_ticker <- aov(Spread ~ Ticker, data = df_anova)
aov_summary <- summary(aov_ticker)[[1]]
SS_between  <- aov_summary["Ticker", "Sum Sq"]
SS_within   <- aov_summary["Residuals", "Sum Sq"]
eta_squared <- SS_between / (SS_between + SS_within)

message("\n── Section 4 : Décomposition de la variance du spread ──")
message("SS Between (inter-émetteur) : ", round(SS_between, 4))
message("SS Within  (intra-émetteur) : ", round(SS_within,  4))
message("η² (part inter) : ", round(eta_squared, 4),
        " — interprétation : >0.5 = shrinkage très utile")

# ── 4B. Idem par maturité bucket (pour voir si la structure varie)
df_var_decomp <- corpo_long |>
  filter(Date == last_date, !is.na(Spread), abs(Spread) < 10) |>
  mutate(Mat_Bucket = cut(Maturity, breaks = c(0, 2, 5, 10, 30),
                          labels = c("0–2Y", "2–5Y", "5–10Y", "10–30Y"))) |>
  group_by(Mat_Bucket) |>
  summarise(
    Var_Total  = var(Spread, na.rm = TRUE),
    Var_Within = mean(tapply(Spread, Ticker, var, na.rm = TRUE), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    Var_Between  = Var_Total - Var_Within,
    Pct_Between  = round(100 * Var_Between / Var_Total, 1),
    Pct_Within   = round(100 * Var_Within  / Var_Total, 1)
  )

message("\n── Décomposition par bucket de maturité ──")
print(df_var_decomp)

p4 <- df_var_decomp |>
  select(Mat_Bucket, Pct_Between, Pct_Within) |>
  pivot_longer(cols = c(Pct_Between, Pct_Within), names_to = "Source", values_to = "Pct") |>
  mutate(Source = recode(Source,
                         "Pct_Between" = "Inter-émetteur (signal)",
                         "Pct_Within"  = "Intra-émetteur (bruit)")) |>
  ggplot(aes(x = Mat_Bucket, y = Pct, fill = Source)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("Inter-émetteur (signal)" = "#1B4F8A",
                               "Intra-émetteur (bruit)"  = "#AED6F1"),
                    name = NULL) +
  geom_text(aes(label = paste0(Pct, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  labs(
    title    = "Section 4 — Décomposition variance spread : inter vs intra-émetteur",
    subtitle = "Part bleue foncée = signal que le modèle peut capturer",
    x = "Bucket de maturité", y = "Part de variance (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
print(p4)

# ==============================================================================
# SECTION 5 — STABILITÉ TEMPORELLE DES SPREADS
# ==============================================================================
# Question : Les spreads sont-ils stables d'une date à l'autre (faible volatilité
#            temporelle) ou très mobiles ?
# Décision : Faible vol → on peut utiliser le spread cross-sectionnel comme a priori.
#            Forte vol → il faut capturer la dynamique temporelle (AR1 ou état-espace).

# ── 5A. Volatilité du spread par émetteur dans le temps
spread_vol_time <- corpo_long |>
  filter(!is.na(Spread), abs(Spread) < 10) |>
  group_by(Ticker, Rating_Bucket) |>
  summarise(
    Spread_Mean  = mean(Spread,   na.rm = TRUE),
    Spread_SD    = sd(Spread,     na.rm = TRUE),
    N_dates      = n_distinct(Date),
    CV_Time      = round(Spread_SD / abs(Spread_Mean), 3),
    .groups = "drop"
  ) |>
  filter(N_dates >= 4)  # Suffisamment de dates pour calculer une vol

message("\n── Section 5 : Volatilité temporelle des spreads ──")
message("Médiane CV temporel : ", round(median(spread_vol_time$CV_Time, na.rm = TRUE), 3))
message("Ratio: si CV < 0.10 → stable, si > 0.30 → très mobile")

p5a <- ggplot(spread_vol_time, aes(x = CV_Time, fill = Rating_Bucket)) +
  geom_histogram(binwidth = 0.02, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "#27AE60", linewidth = 0.8) +
  geom_vline(xintercept = 0.30, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 0.11, y = Inf, vjust = 1.5, hjust = 0,
           label = "Stable", color = "#27AE60", size = 3.2, fontface = "italic") +
  annotate("text", x = 0.31, y = Inf, vjust = 1.5, hjust = 0,
           label = "Mobile", color = "#E74C3C", size = 3.2, fontface = "italic") +
  scale_fill_manual(values = rating_colors, name = "Rating") +
  labs(
    title    = "Section 5A — Coefficient de variation temporel du spread par émetteur",
    subtitle = "CV = SD_temporelle / |Spread_moyen|",
    x = "Coefficient de variation (CV)", y = "Nombre d'émetteurs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
print(p5a)

# ── 5B. Evolution temporelle du spread moyen par rating bucket
p5b <- corpo_long |>
  filter(!is.na(Spread), abs(Spread) < 10) |>
  group_by(Date, Rating_Bucket) |>
  summarise(Spread_Mean = mean(Spread, na.rm = TRUE),
            Spread_SD   = sd(Spread,   na.rm = TRUE),
            .groups = "drop") |>
  mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA","A","BBB","HY"))) |>
  ggplot(aes(x = Date, y = Spread_Mean, color = Rating_Bucket, fill = Rating_Bucket)) +
  geom_ribbon(aes(ymin = Spread_Mean - Spread_SD,
                  ymax = Spread_Mean + Spread_SD), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = rating_colors, name = "Rating") +
  scale_fill_manual(values = rating_colors,  name = "Rating") +
  scale_x_date(date_labels = "%d/%m", date_breaks = "2 weeks") +
  labs(
    title    = "Section 5B — Évolution temporelle du spread moyen par rating",
    subtitle = "Ruban = ± 1σ | Indique la dynamique à modéliser",
    x = NULL, y = "Spread moyen (%)"
  ) +
  theme_minimal() +
  theme(plot.title   = element_text(face = "bold"),
        axis.text.x  = element_text(angle = 45, hjust = 1))
print(p5b)

# ==============================================================================
# SECTION 6 — CORRÉLATIONS INTER-PARAMÈTRES NSS ET QUESTION DU τ
# ==============================================================================
# Question : Faut-il fixer τ1 et τ2 globalement (comme dans Diebold-Li)
#            ou les laisser libres par émetteur ?
# Décision : Si τ est très variable entre émetteurs → le laisser libre.
#            Si τ est stable → le fixer pour réduire la dimensionnalité.

# On estime un NSS sur chaque ticker riche et on regarde la distribution de τ
nss_params_rich <- corpo_long |>
  filter(Date == last_date, Ticker %in% rich_tickers, !is.na(Spread)) |>
  group_by(Ticker, Rating_Bucket, Sector_L1) |>
  do({
    df <- .
    tryCatch({
      cost <- function(p) {
        if (p[5] <= 0 || p[6] <= 0) return(1e9)
        fitted <- nss_func(df$Maturity, p)
        sqrt(mean((df$Spread - fitted)^2))
      }
      opt <- optim(c(mean(df$Spread), -0.1, 0, 0, 2, 8), cost,
                   method = "Nelder-Mead", control = list(maxit = 5000))
      p <- opt$par
      data.frame(b0=p[1], b1=p[2], b2=p[3], b3=p[4],
                 tau1=p[5], tau2=p[6], rmse=opt$value, n=nrow(df))
    }, error = function(e) data.frame(b0=NA,b1=NA,b2=NA,b3=NA,
                                      tau1=NA,tau2=NA,rmse=NA,n=nrow(df)))
  }) |>
  ungroup() |>
  filter(!is.na(tau1))

message("\n── Section 6 : Distribution des τ estimés sur émetteurs riches ──")
message("τ1 : médiane=", round(median(nss_params_rich$tau1), 2),
        " | CV=", round(sd(nss_params_rich$tau1)/mean(nss_params_rich$tau1), 2))
message("τ2 : médiane=", round(median(nss_params_rich$tau2), 2),
        " | CV=", round(sd(nss_params_rich$tau2)/mean(nss_params_rich$tau2), 2))

p6a <- nss_params_rich |>
  select(Ticker, tau1, tau2) |>
  pivot_longer(cols = c(tau1, tau2), names_to = "Param", values_to = "Value") |>
  ggplot(aes(x = Value, fill = Param)) +
  geom_histogram(binwidth = 0.5, color = "white", alpha = 0.8, position = "identity") +
  facet_wrap(~Param, scales = "free_x") +
  scale_fill_manual(values = c("tau1" = "#2E86C1", "tau2" = "#E74C3C")) +
  labs(
    title    = "Section 6A — Distribution de τ1 et τ2 sur les émetteurs riches",
    subtitle = "CV élevé → τ varie beaucoup entre émetteurs → ne pas le fixer globalement",
    x = "Valeur de τ", y = "Nombre d'émetteurs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
print(p6a)

# ── 6B. Matrice de corrélation des paramètres estimés
if (nrow(nss_params_rich) >= 5) {
  cor_mat <- cor(nss_params_rich[, c("b0","b1","b2","b3","tau1","tau2")],
                 use = "complete.obs")
  message("\n── Matrice de corrélation des paramètres NSS ──")
  print(round(cor_mat, 2))
  
  # Visualisation
  cor_long <- as.data.frame(as.table(cor_mat)) |>
    rename(Param1 = Var1, Param2 = Var2, Correlation = Freq)
  
  p6b <- ggplot(cor_long, aes(x = Param1, y = Param2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = c(-1, 1), name = "Corrélation") +
    labs(
      title    = "Section 6B — Corrélations entre paramètres NSS (émetteurs riches)",
      subtitle = "Corrélations fortes → contraintes à imposer dans le modèle hiérarchique",
      x = NULL, y = NULL
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          axis.text  = element_text(size = 11, face = "bold"))
  print(p6b)
}

# ==============================================================================
# SECTION 7 — SYNTHÈSE : TABLEAU DE DÉCISION
# ==============================================================================
# Ce tableau récapitule les signaux obtenus et leur implication méthodologique.

decision_inputs <- list(
  pct_sparse_lt6   = round(100 * mean(sparsity_per_ticker$N_bonds < 6), 1),
  pct_rich_gte10   = round(100 * mean(sparsity_per_ticker$N_bonds >= 10), 1),
  eta_sq           = round(eta_squared, 3),
  tau1_cv          = round(sd(nss_params_rich$tau1)/mean(nss_params_rich$tau1), 2),
  tau2_cv          = round(sd(nss_params_rich$tau2)/mean(nss_params_rich$tau2), 2),
  median_cv_time   = round(median(spread_vol_time$CV_Time, na.rm=TRUE), 3),
  n_groups_non_empty = sum(group_stats$N_bonds >= 3)
)

message("\n", paste(rep("=", 60), collapse = ""))
message("SECTION 7 — TABLEAU DE SYNTHÈSE POUR LE CHOIX DE MÉTHODE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("%-45s %s", "% émetteurs < 6 bonds (sparse):",
                paste0(decision_inputs$pct_sparse_lt6, "%")))
message(sprintf("%-45s %s", "% émetteurs ≥ 10 bonds (riches):",
                paste0(decision_inputs$pct_rich_gte10, "%")))
message(sprintf("%-45s %s", "η² inter-émetteur (signal groupes):",
                decision_inputs$eta_sq))
message(sprintf("%-45s %s", "CV temporel médian (stabilité spreads):",
                decision_inputs$median_cv_time))
message(sprintf("%-45s %s", "CV de τ1 (fixe ou libre?):",
                decision_inputs$tau1_cv))
message(sprintf("%-45s %s", "CV de τ2 (fixe ou libre?):",
                decision_inputs$tau2_cv))
message(sprintf("%-45s %s", "Groupes rating×secteur avec ≥ 3 bonds:",
                decision_inputs$n_groups_non_empty))
message(paste(rep("-", 60), collapse = ""))
message("IMPLICATIONS :")
if (decision_inputs$pct_sparse_lt6 > 50) {
  message("  → MAJORITÉ sparse : le shrinkage est INDISPENSABLE")
} else {
  message("  → Majorité modérée/riche : la Ridge simple peut suffire")
}
if (decision_inputs$eta_sq > 0.5) {
  message("  → Fort signal inter-émetteur : les groupes sont informatifs")
} else {
  message("  → Signal inter-émetteur modéré : a priori groupes à valider")
}
if (decision_inputs$median_cv_time < 0.15) {
  message("  → Spreads stables : modèle statique OK (pas besoin d'AR1)")
} else {
  message("  → Spreads mobiles : envisager une composante dynamique")
}
if (max(decision_inputs$tau1_cv, decision_inputs$tau2_cv) > 0.5) {
  message("  → τ très variable : NE PAS fixer globalement")
} else {
  message("  → τ stable : fixer globalement pour réduire la dimensionnalité")
}
message(paste(rep("=", 60), collapse = ""))