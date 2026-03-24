# ==============================================================================
# IMPLÉMENTATION DE LA MÉTHODOLOGIE REN (2008)
# "Corporate Bond Yield Curve Estimation using Parametric Approach"
# ==============================================================================
#
# Logique : Estimation STATIQUE — une courbe par date d'observation.
# Décomposition : Taux Corporate(t) = Taux Swap Risk-Free(t) + Spread(t)
#
# Étape 1 : Calibrer la courbe Svensson sur les taux swap (risk-free)
#           (L'article utilise la courbe souveraine, mais ici les corpo
#            viennent de plusieurs pays, donc j'utilise la courbe Swap EUR
#            comme référence pour la courbe Risk Free)
# Étape 2 : Pour chaque date, calibrer une courbe Svensson sur le SPREAD
#           (par rapport à la courbe Swap) en minimisant la RMSYE.
# Étape 3 : Reconstruction de la courbe corporate totale + diagnostics
#
# Note sur les adaptations par rapport à l'article :
#   - Input = yields observés (pas des prix dirty) car nos données Bloomberg
#     fournissent directement le "Mid Yield to Convention". On minimise donc
#     la RMSYE directement sur les yields (équivalent asymptotique, cf. note
#     méthodologique dans le code).
#   - Courbe risk-free : taux swap EUR (unique référence pour obligations
#     libellées en EUR multi-pays), cohérent avec le principe de Ren.
#   - Pas de filtre sur une seule notation : on stratifie par rating bucket
#     (IG / HY) pour préserver la richesse des données tout en respectant
#     l'homogénéité du risque recommandée par l'article.
# ==============================================================================

rm(list = ls())

# ── 0. PACKAGES ET UTILITAIRES ────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(scales)
source("utils.R")   # nss_func(), nss_objective(), calculate_duration(), etc.

set.seed(42)  # Reproductibilité des optimisations

# Palette cohérente avec le reste du projet
rating_colors <- c(
  "AA"  = "#1B4F8A",
  "A"   = "#2E86C1",
  "BBB" = "#F39C12",
  "HY"  = "#C0392B"
)

# ==============================================================================
# SECTION 1 : CHARGEMENT ET PRE-TRAITEMENT DES DONNEES
# ==============================================================================

# ── 1.1 Données Corporate ─────────────────────────────────────────────────────

corpo_raw <- read_csv("r_data/NewData_Corpo_clean.csv", show_col_types = FALSE) |>
  mutate(across(starts_with("YLD_"), ~ suppressWarnings(as.numeric(gsub(",", ".", .x)))))

# Colonnes yield disponibles (labels de date au format DDMMYYYY)
yld_cols <- grep("^YLD_", names(corpo_raw), value = TRUE)

# Conversion en format long : 1 ligne = 1 obligation × 1 date
corpo_long <- corpo_raw |>
  pivot_longer(
    cols      = all_of(yld_cols),
    names_to  = "Date_Label",
    values_to = "Yield_Corp"
  ) |>
  mutate(
    # Suppression du préfixe YLD_ et conversion directe avec dmy() de lubridate
    Date = dmy(gsub("^YLD_", "", Date_Label)),
    
    Maturity     = as.numeric(Residual_Maturity),
    Amt_Out      = suppressWarnings(as.numeric(`Amt Out`)),
    Cpn          = suppressWarnings(as.numeric(Cpn)),
    Yield_Corp   = as.numeric(Yield_Corp)
  ) |>
  filter(!is.na(Yield_Corp), !is.na(Maturity), Maturity > 0.25, Maturity <= 30) |>
  rename(Rating = `BBG Composite`, Country_Risk = `Cntry of Risk`,
         Sector_L1 = `BICS Level 1`, Sector_L2 = `BICS Level 2`) |>
  select(ISIN, `Issuer Name`, Ticker, Rating, Sector_L1, Sector_L2,
         Country_Risk, Currency, Date, Maturity, Yield_Corp, Amt_Out, Cpn)

message("Dataset corporate long : ", nrow(corpo_long), " obs | ",
        length(unique(corpo_long$Date)), " dates | ",
        length(unique(corpo_long$ISIN)), " obligations")

# ── 1.2 Segmentation par Rating Bucket ────────────────────────────────────────
# Ren (2008) travaille sur les AA uniquement. Nos données sont principalement IG.
# On crée des buckets homogènes pour respecter l'esprit de l'article
# (analyser des titres comparables ensemble) tout en conservant la richesse des données.
#
# Justification : avec seulement 1 obligation AA et 673 au total, appliquer
# le filtre strict de Ren viderait l'analyse. On opère donc par strate.

corpo_long <- corpo_long |>
  mutate(
    Rating_Bucket = case_when(
      Rating %in% c("AAA", "AA+", "AA", "AA-")  ~ "AA",
      Rating %in% c("A+",  "A",  "A-")           ~ "A",
      Rating %in% c("BBB+", "BBB", "BBB-")        ~ "BBB",
      TRUE                                         ~ "HY"   # BB et inférieurs
    )
  )

# Vue d'ensemble de la segmentation
rating_summary <- corpo_long |>
  filter(Date == max(Date)) |>
  count(Rating_Bucket, Rating, name = "N_bonds") |>
  arrange(Rating_Bucket, Rating)
message("\n── Segmentation par rating bucket (dernière date) ──")
print(rating_summary)

# ── 1.3 Données Swap (Risk-Free) ───────────────────────────────────────────────

swap_raw <- read_csv("r_data/StatApp_Data1_Swap.csv", show_col_types = FALSE)

# Nettoyage du swap (repris de spread_analysis.R avec corrections déjà effectuées)
clean_swap_data <- function(df) {
  convert_tenor <- function(t) {
    if (grepl("M", t))  return(as.numeric(gsub("M", "", t)) / 12)
    if (grepl("Y", t))  return(as.numeric(gsub("Y", "", t)))
    return(NA)
  }
  df_clean <- df |>
    select(TenorStr = 1, SwapRate = 3) |>
    mutate(
      Maturity = sapply(TenorStr, convert_tenor),
      Yield    = suppressWarnings(as.numeric(as.character(SwapRate)))
    ) |>
    filter(!is.na(Maturity), !is.na(Yield)) |>
    arrange(Maturity)
  
  # Auto-correction d'échelle
  if (mean(df_clean$Yield, na.rm = TRUE) < 0.5)
    df_clean$Yield <- df_clean$Yield * 100
  
  df_clean
}

df_swap <- clean_swap_data(swap_raw)
message("\nCourbe swap : ", nrow(df_swap), " points, maturités de ",
        min(df_swap$Maturity), " à ", max(df_swap$Maturity), " ans")

# ==============================================================================
# SECTION 2 : CALIBRATION DE LA COURBE RISK-FREE (SWAP)
# ==============================================================================
# Cette étape correspond à la "Couche 1" de Ren (2008) :
# on estime d'abord la courbe Svensson sur l'actif de référence (swap EUR),
# qui joue ici le rôle des Gilts/OATs dans l'article original.
# La courbe swap est unique pour toutes les dates car nos données swap
# ne sont pas time-series (snapshot statique de la courbe).

# ── 2.1 Paramètres initiaux (Table III de Ren : 4%, -0.5%, 10%, -5%, 2.5, 0.7) ──
swap_start <- c(b0 = 4, b1 = -0.5, b2 = 10, b3 = -5, tau1 = 2.5, tau2 = 0.7)

# ── 2.2 Optimisation (L-BFGS-B avec contraintes explicites) ──
# Contraintes : b0 > 0 (taux long terme positif), tau > 0
swap_opt <- optim(
  par     = swap_start,
  fn      = nss_objective,
  maturities = df_swap$Maturity,
  yields     = df_swap$Yield,
  weights    = rep(1, nrow(df_swap)),   # Poids égaux sur le swap (Ren, Section III)
  method  = "L-BFGS-B",
  lower   = c(0,  -15, -15, -15, 0.1, 0.1),
  upper   = c(20,  15,  15,  15, 15,  30),
  control = list(maxit = 5000, factr = 1e6)
)

p_swap <- swap_opt$par
message("\n── Paramètres Svensson — Courbe Swap (Risk-Free) ──")
message("b0=", round(p_swap[1],3), " b1=", round(p_swap[2],3),
        " b2=", round(p_swap[3],3), " b3=", round(p_swap[4],3),
        " tau1=", round(p_swap[5],3), " tau2=", round(p_swap[6],3))
message("RMSE swap : ",
        round(sqrt(mean((df_swap$Yield - nss_func(df_swap$Maturity, p_swap))^2)), 4), "%")

# ── 2.3 Visualisation du fit swap ─────────────────────────────────────────────
grid_mat <- seq(0.25, 30, length.out = 300)

p_swap_fit <- ggplot() +
  geom_point(data = df_swap, aes(x = Maturity, y = Yield),
             size = 3, color = "black", shape = 16, alpha = 0.8) +
  geom_line(data = data.frame(Maturity = grid_mat,
                              Yield    = nss_func(grid_mat, p_swap)),
            aes(x = Maturity, y = Yield), color = "#2E86C1", linewidth = 1.2) +
  labs(
    title    = "Étape 1 (Ren 2008) — Courbe Svensson sur le Swap EUR (Risk-Free)",
    subtitle = paste0("RMSE = ",
                      round(sqrt(mean((df_swap$Yield - nss_func(df_swap$Maturity, p_swap))^2)), 4),
                      "% | Points = données swap observées | Ligne = modèle Svensson"),
    x = "Maturité (Années)", y = "Taux (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(p_swap_fit)

# ==============================================================================
# SECTION 3 : CALIBRATION DE LA COURBE DE SPREAD CORPORATE (COUCHE 2)
# ==============================================================================
# Correspond directement à la "Couche 2" de Ren (2008) :
# on applique une DEUXIÈME fonction Svensson au spread résiduel
#   Spread(t) = Yield_Corporate(t) - Yield_Swap(t)
# puis on minimise la RMSYE sur ce spread.
#
# Note méthodologique sur la RMSYE vs prix :
# Ren calcule théoriquement un prix modélisé puis reconvertit en yield.
# Pour des obligations à coupon non nul, la différence par rapport à une
# minimisation directe sur les yields est de l'ordre de Duration × RMSE_yield,
# soit typiquement < 1bp pour nos niveaux de RMSE. On minimise directement
# sur les yields, ce qui est standard dans la littérature empirique
# (cf. BIS Paper No. 25, 2005, méthode Bundesbank) et évite le problème
# d'identification numérique lié au calcul du prix dirty à partir d'un coupon
# annuel alors que nos obligations peuvent avoir diverses fréquences.

# ── 3.1 Calcul du spread observé pour chaque obligation × date ────────────────
corpo_long <- corpo_long |>
  mutate(
    # Taux swap interpolé à la maturité exacte de chaque obligation
    Yield_Swap   = nss_func(Maturity, p_swap),
    # Spread = ce que Ren appelle la "corporate spread" à modéliser
    Spread_Obs   = Yield_Corp - Yield_Swap
  )

# Diagnostic : vérifier la cohérence économique du spread
spread_check <- corpo_long |>
  filter(Date == max(Date)) |>
  group_by(Rating_Bucket) |>
  summarise(
    N            = n(),
    Spread_Mean  = round(mean(Spread_Obs, na.rm = TRUE), 3),
    Spread_Min   = round(min(Spread_Obs,  na.rm = TRUE), 3),
    Spread_Max   = round(max(Spread_Obs,  na.rm = TRUE), 3),
    .groups = "drop"
  )
message("\n── Spreads observés par bucket (dernière date) ──")
print(spread_check)
# Attendu : AA < A < BBB < HY (hiérarchie de crédit)

# ── 3.2 Fonction de filtrage outliers (Ren recommande une sélection rigoureuse) ──
# L'article insiste sur la sélection de l'échantillon comme "vital issue".
# On réutilise filter_ecb_criteria() de utils.R adapté au spread.
# Cette fonction retire les obligations dont le spread est à > 2σ du spread
# médian dans leur bucket de maturité — cohérent avec l'approche de Ren
# qui exclut les obligations "atypiques" de son échantillon AA.

filter_spread_outliers <- function(df, yield_col = "Spread_Obs") {
  df <- df |> rename(.y = all_of(yield_col))
  df_filtered <- df |> filter(.y > -2, .y < 20)  # Spread corporates : [−200bps, 2000bps]
  
  for (i in 1:2) {  # Deux passages, comme la BCE
    n_buckets <- max(1, floor(nrow(df_filtered) / 12))
    df_filtered <- df_filtered |>
      mutate(Bucket = ggplot2::cut_number(Maturity, n = n_buckets)) |>
      group_by(Bucket) |>
      mutate(
        mu    = mean(.y, na.rm = TRUE),
        sigma = sd(.y,   na.rm = TRUE),
        Lower = mu - 2.5 * sigma,
        Upper = mu + 2.5 * sigma
      ) |>
      ungroup() |>
      filter(.y >= Lower & .y <= Upper) |>
      select(-Bucket, -mu, -sigma, -Lower, -Upper)
  }
  df_filtered |> rename(!!yield_col := .y)
}

# ── 3.3 Fonction d'estimation de la courbe de spread (cœur de Ren 2008) ───────
#
# Paramètres initiaux pour le spread selon Ren (Table : spread moyen ~80bps,
# légèrement décroissant avec la maturité pour l'IG).
# On adapte dynamiquement à partir de la moyenne empirique.

fit_spread_curve <- function(df_spread_obs, rating_bucket) {
  
  # Sécurité : il faut au moins 6 points pour identifier 6 paramètres
  if (nrow(df_spread_obs) < 6) return(NULL)
  
  # Filtrage des outliers de spread (sélection d'échantillon de Ren)
  df_clean <- filter_spread_outliers(df_spread_obs)
  if (nrow(df_clean) < 6) return(NULL)
  
  mats    <- df_clean$Maturity
  spreads <- df_clean$Spread_Obs
  wts     <- pmax(df_clean$Amt_Out, 0, na.rm = TRUE)
  if (sum(wts, na.rm = TRUE) == 0) wts <- rep(1, length(spreads))
  wts <- wts / sum(wts)
  
  # Paramètres initiaux économiquement fondés (Ren Section IV)
  mean_spread <- mean(spreads, na.rm = TRUE)
  start <- c(
    b0   = max(mean_spread, 0.1),  # niveau long terme du spread (> 0 pour l'IG)
    b1   = -0.5,                   # spread légèrement décroissant (Ren: β1 < 0)
    b2   =  10,                    # courbure (Ren: initialisation large)
    b3   = -5,                     # contre-courbure (Ren)
    tau1 =  2.5,                   # même initialisation que dans l'article
    tau2 =  0.7
  )
  
  # Optimisation : minimisation de la RMSYE sur le spread
  # (équivalent de la procédure de Ren, Section III, step 6)
  result <- tryCatch({
    optim(
      par     = start,
      fn      = function(params) {
        # Contraintes souples (Ren impose b0 > 0, tau > 0)
        if (params[1] < -2 || params[5] <= 0.05 || params[6] <= 0.05) return(1e9)
        fitted <- nss_func(mats, params)
        sqrt(sum(wts * (spreads - fitted)^2, na.rm = TRUE) / sum(wts))
      },
      method  = "L-BFGS-B",
      lower   = c(-2,  -15, -15, -15, 0.05, 0.05),
      upper   = c(20,   15,  15,  15, 15,   30),
      control = list(maxit = 3000, factr = 1e7)
    )
  }, error = function(e) NULL)
  
  if (is.null(result) || result$convergence > 1) return(NULL)
  
  list(
    params        = result$par,
    rmsye         = result$value,
    n_bonds       = nrow(df_clean),
    n_bonds_raw   = nrow(df_spread_obs),
    data_clean    = df_clean
  )
}

# ── 3.4 Boucle principale : estimation par date × rating bucket ───────────────

dates_obs   <- sort(unique(corpo_long$Date))
buckets_obs <- c("AA", "A", "BBB", "HY")

results_all <- list()
pb_total    <- length(dates_obs) * length(buckets_obs)
counter     <- 0

message("\n── Estimation Ren (2008) : ", length(dates_obs), " dates × ",
        length(buckets_obs), " buckets ──")

for (d in dates_obs) {
  for (bucket in buckets_obs) {
    
    counter <- counter + 1
    
    df_sub <- corpo_long |>
      filter(Date == d, Rating_Bucket == bucket) |>
      filter(!is.na(Spread_Obs))
    
    fit <- fit_spread_curve(df_sub, bucket)
    
    key <- paste(as.character(d), bucket, sep = "_")
    results_all[[key]] <- list(
      Date          = d,
      Rating_Bucket = bucket,
      fit           = fit,
      n_raw         = nrow(df_sub)
    )
  }
}

message("Estimation terminée. ", sum(sapply(results_all, function(x) !is.null(x$fit))),
        " / ", pb_total, " modèles convergés.")

# ==============================================================================
# SECTION 4 : EXTRACTION DES RÉSULTATS ET DIAGNOSTICS
# ==============================================================================

# ── 4.1 Table des paramètres estimés (équivalent Table II de Ren) ─────────────

params_df <- bind_rows(lapply(results_all, function(x) {
  if (is.null(x$fit)) {
    return(data.frame(Date = x$Date, Rating_Bucket = x$Rating_Bucket,
                      b0 = NA, b1 = NA, b2 = NA, b3 = NA,
                      tau1 = NA, tau2 = NA, RMSYE = NA,
                      N_bonds_clean = NA, N_bonds_raw = x$n_raw))
  }
  p <- x$fit$params
  data.frame(
    Date          = x$Date,
    Rating_Bucket = x$Rating_Bucket,
    b0   = p[1], b1 = p[2], b2 = p[3], b3 = p[4],
    tau1 = p[5], tau2 = p[6],
    RMSYE         = x$fit$rmsye,
    N_bonds_clean = x$fit$n_bonds,
    N_bonds_raw   = x$n_raw
  )
})) |>
  arrange(Rating_Bucket, Date)

message("\n── Paramètres moyens par Rating Bucket (Tableau de Ren) ──")
print(
  params_df |>
    group_by(Rating_Bucket) |>
    summarise(across(c(b0, b1, b2, b3, tau1, tau2, RMSYE),
                     ~ round(mean(.x, na.rm = TRUE), 4)),
              .groups = "drop")
)

# ── 4.2 Reconstruction des courbes de spread et de taux total ─────────────────

# Pour chaque date × bucket avec un fit valide, on calcule :
# Spread(t) = NSS(t, params_spread)
# Yield_Corp_Model(t) = Yield_Swap(t) + Spread(t)

reconstruct_curves <- function(params, p_swap_ref, mats = seq(0.25, 30, by = 0.1)) {
  spread_curve <- nss_func(mats, params)
  rf_curve     <- nss_func(mats, p_swap_ref)
  data.frame(
    Maturity     = mats,
    Spread_Model = spread_curve,
    Yield_RF     = rf_curve,
    Yield_Total  = rf_curve + spread_curve
  )
}

curves_all <- bind_rows(lapply(results_all, function(x) {
  if (is.null(x$fit)) return(NULL)
  curves <- reconstruct_curves(x$fit$params, p_swap)
  curves$Date          <- x$Date
  curves$Rating_Bucket <- x$Rating_Bucket
  curves
}))

# ── 4.3 Calcul des résidus in-sample (RMSYE par obligation) ──────────────────

residuals_df <- bind_rows(lapply(results_all, function(x) {
  if (is.null(x$fit)) return(NULL)
  df_c <- x$fit$data_clean
  df_c$Spread_Fitted <- nss_func(df_c$Maturity, x$fit$params)
  df_c$Yield_Fitted  <- df_c$Yield_Swap + df_c$Spread_Fitted
  df_c$Residual_Yld  <- df_c$Yield_Corp - df_c$Yield_Fitted
  df_c$Residual_Spd  <- df_c$Spread_Obs - df_c$Spread_Fitted
  df_c$Date          <- x$Date
  df_c$Rating_Bucket <- x$Rating_Bucket
  df_c
}))

# Performance globale par bucket (Table I équivalent)
perf_summary <- residuals_df |>
  group_by(Rating_Bucket) |>
  summarise(
    N_obs         = n(),
    RMSYE_Yield   = round(sqrt(mean(Residual_Yld^2, na.rm = TRUE)), 4),
    RMSYE_Spread  = round(sqrt(mean(Residual_Spd^2, na.rm = TRUE)), 4),
    Bias_Yield    = round(mean(Residual_Yld,          na.rm = TRUE), 4),
    .groups = "drop"
  ) |>
  arrange(match(Rating_Bucket, c("AA", "A", "BBB", "HY")))

message("\n── Performance in-sample (RMSYE en %) — Ren Table I equivalent ──")
print(perf_summary)

# ── 4.3 Calcul des métriques de performance in-sample (Via utils.R) ──────────

# 1. Évaluation détaillée pour chaque modèle (Date x Bucket)
scorecards_list <- lapply(results_all, function(x) {
  if (is.null(x$fit)) return(NULL)
  
  # Préparation : on renomme Yield_Corp en Yield pour la fonction
  df_c <- x$fit$data_clean |>
    mutate(Yield = Yield_Corp)
  
  # Construction de la fonction de prédiction CONJOINTE (Svensson Swap + Svensson Spread)
  # C'est indispensable pour évaluer le "Yield_Total" et calculer la Wiggliness correctement
  pred_joint <- function(m) {
    nss_func(m, p_swap) + nss_func(m, x$fit$params)
  }
  
  # Calcul des métriques via ta fonction
  score <- suppressWarnings(
    evaluate_model_performance(
      model_name   = paste(x$Date, x$Rating_Bucket, sep = "_"),
      data         = df_c,
      predict_func = pred_joint,
      col_weight   = "Amt_Out"
    )
  )
  
  # Ajout des métadonnées pour l'agrégation
  score$Date <- x$Date
  score$Rating_Bucket <- x$Rating_Bucket
  
  return(score)
})

# 2. Assemblage de tous les tableaux de bord
all_scorecards <- bind_rows(scorecards_list)

# 3. Agrégation pour le rapport final (Moyenne par Rating Bucket)
perf_summary <- all_scorecards |>
  group_by(Rating_Bucket) |>
  summarise(
    N_Models        = n(),
    W_RMSE          = round(mean(W_RMSE, na.rm = TRUE), 4),
    Price_RMSE      = round(mean(Price_RMSE, na.rm = TRUE), 4),
    Wiggliness      = round(mean(Wiggliness, na.rm = TRUE), 4),
    Tail_Stab_30_50 = round(mean(Tail_Stab_30_50, na.rm = TRUE), 4),
    # On utilise any_of au cas où certains buckets de maturité n'existent pas sur toutes les dates
    across(any_of(c("Bias_0-2Y", "Bias_2-5Y", "Bias_5-10Y", "Bias_10-30Y", "Bias_30Y+")), 
           ~ round(mean(.x, na.rm = TRUE), 4)),
    .groups = "drop"
  ) |>
  arrange(match(Rating_Bucket, c("AA", "A", "BBB", "HY")))

message("\n── Performance in-sample (Agrégée par Rating Bucket) ──")
print(perf_summary)

# 4. (Optionnel mais requis pour ta VIZ 5) Recréation du residuals_df
residuals_df <- bind_rows(lapply(results_all, function(x) {
  if (is.null(x$fit)) return(NULL)
  df_c <- x$fit$data_clean
  df_c$Spread_Fitted <- nss_func(df_c$Maturity, x$fit$params)
  df_c$Yield_Fitted  <- df_c$Yield_Swap + df_c$Spread_Fitted
  df_c$Residual_Yld  <- df_c$Yield_Corp - df_c$Yield_Fitted
  df_c$Residual_Spd  <- df_c$Spread_Obs - df_c$Spread_Fitted
  df_c$Date          <- x$Date
  df_c$Rating_Bucket <- x$Rating_Bucket
  df_c
}))

# ==============================================================================
# SECTION 5 : VISUALISATIONS POUR LE RAPPORT
# ==============================================================================

# --- CORRECTION DE LA CLASSE DES DATES ---
curves_all   <- curves_all   |> mutate(Date = as.Date(Date, origin = "1970-01-01"))
params_df    <- params_df    |> mutate(Date = as.Date(Date, origin = "1970-01-01"))
residuals_df <- residuals_df |> mutate(Date = as.Date(Date, origin = "1970-01-01")) # Ligne ajoutée !

last_date  <- max(dates_obs)
first_date <- min(dates_obs)
mid_date   <- dates_obs[round(length(dates_obs) / 2)]

# ── VIZ 1 : Courbes de taux corporate totales par bucket (snapshot) ───────────
# Équivalent de la Figure principale de Ren : Corporate Yield Curve estimée

df_snapshot <- curves_all |>
  filter(Date == last_date) |>
  mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY")))

df_rf_line <- data.frame(
  Maturity = grid_mat,
  Yield    = nss_func(grid_mat, p_swap),
  label    = "Swap (Risk-Free)"
)

p_yield_curves <- ggplot() +
  # Courbe risk-free en référence
  geom_line(data = df_rf_line, aes(x = Maturity, y = Yield, linetype = label),
            color = "black", linewidth = 0.8) +
  # Courbes corporate par bucket
  geom_line(data = df_snapshot,
            aes(x = Maturity, y = Yield_Total, color = Rating_Bucket),
            linewidth = 1.1) +
  # Points observés sur la dernière date
  geom_point(
    data = residuals_df |> 
      filter(Date == last_date) |>
      mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY"))),
    aes(x = Maturity, y = Yield_Corp, color = Rating_Bucket),
    alpha = 0.35, size = 1.5, shape = 16
  ) +
  scale_color_manual(values = rating_colors, name = "Rating Bucket") +
  scale_linetype_manual(values = c("Swap (Risk-Free)" = "dashed"), name = NULL) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title    = paste0("Courbes de Taux Corporate"),
    subtitle = paste0("Date : ", format(last_date, "%d/%m/%Y"),
                      " | Points = obligations observées (outliers retirés) | Lignes = modèle NSS"),
    x = "Maturité (Années)", y = "Taux Corporate (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_yield_curves)

# ── VIZ 2 : Courbes de SPREAD par bucket (cœur de la contribution de Ren) ────
# C'est la "Spread Curve" que Ren paramétrise par une 2ème Svensson

p_spread_curves <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_line(
    data = df_snapshot |> filter(!is.na(Spread_Model)),
    aes(x = Maturity, y = Spread_Model, color = Rating_Bucket),
    linewidth = 1.1
  ) +
  # Points observés (spread brut)
  geom_point(
    data = residuals_df |>
      filter(Date == last_date) |>
      mutate(Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY"))),
    aes(x = Maturity, y = Spread_Obs, color = Rating_Bucket),
    alpha = 0.3, size = 1.5
  ) +
  scale_color_manual(values = rating_colors, name = "Rating Bucket") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title    = "Structure par Terme du Spread Corporate vs Swap",
    subtitle = paste0("Points = spreads observés | Lignes = NSS sur le spread"),
    x = "Maturité (Années)", y = "Spread vs Swap (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_spread_curves)

# ── VIZ 3 : Évolution temporelle des paramètres de spread (Table II de Ren) ───
# C'est LE graphique qui montre la richesse de l'approche : les paramètres
# capturent l'évolution du crédit dans le temps

p_params <- params_df |>
  filter(!is.na(b0)) |>
  select(Date, Rating_Bucket, b0, b1, b2, b3, RMSYE) |>
  pivot_longer(cols = c(b0, b1, b2, b3, RMSYE), names_to = "Param", values_to = "Value") |>
  mutate(
    Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY")),
    Param_Label   = recode(Param,
                           b0 = "β0 — Niveau long terme",
                           b1 = "β1 — Pente",
                           b2 = "β2 — Courbure 1",
                           b3 = "β3 — Courbure 2",
                           RMSYE = "RMSYE (%)"
    )
  ) |>
  ggplot(aes(x = Date, y = Value, color = Rating_Bucket)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_wrap(~Param_Label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = rating_colors, name = "Rating Bucket") +
  scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") +
  labs(
    title    = "Dynamique des Paramètres NSS sur le Spread Corporate",
    x = NULL, y = "Valeur du paramètre"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold", size = 9)
  )

print(p_params)

# ── VIZ 4 : Hiérarchie des courbes de spread à 3 dates clés ──────────────────
# Reproduit l'analyse de Ren (2008) section IV : évolution cross-temporelle

df_three_dates <- curves_all |>
  filter(Date %in% c(first_date, mid_date, last_date)) |>
  mutate(
    Date_Label    = format(Date, "%d/%m/%Y"),
    Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY"))
  )

p_cross_time <- ggplot(df_three_dates, aes(x = Maturity, y = Spread_Model,
                                           color = Rating_Bucket,
                                           linetype = Date_Label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dotted") +
  scale_color_manual(values = rating_colors, name = "Rating") +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    name   = "Date"
  ) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  facet_wrap(~Rating_Bucket, scales = "free_y", ncol = 2) +
  labs(
    title    = "Évolution Cross-Temporelle du Spread sur 3 Dates",
    x = "Maturité (Années)", y = "Spread vs Swap (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", color = "gray30")
  )

print(p_cross_time)

# ── VIZ 5 : Résidus in-sample par bucket (diagnostic de qualité) ──────────────

p_residuals <- residuals_df |>
  mutate(
    Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY")),
    Maturity_Bucket = cut(Maturity,
                          breaks = c(0, 2, 5, 10, 30),
                          labels = c("0–2Y", "2–5Y", "5–10Y", "10–30Y"))
  ) |>
  filter(!is.na(Residual_Yld)) |>
  ggplot(aes(x = Maturity, y = Residual_Yld, color = Rating_Bucket)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.5) +
  facet_wrap(~Rating_Bucket, ncol = 2, scales = "free_y") +
  scale_color_manual(values = rating_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title    = "Résidus In-Sample (Yield observé − Yield modélisé)",
    subtitle = "Courbe = tendance LOESS",
    x = "Maturité (Années)", y = "Résidu (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

print(p_residuals)

# ── VIZ 6 : Taux spot à maturités clés — Table I de Ren (2008) ───────────────
# Équivalent des Y1, Y3, Y5, Y7, Y10, Y15 de l'article, ici Y2, Y5, Y10

key_mats <- c(2, 5, 10)

key_rates <- bind_rows(lapply(results_all, function(x) {
  if (is.null(x$fit)) return(NULL)
  sapply(key_mats, function(m) {
    spread_val <- nss_func(m, x$fit$params)
    rf_val     <- nss_func(m, p_swap)
    rf_val + spread_val
  }) |>
    setNames(paste0("Y", key_mats)) |>
    as.list() |>
    as.data.frame() |>
    mutate(Date = x$Date, Rating_Bucket = x$Rating_Bucket)
}))

message("\n── Taux Corporate à maturités clés (Tableau I de Ren) ──")
print(
  key_rates |>
    mutate(across(starts_with("Y"), ~ round(.x, 3))) |>
    arrange(Rating_Bucket, Date)
)

p_key_rates <- key_rates |>
  pivot_longer(starts_with("Y"), names_to = "Tenor", values_to = "Yield") |>
  mutate(
    Rating_Bucket = factor(Rating_Bucket, levels = c("AA", "A", "BBB", "HY")),
    Tenor         = factor(Tenor, levels = paste0("Y", key_mats))
  ) |>
  ggplot(aes(x = Date, y = Yield, color = Rating_Bucket, linetype = Tenor)) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = rating_colors, name = "Rating") +
  scale_linetype_manual(
    values = c("Y2" = "solid", "Y5" = "dashed", "Y10" = "dotted"),
    name   = "Maturité"
  ) +
  scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") +
  facet_wrap(~Rating_Bucket, ncol = 2, scales = "free_y") +
  labs(
    title    = "Taux Spot Corporate aux Maturités 2Y, 5Y et 10Y",
    x = NULL, y = "Taux (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold")
  )

print(p_key_rates)