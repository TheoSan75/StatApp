rm(list = ls())

# 1. Setup ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(gridExtra)

colors_flags <- c(
  "France" = "#0055A4",  # Bleu officiel Français
  "Germany" = "#BA0C2F", # Rouge officiel allemand
  "Italy" = "#008C45"    # Vert officiel italien
)
# 2. Load Data -----------------------------------------------------------------
data_france <- read_csv("r_data/StatApp_Data1_France.csv", show_col_types = FALSE)
data_allemagne <- read_csv("r_data/StatApp_Data1_Allemagne.csv", show_col_types = FALSE)
data_italie <- read_csv("r_data/StatApp_Data1_Italie.csv", show_col_types = FALSE)
data_swaps <- read_csv("r_data/StatApp_Data1_Swap.csv", show_col_types = FALSE)


# 3. Fonctions utilitaires --------------------------------------

# Fonction pour calculer la "Wiggliness" (Rugosité) d'une courbe
calc_wiggliness <- function(predict_func, min_mat = 0, max_mat = 50, step = 0.25) {
  # On crée une grille de points artificiels (tous les 3 mois)
  grid_mat <- seq(min_mat, max_mat, by = step)
  # On prédit les taux sur cette grille
  grid_yields <- predict_func(grid_mat)
  # On calcule la dérivée seconde numérique (différences secondes)
  second_deriv <- diff(diff(grid_yields))
  # Métrique : Somme des carrés des dérivées secondes (standardisée par la taille)
  wiggliness_score <- sum(second_deriv^2) * 1000 # *1000 pour lisibilité
  return(wiggliness_score)
}

# Fonction générale pour calculer toutes les métriques
evaluate_model_performance <- function(model_name, data, predict_func, col_weight = "Amt_Out") {
  # Prédictions sur les données réelles
  preds <- predict_func(data$Maturity)
  residuals <- data$Yield - preds
  # Gestion des Poids (Si pas de poids ou poids nuls, on met 1)
  weights <- data[[col_weight]]
  if(is.null(weights) | all(is.na(weights))) { weights <- rep(1, nrow(data)) }
  weights[is.na(weights)] <- 0 # Sécurité
  
  # --- MÉTRIQUE 1 : RMSE Pondéré ---
  w_rmse <- sqrt(sum(weights * residuals^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))
  # --- MÉTRIQUE 2 : Wiggliness ---
  wiggliness <- calc_wiggliness(predict_func)
  # --- MÉTRIQUE 3 : Test de Queues (Ecart 30Y vs 50Y) ---
  pred_30 <- predict_func(30)
  pred_50 <- predict_func(50)
  tail_diff <- abs(pred_50 - pred_30)
  # --- MÉTRIQUE 4 : Analyse des Résidus par Bucket ---
  # On regarde si le modèle se trompe toujours dans le même sens sur certaines zones
  data_res <- data %>%
    mutate(
      Prediction = preds,
      Residual = residuals,
      Bucket = cut(Maturity, 
                   breaks = c(0, 2, 5, 10, 30, 100), 
                   labels = c("0-2Y", "2-5Y", "5-10Y", "10-30Y", "30Y+"))
    )
  # On calcule le Biais Moyen (Erreur moyenne) par bucket
  # Si Bias > 0 : Le modèle sous-estime le taux (Marché > Modèle)
  # Si Bias < 0 : Le modèle surestime le taux
  bucket_bias <- data_res %>%
    group_by(Bucket) %>%
    summarise(Bias = mean(Residual, na.rm = TRUE), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = Bucket, values_from = Bias, names_prefix = "Bias_")
  
  # 3. Assemblage du "Scorecard" final
  scorecard <- tibble(
    Model = model_name,
    W_RMSE = round(w_rmse, 4),     # Précision globale
    Wiggliness = round(wiggliness, 4), # Stabilité
    Tail_Stab_30_50 = round(tail_diff, 4) # Cohérence long terme
  ) %>%
    bind_cols(bucket_bias) # Ajoute les colonnes de biais (Bias_0-2Y, etc.)
  
  return(scorecard)
}

filter_ecb_criteria <- function(df) {
  # 1. Maturity Spectrum: 3 months (0.25Y) to 30 Years
  df_filtered <- df %>%
    filter(Maturity >= 0.25 & Maturity <= 30)
  
  # 2. Iterative Outlier Removal (2 Sigma Rule per bucket)
  clean_data <- df_filtered %>%
    mutate(Bucket = cut(Maturity, breaks = c(0, 2, 5, 10, 30), include.lowest = TRUE))
  
  for(i in 1:2) {
    clean_data <- clean_data %>%
      group_by(Bucket) %>%
      mutate(mu = mean(Yield, na.rm=TRUE), sigma = sd(Yield, na.rm=TRUE)) %>%
      ungroup() %>%
      filter(abs(Yield - mu) <= 2*sigma | sigma < 0.05) %>%
      select(-mu, -sigma)
  }
  return(clean_data)
}

# 4. Data Cleaning & Filtering -------------------------------------------------
clean_bond_data <- function(df, country_name, ref_date) {
  
  # Filtrage des instruments (Bonds nominaux uniquement)
  if (country_name == "France") {
    df_clean <- df %>% filter(Series == "OAT")
  } else if (country_name == "Germany") {
    df_clean <- df %>% filter(!Series %in% c("G", "TWIN")) 
  } else if (country_name == "Italy") {
    df_clean <- df %>% filter(!Series %in% c("CPI", "ICPI", "PIU", "VALR", "FUT"))
  }
  
  # Renommage et Conversion
  df_clean <- df_clean %>%
    rename(
      RawYield = `Mid Yield to Convention`, 
      MaturityStr = Maturity,
      RawAmt = `Amt Out`,
      RawCpn = Cpn
    ) %>%
    mutate(
      Yield = suppressWarnings(as.numeric(as.character(RawYield))),
      Amt_Out = suppressWarnings(as.numeric(as.character(RawAmt))),
      Cpn = suppressWarnings(as.numeric(as.character(RawCpn))),
      MaturityDate = as.Date(MaturityStr, format = "%d/%m/%Y"),
      Country = country_name,
      # Calcul de la maturité en années
      Maturity = as.numeric(MaturityDate - ref_date) / 365.25 
    ) %>%
    # Nettoyage des NA et valeurs aberrantes
    filter(!is.na(Yield), !is.na(Maturity), !is.na(Cpn), Maturity > 0) %>%
    select(Country, MaturityDate, Maturity, Yield, Amt_Out, Cpn) 
  
  return(df_clean)
}

# --- Préparation de la date d'analyse ---
# On convertit temporairement en date pour trouver le min global
dates_vec <- c(
  as.Date(data_france$Maturity, format="%d/%m/%Y"), 
  as.Date(data_allemagne$Maturity, format="%d/%m/%Y"),
  as.Date(data_italie$Maturity, format="%d/%m/%Y")
)
analysis_date <- min(dates_vec, na.rm=TRUE) - 2
message("Date d'analyse fixée au : ", analysis_date)

# --- Création du dataset global nettoyé ---
df_all <- bind_rows(
  clean_bond_data(data_france, "France", analysis_date),
  clean_bond_data(data_allemagne, "Germany", analysis_date),
  clean_bond_data(data_italie, "Italy", analysis_date)
)


# 5. Modèle naïf: Regression linéaire ------------------------------------------

# Liste pour stocker les résultats
plots_reglin <- list()
scorecards_reglin <- list()

data_list <- split(df_all, df_all$Country)

for (country in names(data_list)) {
  df <- data_list[[country]]
  model_lin <- lm(Yield ~ Maturity, data = df)

  # L'évaluateur attend une fonction func(maturité) -> rendement
  predict_lin <- function(m) {
    predict(model_lin, newdata = data.frame(Maturity = m))
  }

  # Calcul des Métriques
  score <- evaluate_model_performance(
    model_name = paste("Linear_", country),
    data = df,
    predict_func = predict_lin,
    col_weight = "Amt_Out"
  )
  scorecards_reglin[[country]] <- score
  
  # Visualisation
  grid_x <- seq(0, 50, length.out = 500)
  grid_y <- predict_lin(grid_x)
  line_data <- data.frame(Maturity = grid_x, Yield = grid_y)
  
  p <- ggplot(df, aes(x = Maturity, y = Yield)) +
    # 1. Les points observés (transparents pour voir la densité)
    geom_point(color = colors_flags[country], alpha = 0.6, size = 2) +
    # 2. Le Modèle (La Droite)
    geom_line(data = line_data, aes(x = Maturity, y = Yield), 
              color = "black", linetype = "dashed", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    # 3. Méta-informations
    labs(
      title = paste(country, "- Modèle Linéaire"),
      subtitle = paste("RMSE Pondéré:", score$W_RMSE, "| Biais 10-30Y:", score$`Bias_10-30Y`),
      x = "Maturité (Années)", 
      y = "Rendement (Yield %)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", color = colors_flags[country]))
  
  # Graphique des résidus pour montrer le fameux "U inversé"
  df$Residus <- df$Yield - predict(model_lin, df)
  
  p_res <- ggplot(df, aes(x = Maturity, y = Residus)) +
    geom_point(color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black") +
    geom_smooth(method = "loess", se = FALSE, color = "gray", linetype = "dotted") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Structure des Résidus", y = "Résidu (bps)", x = "Maturité") +
    theme_minimal()
  
  # On combine les deux graphes (Haut : Fit, Bas : Résidus)
  plots_reglin[[country]] <- grid.arrange(p, p_res, nrow = 2, heights = c(2, 1))
}

# Affichage du tableau récapitulatif
final_scorecard_reglin <- bind_rows(scorecards_reglin)
print(final_scorecard_reglin)


# 6. Modèle NSS OLS ------------------------------------------------------------

# Définition de la fonction NSS
# b0 : Taux long terme (Level)
# b1 : Taux court - Taux long (Slope)
# b2 : Courbure à moyen terme (Hump 1)
# b3 : Courbure à long terme (Hump 2)
# tau1, tau2 : Facteurs d'échelle temporelle
nss_func <- function(m, params) {
  b0 <- params[1]
  b1 <- params[2]
  b2 <- params[3]
  b3 <- params[4]
  tau1 <- params[5]
  tau2 <- params[6]
  term1 <- (1 - exp(-m/tau1)) / (m/tau1)
  term2 <- term1 - exp(-m/tau1)
  term3 <- (1 - exp(-m/tau2)) / (m/tau2) - exp(-m/tau2)
  
  return(b0 + b1*term1 + b2*term2 + b3*term3)
}

# Fonction d'optimisation (OLS = Minimiser la somme des carrés des erreurs)

fit_nss_ols <- function(df) {
  # Fonction Objective : Somme des Carrés des Résidus (SSR)
  ssr_cost <- function(params) {
    # Contraintes : tau doit être positif
    if(params[5] <= 0 || params[6] <= 0) {
      return(1e9)
    }
    
    preds <- nss_func(df$Maturity, params)
    sum((df$Yield - preds)^2) # OLS : Pas de pondération par Amt_Out
  }
  
  # Paramètres initiaux (Guessing intelligent pour aider l'algo)
  start_params <- c(b0=3, b1=-1, b2=0, b3=0, tau1=2, tau2=10)
  
  # Optimisation via Nelder-Mead (Robuste)
  opt <- optim(par = start_params, fn = ssr_cost, control = list(maxit = 5000))
  return(opt$par)
}

# Exécution de a méthode
plots_nssols <- list()
scorecards_nssols <- list()
models_nssols_params <- list() # On stocke les paramètres pour plus tard

for (country in names(data_list)) {
  df <- data_list[[country]]
  best_params <- fit_nss_ols(df)
  models_nssols_params[[country]] <- best_params
  
  # Wrapper de prédiction
  predict_nss <- function(m){
    nss_func(m, best_params)
    }
  
  # Calcul des métriques
  score <- evaluate_model_performance(
    model_name = paste("NSS_OLS_", country),
    data = df,
    predict_func = predict_nss,
    col_weight = "Amt_Out"
  )
  scorecards_nssols[[country]] <- score
  
  # Visualisation
  grid_x <- seq(0, 50, length.out = 500)
  line_data <- data.frame(Maturity = grid_x, Yield = predict_nss(grid_x))
  
  # Graphe Principal
  p <- ggplot(df, aes(x = Maturity, y = Yield)) +
    # Points : La taille dépend du Volume (Amt_Out) pour voir les "Trash Bonds"
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    # La Courbe NSS
    geom_line(data = line_data, aes(x = Maturity, y = Yield), 
              color = "black", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(
      title = paste(country, "- NSS Classique (OLS)"),
      subtitle = paste("Le modèle traite le point petit (Trash) comme le gros (Benchmark)."),
      x = "Maturité", y = "Yield (%)", size = "Volume"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", color = colors_flags[country]))
  
  # Graphe des Résidus
  df$Residus <- df$Yield - predict_nss(df$Maturity)
  
  p_res <- ggplot(df, aes(x = Maturity, y = Residus)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0) +
    # On ajoute une ligne de tendance locale pour voir si le biais persiste
    geom_smooth(method = "loess", se = FALSE, color = "gray", linetype = "dotted") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Résidus : Hétéroscédasticité visible", y = "Bps", x = "Maturité", size = "Volume") +
    theme_minimal()
  
  plots_nssols[[country]] <- grid.arrange(p, p_res, nrow = 2, heights = c(2, 1))
}

# Affichage des Résultats
final_scorecard_nssols <- bind_rows(scorecards_nssols)
print(final_scorecard_nssols)

# Comparaison RegLin vs NSS_OLS
print("Comparaison rapide RMSE (Linear vs NSS OLS):")
print(bind_rows(final_scorecard_reglin %>% select(Model, W_RMSE), 
                final_scorecard_nssols %>% select(Model, W_RMSE)))


# ==============================================================================
# 6.2 INTERLUDE: GLOBAL ECB DATA PRE-PROCESSING & VISUALIZATION
# ==============================================================================

# 1. Apply Filtering Globally
# We create a temporary cleaned list to compare
data_list_clean <- lapply(data_list, filter_ecb_criteria)

# 2. Visualize Removed Data (Before vs After)
plot_cleaning_impact <- function(country) {
  raw <- data_list[[country]] %>% mutate(Status = "Removed")
  clean <- data_list_clean[[country]] %>% mutate(Status = "Kept")
  
  # We identify points that are in 'raw' but NOT in 'clean'
  combined <- bind_rows(clean, anti_join(raw, clean, by = c("Maturity", "Yield")))
  
  ggplot(combined, aes(x = Maturity, y = Yield, color = Status, shape = Status)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Kept" = "#0055A4", "Removed" = "red")) +
    scale_shape_manual(values = c("Kept" = 19, "Removed" = 4)) +
    labs(title = paste("ECB Cleaning Impact -", country),
         subtitle = "Red crosses = Outliers or maturities outside 3M-30Y range") +
    theme_minimal()
}

# Generate plots
p_clean_fr <- plot_cleaning_impact("France")
p_clean_de <- plot_cleaning_impact("Germany")
p_clean_it <- plot_cleaning_impact("Italy")

grid.arrange(p_clean_fr, p_clean_de, p_clean_it, nrow = 2)

# 3. OVERWRITE GLOBAL DATA
# From this point on, 'data_list' contains ONLY the cleaned data.
# Models 7, 8, and 9 will automatically use this cleaned data.
data_list <- data_list_clean 
message("⚠️ DATA UPDATE: 'data_list' has been overwritten with ECB-filtered data.")

# 6,5 . Modèle NSS OLS ECB guidelines compatible -------------------------------


# --- Objective Function (Weighted SSE) ---
nss_objective <- function(params, maturities, yields, weights) {
  # Constraints check (Soft boundaries via penalty)
  # b0 > 0, tau1 > 0, tau2 > 0
  if(params[1] < 0 || params[5] <= 0 || params[6] <= 0) {
    return(1e9) # Penalty for violating constraints
  }
  
  fitted <- nss_func(maturities, params)
  residuals <- yields - fitted
  
  # Weighted Sum of Squared Errors
  # Using Amt_Out as weight (Proxy for Liquidity preference mentioned in ECB paper)
  sse <- sum(weights * residuals^2, na.rm=TRUE)
  return(sse)
}

# --- D. Analysis Loop ---

plots_nss <- list()
scorecards_nss <- list()
params_nss <- list()

# Starting values (Heuristic: b0=long term, b1=spread, others=small)
# These are crucial for NSS as the function is non-convex
start_params <- c(b0=2, b1=-1, b2=0, b3=0, tau1=1, tau2=5) 

data_list <- split(df_all, df_all$Country)

for (country in names(data_list)) {
  
  # 1. Apply ECB Filters
  df_ecb <- data_list[[country]]
  
  # 2. Optimization
  # Using L-BFGS-B to enforce bounds: Tau > 0, B0 > 0
  opt_res <- optim(
    par = start_params,
    fn = nss_objective,
    maturities = df_ecb$Maturity,
    yields = df_ecb$Yield,
    weights = df_ecb$Amt_Out,
    method = "L-BFGS-B",
    lower = c(0, -15, -15, -15, 0.1, 0.1), # b0>0, tau>0.1
    upper = c(15, 15, 15, 15, 10, 30),
    control = list(maxit = 2000)
  )
  
  best_params <- opt_res$par
  params_nss[[country]] <- best_params
  
  # 3. Create Prediction Function closure
  predict_nss <- function(m) {
    nss_func(m, best_params)
  }
  
  # 4. Evaluation using your scorecard function
  # Note: We evaluate on the Filtered "High Quality" dataset, 
  # as ECB implies fitting errors should be measured on the selected sample.
  score <- evaluate_model_performance(
    model_name = paste("NSS_ECB_", country),
    data = df_ecb,
    predict_func = predict_nss,
    col_weight = "Amt_Out"
  )
  scorecards_nss[[country]] <- score
  
  # 5. Visualisation
  # On trace jusqu'à 30 ans (limite ECB)
  grid_x <- seq(0, 30, length.out = 300)
  grid_y <- predict_nss(grid_x)
  line_data <- data.frame(Maturity = grid_x, Yield = grid_y)
  
  # Graphique Principal : Courbe + Points (taille = volume)
  p <- ggplot(df_ecb, aes(x = Maturity, y = Yield)) +
    # Points observés (Nettoyés)
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.6) +
    # Courbe Modélisée
    geom_line(data = line_data, aes(x = Maturity, y = Yield), color = "black", size = 1) +
    
    scale_x_continuous(limits = c(0, 32), breaks = seq(0, 30, by = 5)) +
    labs(
      title = paste(country, "- NSS ECB (Cleaned Data)"),
      subtitle = paste("RMSE Pondéré:", score$W_RMSE, "| Stabilité (Wiggliness):", score$Wiggliness),
      y = "Yield (%)",
      size = "Volume (Mds)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", color = colors_flags[country]))
  
  # Graphique Secondaire : Résidus
  df_ecb$Residus <- df_ecb$Yield - predict_nss(df_ecb$Maturity)
  
  p_res <- ggplot(df_ecb, aes(x = Maturity, y = Residus)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.6) +
    geom_hline(yintercept = 0, color = "black") +
    scale_x_continuous(limits = c(0, 32), breaks = seq(0, 30, by = 5)) +
    labs(title = "Distribution des Résidus", y = "Ecart (bps)", size = "Volume") +
    theme_minimal()
  
  # Stockage dans la liste correspondante (plots_nss)
  plots_nss[[country]] <- grid.arrange(p, p_res, nrow = 2, heights = c(2, 1))
  
}

# --- Results ---
final_scorecard_nss <- bind_rows(scorecards_nss)
print(final_scorecard_nss)

# 7. Modèle NSS WLS ------------------------------------------------------------
# Hypothèse : On donne plus de poids aux obligations liquides (Gros Amt_Out).
# Cela permet d'ignorer les "Trash Bonds" qui faussent la courbe.

fit_nss_wls <- function(df, x_col = "Maturity") {
  weights <- df$Amt_Out
  weights[is.na(weights)] <- 0
  if(sum(weights) == 0) weights <- rep(1, length(weights))
  weights <- weights / sum(weights)
  
  # Détection de l'échelle pour le "sanity check"
  is_duration <- x_col == "Duration"
  check_point <- if(is_duration) 20 else 50 
  
  # Fonction Coût
  ssr_cost_wls <- function(params) {
    # Contrainte: tau doit être positif
    if(params[5] <= 0.1 || params[6] <= 0.5){
      return(1e9)
    }
    
    # Calcul des prédictions
    preds <- nss_func(df[[x_col]], params)
    # Erreur pondérée
    weighted_sq_err <- sum(weights * (df$Yield - preds)^2)
    
    # Pénalité de stabilité (Tail Check)
    # Si la courbe explose à 20 ans (Duration) ou 50 ans (Maturité), on pénalise
    tail_rate <- nss_func(check_point, params)
    penalty <- 0
    if(tail_rate < -2 || tail_rate > 15){
      penalty <- 1000
    }
    
    return(weighted_sq_err + penalty)
  }
  # Paramètres initiaux
  mean_yield <- mean(df$Yield, na.rm=TRUE)
  start_params <- c(b0=mean_yield, b1=0, b2=0, b3=0, tau1=1, tau2=5)
  
  # Optimisation
  opt <- optim(par = start_params, 
               fn = ssr_cost_wls, 
               method = "Nelder-Mead",
               control = list(maxit = 10000))
  
  return(opt$par)
}

plots_nsswls <- list()
scorecards_nsswls <- list()

for (country in names(data_list)) {
  df <- data_list[[country]]
  best_params <- fit_nss_wls(df, x_col = "Maturity")
  
  # Prediction
  predict_nss_wls <- function(m){
    nss_func(m, best_params)
    }
  
  score <- evaluate_model_performance(paste("NSS_WLS_", country), df, predict_nss_wls, "Amt_Out")
  scorecards_nsswls[[country]] <- score
  
  # Visualisation
  grid_x <- seq(0, 50, length.out = 500)
  line_data <- data.frame(Maturity = grid_x, Yield = predict_nss_wls(grid_x))
  
  p <- ggplot(df, aes(x = Maturity, y = Yield)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_line(data = line_data, color = "black", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = paste(country, "- NSS Pondéré par le Volume (WLS)"), 
         subtitle = "Fit Pondéré (Algorithme Nelder-Mead)", y="Yield") +
    theme_minimal() + theme(plot.title = element_text(face="bold", color=colors_flags[country]))
  
  df$Residus <- df$Yield - predict_nss_wls(df$Maturity)
  p_res <- ggplot(df, aes(x = Maturity, y = Residus)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Résidus", y="Bps", size="Volume") + theme_minimal()
  
  plots_nsswls[[country]] <- grid.arrange(p, p_res, nrow=2, heights=c(2,1))
}

# Affichage des Résultats
final_scorecard_nsswls <- bind_rows(scorecards_nsswls)
print(final_scorecard_nsswls)

# 8. Modèle NSS WLS en Duration ------------------------------------------------

# Fonction Macaulay (formule fermée)
calculate_duration <- function(yield_pct, coupon_pct, maturity) {
  y <- yield_pct / 100; c <- coupon_pct / 100; n <- maturity
  res <- numeric(length(n))
  for(i in seq_along(n)) {
    if(c[i] == 0) { res[i] <- n[i] } else {
      res[i] <- (1+y[i])/y[i] - (1+y[i]+n[i]*(c[i]-y[i])) / (c[i]*((1+y[i])^n[i]-1) + y[i])
    }
  }
  return(res)
}

plots_nsswlsduration <- list()
scorecards_nsswlsduration <- list()

for (country in names(data_list)) {
  df <- data_list[[country]]
  
  # 1. Calcul de la Duration
  df$Duration <- calculate_duration(df$Yield, df$Cpn, df$Maturity)
  
  # 2. Fitting (WLS sur DURATION)
  best_params_dur <- fit_nss_wls(df, x_col = "Duration")
  
  # 3. Prediction
  predict_nss_dur <- function(d) { nss_func(d, best_params_dur) }
  
  # Hack pour l'évaluation
  df_eval <- df
  df_eval$Maturity <- df$Duration 
  
  score <- evaluate_model_performance(paste("NSS_Duration_", country), df_eval, predict_nss_dur, "Amt_Out")
  scorecards_nsswlsduration[[country]] <- score
  
  # 4. Visu
  grid_d <- seq(0, 50, length.out = 500)
  line_data <- data.frame(Duration = grid_d, Yield = predict_nss_dur(grid_d))
  
  p <- ggplot(df, aes(x = Duration, y = Yield)) + 
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_line(data = line_data, color = "purple", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = paste(country, "- NSS sur Duration"),
         subtitle = "L'axe Duration comprime l'échelle temporelle (Nelder-Mead)",
         x = "Duration", y = "Yield") +
    theme_minimal() + theme(plot.title = element_text(face="bold", color=colors_flags[country]))
  
  df$Residus <- df$Yield - predict_nss_dur(df$Duration)
  p_res <- ggplot(df, aes(x = Duration, y = Residus)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Résidus vs Duration", y="Bps", x="Duration", size="Volume") + theme_minimal()
  
  plots_nsswlsduration[[country]] <- grid.arrange(p, p_res, nrow=2, heights=c(2,1))
}

# Affichage des Résultats
final_scorecard_nsswlsduration <- bind_rows(scorecards_nsswlsduration)
print(final_scorecard_nsswlsduration)

# 9. Modèle Splines cubiques ---------------------------------------------------

# Fonction d'optimisation du paramètre de lissage (spar) via validation croisée
fit_cubic_spline_optimized <- function(df) {
  # On teste une grille de paramètres de lissage
  # spar ~ 0 : Interpolation quasi-parfaite (très wiggly)
  # spar ~ 1 : Régression linéaire (très rigide)
  spar_grid <- seq(0.1, 1.5, by = 0.05) 
  
  best_score <- Inf
  best_model <- NULL

  weights <- df$Amt_Out
  weights[is.na(weights)] <- 0

  for(s in spar_grid) {
    # smooth.spline effectue une régression spline cubique pénalisée
    # cv = TRUE force la validation croisée généralisée (GCV) si dispo
    fit <- tryCatch({
      smooth.spline(x = df$Maturity, y = df$Yield, w = weights, spar = s, cv = TRUE)
    }, error = function(e) NULL)
    
    if(!is.null(fit)) {
      # On cherche à minimiser le critère de validation croisée (cv.crit)
      score <- fit$cv.crit
      if(score < best_score) {
        best_score <- score
        best_model <- fit
      }
    }
  }
  return(best_model)
}

plots_spline <- list()
scorecards_spline <- list()

print("--- Running Cubic Spline Analysis ---")

for (country in names(data_list)) {
  df <- data_list[[country]]
  
  # 1. Fitting
  best_spline <- fit_cubic_spline_optimized(df)
  
  # 2. Wrapper de prédiction
  # smooth.spline retourne une liste, on extrait $y
  predict_spline <- function(m) {
    predict(best_spline, x = m)$y
  }
  
  # 3. Évaluation
  score <- evaluate_model_performance(
    model_name = paste("Cubic_Spline_", country), 
    data = df, 
    predict_func = predict_spline, 
    col_weight = "Amt_Out"
  )
  scorecards_spline[[country]] <- score
  
  # 4. Visualisation
  grid_x <- seq(0, 50, length.out = 500)
  line_data <- data.frame(Maturity = grid_x, Yield = predict_spline(grid_x))
  
  # Graphe Principal : Fit
  p <- ggplot(df, aes(x = Maturity, y = Yield)) + 
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_line(data = line_data, color = "darkgreen", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(
      title = paste(country, "- Cubic Smoothing Spline"),
      subtitle = paste("Optimized Spar =", round(best_spline$spar, 3), "| CV Score =", round(best_spline$cv.crit, 4)),
      x = "Maturité (Années)", 
      y = "Yield (%)",
      size = "Volume"
    ) +
    theme_minimal() + 
    theme(plot.title = element_text(face = "bold", color = colors_flags[country]))
  
  # Graphe Secondaire : Résidus
  df$Residus <- df$Yield - predict_spline(df$Maturity)
  
  p_res <- ggplot(df, aes(x = Maturity, y = Residus)) +
    geom_point(aes(size = Amt_Out), color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black") +
    geom_smooth(method = "loess", se = FALSE, color = "gray", linetype = "dotted") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Structure des Résidus (Spline)", y = "Bps", x = "Maturité") +
    theme_minimal()
  
  # Stockage du plot combiné
  plots_spline[[country]] <- grid.arrange(p, p_res, nrow = 2, heights = c(2, 1))
}

# Affichage des Résultats
final_scorecard_spline <- bind_rows(scorecards_spline)
print(final_scorecard_spline)


































# ==============================================================================
# INTERLUDE : VÉRIFICATION DE L'HYPOTHÈSE (OLS PUR) - CORRIGÉ
# ==============================================================================
library(minpack.lm)
library(dplyr)
library(ggplot2)
library(gridExtra)

# 1. Fonction d'optimisation OLS PURE (Sans Poids)
fit_nss_ols_robust <- function(df, x_col = "Maturity") {
  
  y_obs <- df$Yield
  x_obs <- df[[x_col]] # C'est ici que ça plantait si x_col n'existe pas
  
  # Initialisation adaptée à l'échelle
  mean_y <- mean(y_obs, na.rm=TRUE)
  if (x_col == "Duration") {
    start_list <- list(b0 = mean_y+1, b1 = -1, b2 = 0, b3 = 0, tau1 = 1.0, tau2 = 3.0)
  } else {
    start_list <- list(b0 = mean_y+1, b1 = -1, b2 = 0, b3 = 0, tau1 = 2.0, tau2 = 10.0)
  }
  
  # Formule NSS explicite pour nlsLM
  nss_formula <- y_obs ~ b0 + 
    b1 * ((1 - exp(-x_obs/tau1)) / (x_obs/tau1)) + 
    b2 * (((1 - exp(-x_obs/tau1)) / (x_obs/tau1)) - exp(-x_obs/tau1)) + 
    b3 * (((1 - exp(-x_obs/tau2)) / (x_obs/tau2)) - exp(-x_obs/tau2))
  
  # Optimisation
  tryCatch({
    fit <- nlsLM(nss_formula, start = start_list,
                 lower = c(b0=0, b1=-30, b2=-100, b3=-100, tau1=0.1, tau2=0.2),
                 upper = c(b0=20, b1=30, b2=100, b3=100, tau1=20, tau2=50),
                 control = nls.lm.control(maxiter = 1000, ftol=1e-8))
    return(coef(fit))
  }, error = function(e) { return(NULL) })
}

# 2. Boucle de Comparaison
results_comparison <- list()
plots_comparison <- list()

print("--- Comparaison OLS : Maturité vs Duration ---")

for (country in names(data_list)) {
  df <- data_list[[country]]
  
  # On recalcule la Duration pour être sûr qu'elle existe
  # On utilise ta fonction calculate_duration (qui est maintenant la Modified Duration)
  df$Duration <- calculate_duration(df$Yield, df$Cpn, df$Maturity)
  
  # A. Modèle Maturité
  params_mat <- fit_nss_ols_robust(df, "Maturity")
  if(is.null(params_mat)) next # Sécurité si le fit échoue
  preds_mat <- nss_func(df$Maturity, params_mat)
  rmse_mat <- sqrt(mean((df$Yield - preds_mat)^2)) 
  
  # B. Modèle Duration 
  params_dur <- fit_nss_ols_robust(df, "Duration")
  if(is.null(params_dur)) next
  preds_dur <- nss_func(df$Duration, params_dur)
  rmse_dur <- sqrt(mean((df$Yield - preds_dur)^2))
  
  # Stockage Résultats
  results_comparison[[country]] <- tibble(
    Country = country,
    RMSE_Maturity_OLS = round(rmse_mat, 5),
    RMSE_Duration_OLS = round(rmse_dur, 5),
    Gain_Pct = round((rmse_mat - rmse_dur)/rmse_mat * 100, 2)
  )
  
  # C. Visualisation des Résidus
  df$Resid_Mat <- df$Yield - preds_mat
  df$Resid_Dur <- df$Yield - preds_dur
  
  p1 <- ggplot(df, aes(x = Maturity, y = Resid_Mat)) +
    geom_point(alpha=0.5, color="blue") + geom_hline(yintercept=0) +
    labs(title=paste(country, "Maturity OLS"), subtitle=paste("RMSE:", round(rmse_mat,4)), y="Residus") + theme_minimal()
  
  p2 <- ggplot(df, aes(x = Duration, y = Resid_Dur)) +
    geom_point(alpha=0.5, color="purple") + geom_hline(yintercept=0) +
    labs(title=paste(country, "Duration OLS"), subtitle=paste("RMSE:", round(rmse_dur,4)), y="Residus") + theme_minimal()
  
  plots_comparison[[country]] <- grid.arrange(p1, p2, nrow=1)
}

final_comparison <- bind_rows(results_comparison)
print(final_comparison)

# 8. Modèle NSS WLS avec pondération par la Duration ---------------------------


calculate_duration <- function(yield_pct, coupon_pct, maturity_years) {
  y <- yield_pct / 100
  c <- coupon_pct / 100
  
  n_bonds <- length(maturity_years)
  duration <- numeric(n_bonds)
  
  for (i in seq_len(n_bonds)) {
    # Nombre de paiements (arrondi standard)
    n <- ceiling(maturity_years[i])
    # Cas zéro-coupon
    if (c[i] == 0) {
      duration[i] <- maturity_years[i]
      next
    }
    
    # Flux de trésorerie
    times <- 1:n
    cashflows <- rep(c[i], n)
    cashflows[n] <- cashflows[n] + 1  # remboursement du nominal
    # Actualisation
    discount_factors <- 1 / (1 + y[i])^times
    pv_cashflows <- cashflows * discount_factors
    # Prix
    price <- sum(pv_cashflows)
    # Duration de Macaulay
    duration[i] <- sum(times * pv_cashflows) / price
  }
  return(duration)
}


plots_nsswls_durweight <- list()
scorecards_nsswls_durweight <- list()

for (country in names(data_list)) {
  df <- data_list[[country]]
  
  # 1. Calcul de la Duration (Macaulay)
  df$Duration <- calculate_duration(df$Yield, df$Cpn, df$Maturity)
  
  # 2. Poids économiques : Duration × Amount Outstanding
  df$W_Dur <- df$Duration * df$Amt_Out
  df$W_Dur[is.na(df$W_Dur)] <- 0
  df$W_Dur <- df$W_Dur / sum(df$W_Dur)  # normalisation
  
  # 3. Fonction coût NSS WLS (pondérée par Duration)
  ssr_cost_dur_wls <- function(params) {
    if(params[5] <= 0.1 || params[6] <= 0.5) return(1e9)
    
    preds <- nss_func(df$Maturity, params)
    sum(df$W_Dur * (df$Yield - preds)^2)
  }
  
  # 4. Estimation
  start_params <- c(
    b0 = mean(df$Yield),
    b1 = -1,
    b2 = 0,
    b3 = 0,
    tau1 = 1,
    tau2 = 5
  )
  
  opt <- optim(
    par = start_params,
    fn = ssr_cost_dur_wls,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  
  best_params <- opt$par
  
  # 5. Prédiction
  predict_nss <- function(m) nss_func(m, best_params)
  
  # 6. Évaluation (inchangée)
  score <- evaluate_model_performance(
    model_name = paste("NSS_WLS_DurationWeight_", country),
    data = df,
    predict_func = predict_nss,
    col_weight = "Amt_Out" #"W_Dur"
  )
  scorecards_nsswls_durweight[[country]] <- score
  
  # 7. Visualisation
  grid_x <- seq(0, 50, length.out = 500)
  line_data <- data.frame(Maturity = grid_x, Yield = predict_nss(grid_x))
  
  p <- ggplot(df, aes(x = Maturity, y = Yield)) +
    geom_point(aes(size = W_Dur), color = colors_flags[country], alpha = 0.5) +
    geom_line(data = line_data, color = "black", size = 1) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(
      title = paste(country, "- NSS WLS pondéré par la Duration"),
      subtitle = "Les obligations longues dominent la perte (logique DV01)",
      x = "Maturité",
      y = "Yield",
      size = "Poids éco"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", color = colors_flags[country]))
  
  df$Residus <- df$Yield - predict_nss(df$Maturity)
  p_res <- ggplot(df, aes(x = Maturity, y = Residus)) +
    geom_point(aes(size = W_Dur), color = colors_flags[country], alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    labs(title = "Résidus pondérés par la Duration", y = "Bps", size = "Poids éco") +
    theme_minimal()
  
  plots_nsswls_durweight[[country]] <- grid.arrange(p, p_res, nrow = 2, heights = c(2,1))
}

final_scorecard_nsswls_durweight <- bind_rows(scorecards_nsswls_durweight)
print(final_scorecard_nsswls_durweight)

# 4. Modeling (NSS & Splines) ---------------------------------------------

# Nelson-Siegel-Svensson Function
nss_func <- function(p, t) {
  if(any(is.na(p))) return(rep(NA, length(t)))
  term1 <- (1 - exp(-t/p[5])) / (t/p[5])
  term2 <- term1 - exp(-t/p[5])
  term3 <- (1 - exp(-t/p[6])) / (t/p[6]) - exp(-t/p[6])
  return(p[1] + p[2]*term1 + p[3]*term2 + p[4]*term3)
}

# Optimization Wrapper
fit_nss <- function(data) {
  obj <- function(p) {
    if(p[5] <= 0 || p[6] <= 0) return(1e9) # Constraints
    sum((data$Yield - nss_func(p, data$TTM))^2)
  }
  # Start: Beta0, Beta1, Beta2, Beta3, Tau1, Tau2
  start <- c(3, -1, -1, 0, 2, 5) 
  tryCatch(optim(par=start, fn=obj, method="BFGS")$par, error=function(e) rep(NA, 6))
}

# Create Grid covering max maturity
max_ttm <- max(df_all$TTM, na.rm=TRUE)
grid_ttm <- seq(0.1, ceiling(max_ttm) + 1, length.out = 200)

results_list <- list()

# Fit models for each country
for (ctry in unique(df_all$Country)) {
  sub_dat <- df_all %>% filter(Country == ctry)
  
  # 1. Spline
  fit_spl <- smooth.spline(sub_dat$TTM, sub_dat$Yield, spar=0.5)
  y_spl <- predict(fit_spl, grid_ttm)$y
  
  # 2. NSS
  p_nss <- fit_nss(sub_dat)
  y_nss <- nss_func(p_nss, grid_ttm)
  
  # On calcule la valeur théorique pour les points existants
  y_fitted <- nss_func(p_nss, sub_dat$TTM) 
  sub_dat$Residual <- sub_dat$Yield - y_fitted

  # Store results
  results_list[[ctry]] <- list(
    params = p_nss,
    preds = data.frame(TTM = grid_ttm, Spline = y_spl, NSS = y_nss, Country = ctry),
    data = sub_dat
  )
}

# 5. Plotting: One Graph Per Country --------------------------------------

# Helper function to generate standardized plots
plot_country_curve <- function(country_name) {
  
  # Get observed data
  obs_data <- df_all %>% filter(Country == country_name)
  
  # Get fitted data and reshape for ggplot
  fit_data <- results_list[[country_name]]$preds %>%
    pivot_longer(cols = c("Spline", "NSS"), names_to = "Model", values_to = "Yield")
  
  # Plot
  p <- ggplot() +
    # 1. Observed Points
    geom_point(data = obs_data, aes(x = TTM, y = Yield), alpha = 0.6, size = 1.5) +
    # 2. Fitted Lines (Spline + NSS)
    geom_line(data = fit_data, aes(x = TTM, y = Yield, color = Model, linetype = Model), size = 1) +
    # Aesthetics
    scale_color_manual(values = c("NSS" = "blue", "Spline" = "red")) +
    scale_linetype_manual(values = c("NSS" = "solid", "Spline" = "dashed")) +
    labs(title = paste(country_name, "- Yield Curve Analysis"),
         subtitle = "Comparison: Observed Points vs. Spline vs. NSS",
         y = "Yield (%)", x = "Maturity (Years)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Generate and Print Plots
p1 <- plot_country_curve("France")
p2 <- plot_country_curve("Germany")
p3 <- plot_country_curve("Italy")

print(p1)
print(p2)
print(p3)

# 6. Spread Analysis ------------------------------------------

tenors <- seq(1, floor(max_ttm), 0.5)
y_fr <- nss_func(results_list[["France"]]$params, tenors)
y_de <- nss_func(results_list[["Germany"]]$params, tenors)
y_it <- nss_func(results_list[["Italy"]]$params, tenors)

spread_df <- data.frame(Maturity = tenors,
                        France_Spread = (y_fr - y_de) * 100,
                        Italy_Spread = (y_it - y_de) * 100) %>%
  pivot_longer(-Maturity, names_to = "Pair", values_to = "Spread_bps")

p_spread <- ggplot(spread_df, aes(x = Maturity, y = Spread_bps, color = Pair)) +
  geom_line(size = 1.2) +
  labs(title = "Sovereign Spreads vs Germany (NSS Model)",
       y = "Spread (bps)", x = "Maturity (Years)") +
  theme_minimal()

print(p_spread)

# 7. Residuals Analysis ------------------------------------------

# Fonction pour tracer les résidus et identifier les outliers
plot_residuals_analysis <- function(country_name) {
  
  # Récupération des données stockées
  dat <- results_list[[country_name]]$data
  
  # Calcul des seuils pour outliers (2 écarts-types)
  sigma <- sd(dat$Residual, na.rm = TRUE)
  threshold <- 2 * sigma
  
  # Identifier les outliers
  dat <- dat %>%
    mutate(IsOutlier = abs(Residual) > threshold,
           Label = ifelse(IsOutlier, as.character(MaturityDate), NA))
  
  # Plot
  p <- ggplot(dat, aes(x = TTM, y = Residual)) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    # Zone de tolérance normale (grise)
    geom_ribbon(aes(ymin = -threshold, ymax = threshold), alpha = 0.1, fill = "blue") +
    # Points (Rouge si outlier, Bleu sinon)
    geom_point(aes(color = IsOutlier), size = 2, alpha = 0.8) +
    # Labels pour les outliers seulement
    geom_text(aes(label = Label), vjust = -0.8, color = "red", size = 3, check_overlap = TRUE) +
    # Esthétique
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
    labs(title = paste("Analyse des Résidus NSS -", country_name),
         subtitle = paste("Outliers détectés (> 2 sigma). RMSE =", round(sqrt(mean(dat$Residual^2)), 4)),
         y = "Résidu (Yield Observé - NSS)",
         x = "Maturité (Années)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}

grid.arrange(
  plot_residuals_analysis("France"),
  plot_residuals_analysis("Germany"),
  plot_residuals_analysis("Italy"),
  ncol = 1
)
# 8. NSS coefficients Analysis ------------------------------------------

# Fonction pour afficher le tableau des paramètres
get_nss_parameters_table <- function() {
  # Extraction
  params_df <- do.call(rbind, lapply(names(results_list), function(nm) {
    p <- results_list[[nm]]$params
    data.frame(
      Country = nm,
      Beta0_LongTerm = p[1],
      Beta1_Slope = p[2],
      Beta2_Hump = p[3],
      Beta3_Hump2 = p[4],
      Tau1 = p[5],
      Tau2 = p[6]
    )
  }))
  
  return(params_df)
}

print("Paramètres du modèle Nelson-Siegel-Svensson :")
print(get_nss_parameters_table())


# 9. Analyse de la courbe swap et comparaison -------------------------------------

# --- 9.1 Nettoyage des données Swaps ---

clean_swap_data <- function(df) {
  # Fonction pour convertir "6M", "2Y" en années numériques
  convert_tenor <- function(t) {
    if (grepl("M", t)) {
      return(as.numeric(gsub("M", "", t)) / 12)
    } else if (grepl("Y", t)) {
      return(as.numeric(gsub("Y", "", t)))
    } else {
      return(NA)
    }
  }
  
  df %>%
    # Adapter les noms de colonnes ici si nécessaire
    # On suppose ici que la colonne des taux s'appelle 'Rate' ou est la 2ème colonne
    rename(TenorStr = 1, SwapRate = 2) %>% 
    mutate(
      TTM = sapply(TenorStr, convert_tenor),
      Yield = as.numeric(SwapRate), # Si les données sont déjà en %, sinon * 100
      Country = "Swap (Risk Free)"
    ) %>%
    filter(!is.na(TTM), !is.na(Yield)) %>%
    arrange(TTM)
}

df_swaps <- clean_swap_data(data_swaps)

# --- 9.2 Fitting du modèle NSS sur la courbe des Swaps ---
# Les Swaps étant des données "propres", le modèle NSS doit fitter parfaitement (Calibration)

p_swap <- fit_nss(df_swaps)
y_swap_nss <- nss_func(p_swap, grid_ttm)

# On stocke pour l'analyse
results_list[["Swap"]] <- list(
  params = p_swap,
  preds = data.frame(TTM = grid_ttm, NSS = y_swap_nss, Country = "Swap"),
  data = df_swaps
)

message("Paramètres NSS Swap calculés.")

# --- 9.3 Analyse 1 : Superposition des Courbes Souveraines et Swap ---

# On récupère les courbes NSS générées précédemment
plot_data_comparison <- bind_rows(
  results_list[["Germany"]]$preds %>% select(TTM, NSS) %>% mutate(Curve = "Allemagne"),
  results_list[["Italy"]]$preds %>% select(TTM, NSS) %>% mutate(Curve = "Italie"),
  results_list[["France"]]$preds %>% select(TTM, NSS) %>% mutate(Curve = "France"),
  results_list[["Swap"]]$preds %>% select(TTM, NSS) %>% mutate(Curve = "Swap")
)

p_superposition <- ggplot(plot_data_comparison, aes(x = TTM, y = NSS, color = Curve)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Allemagne" = "darkgreen", 
                                "Swap" = "black", 
                                "Italie" = "red",
                                "France" = "blue")) +
  labs(title = "Hiérarchie des Taux en Zone Euro",
       y = "Yield (%)", x = "Maturité (Années)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_superposition)

# --- 9.4 Analyse 2 : Asset Swap Spreads (ASW) ---
# Objectif : Calculer l'écart exact entre chaque obligation et le taux swap interpolé

calculate_asw <- function(country_name) {
  bond_data <- results_list[[country_name]]$data
  
  # Calcul du taux swap théorique pour chaque maturité d'obligation
  # On utilise le modèle NSS calibré sur les swaps pour interpoler
  swap_yields_at_bond_maturities <- nss_func(p_swap, bond_data$TTM)
  
  bond_data %>%
    mutate(
      SwapYield = swap_yields_at_bond_maturities,
      ASW_Spread = (Yield - SwapYield) * 100, # En points de base (bps)
      Country = country_name
    )
}

asw_germany <- calculate_asw("Germany")
asw_italy <- calculate_asw("Italy")
asw_france <- calculate_asw("France")
asw_data <- bind_rows(asw_germany, asw_italy, asw_france)

p_asw <- ggplot(asw_data, aes(x = TTM, y = ASW_Spread, color = Country)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, size = 0.8, linetype="dashed") + # Tendance visuelle par LOESS
  labs(title = "Asset Swap Spreads (ASW)",
       y = "Spread vs Swap (bps)", x = "Maturité (Années)") +
  scale_color_manual(values = c("Germany" = "darkgreen", "Italy" = "red", "France" = "blue")) +
  theme_minimal()

print(p_asw)

# --- 9.5 Analyse 3 : Analyse de Pente (2Y - 10Y) ---
# Objectif : Voir si la courbe italienne est plus pentue (=prime de terme excessive) que le Swap

calc_slope_2_10 <- function(params) {
  y2 <- nss_func(params, 2)
  y10 <- nss_func(params, 10)
  return((y10 - y2) * 100) # Pente en bps
}

slope_swap <- calc_slope_2_10(p_swap)
slope_italy <- calc_slope_2_10(results_list[["Italy"]]$params)
slope_germany <- calc_slope_2_10(results_list[["Germany"]]$params)
slope_france <- calc_slope_2_10(results_list[["France"]]$params)

slope_df <- data.frame(
  Curve = c("Swap", "Italy", "Germany", "France"),
  Slope_2s10s = c(slope_swap, slope_italy, slope_germany, slope_france)
)

p_slope <- ggplot(slope_df, aes(x = reorder(Curve, Slope_2s10s), y = Slope_2s10s, fill = Curve)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(Slope_2s10s, 1)), vjust = -0.5) +
  scale_fill_manual(values = c("Germany" = "darkgreen", "Swap" = "gray50", "Italy" = "red", "France" = "blue")) +
  labs(title = "Pente de la courbe (2 ans - 10 ans)",
       subtitle = "Une pente plus forte indique une prime de risque à long terme plus élevée",
       y = "Pente (bps)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_slope)