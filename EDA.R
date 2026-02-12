rm(list = ls())

# 1. Setup ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(gridExtra)
source("utils.R")

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

# 3. Manipulation de la date ---------------------------------------------------
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

results_lin <- list() # Conteneur pour stocker les modèles
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
  # On stocke le modèle brut (lm) et le score
  results_lin[[country]] <- list(
    model = model_lin,
    score = score,
    data = df # Optionnel : si vous voulez garder les données associées
  )
}
saveRDS(results_lin, file = "models/linear_results.rds")
# Affichage du tableau récapitulatif
final_scorecard_reglin <- bind_rows(scorecards_reglin)
print(final_scorecard_reglin)


# 6. Modèle NSS OLS ------------------------------------------------------------

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
results_nssols <- list()

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
  results_nssols[[country]] <- list(
    params = best_params,   # Les 6 paramètres NSS (b0, b1...)
    score = score,          # La scorecard
    data = df               # Les données utilisées
  )
}
saveRDS(results_nssols, file = "models/nss_ols_results.rds")
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
data_list_clean <- lapply(data_list, filter_ecb_criteria)

# 2. Visualize with Dynamic Buckets (Green Tunnel Style)
plot_cleaning_impact <- function(country) {
  
  # A. Prepare Data
  raw <- data_list[[country]] %>% mutate(Status = "Removed")
  clean <- data_list_clean[[country]] %>% mutate(Status = "Kept")
  
  removed_points <- anti_join(raw, clean, by = c("Maturity", "Yield"))
  combined_points <- bind_rows(clean, removed_points)
  
  # B. Create Dynamic Buckets
  n_obs_total <- nrow(clean)
  target_n <- 10 
  n_buckets <- max(1, floor(n_obs_total / target_n))
  
  stats_per_bucket <- clean %>%
    mutate(Bucket = ggplot2::cut_number(Maturity, n = n_buckets)) %>%
    group_by(Bucket) %>%
    summarise(
      mu = mean(Yield, na.rm=TRUE), 
      sigma = sd(Yield, na.rm=TRUE),
      raw_interval = unique(as.character(Bucket)),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mu - 2*sigma, 
      ymax = mu + 2*sigma,
      xmin = as.numeric(sub("^\\((.*),.*", "\\1", raw_interval)),
      xmin = ifelse(is.na(xmin), as.numeric(sub("^\\[(.*),.*", "\\1", raw_interval)), xmin),
      xmax = as.numeric(sub(".*,(.*)\\]", "\\1", raw_interval))
    ) %>%
    arrange(xmin)
  
  # C. Plotting
  p <- ggplot() +
    
    # --- 0. BACKGROUND FILL (Green Area) --- 
    # Must be first to be behind the points
    geom_rect(data = stats_per_bucket,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "darkgreen", alpha = 0.15) + # alpha gère la transparence
    
    # --- 1. The Bounds (Step Lines - Deep Green Dashed) ---
    # Upper Bound
    geom_segment(data = stats_per_bucket, 
                 aes(x = xmin, xend = xmax, y = ymax, yend = ymax), 
                 color = "darkgreen", linetype = "dashed", size = 0.6) +
    # Lower Bound
    geom_segment(data = stats_per_bucket, 
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymin), 
                 color = "darkgreen", linetype = "dashed", size = 0.6) +
    
    # --- 2. Vertical Connectors (Deep Green Dashed) ---
    geom_segment(data = stats_per_bucket[-nrow(stats_per_bucket),],
                 aes(x = xmax, xend = xmax, y = ymax, yend = stats_per_bucket$ymax[-1]),
                 color = "darkgreen", linetype = "dashed", size = 0.4) +
    geom_segment(data = stats_per_bucket[-nrow(stats_per_bucket),],
                 aes(x = xmax, xend = xmax, y = ymin, yend = stats_per_bucket$ymin[-1]),
                 color = "darkgreen", linetype = "dashed", size = 0.4) +
    
    # --- 3. Mean Line (Black Dotted - kept distinct for readability) ---
    geom_segment(data = stats_per_bucket, 
                 aes(x = xmin, xend = xmax, y = mu, yend = mu), 
                 color = "black", linetype = "dotted", size = 0.5) +
    
    # --- 4. Maturity Limits ---
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "darkred", size = 0.5) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "darkred", size = 0.5) +
    
    # --- 5. Data Points ---
    geom_point(data = combined_points, 
               aes(x = Maturity, y = Yield, color = Status, shape = Status), 
               size = 2, alpha = 0.7) +
    
    # --- 6. Styling ---
    scale_color_manual(values = c("Kept" = "#0055A4", "Removed" = "red")) +
    scale_shape_manual(values = c("Kept" = 19, "Removed" = 4)) + 
    scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
    
    labs(
      title = paste("Outlier Detection (Dynamic Buckets) -", country),
      subtitle = paste0("Green Area = +/- 2 Sigma (", n_buckets, " buckets) | Red Cross = Outlier"),
      x = "Residual Maturity", y = "Yield"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", color = colors_flags[country])
    )
  
  return(p)
}

# Generate plots
p_clean_fr <- plot_cleaning_impact("France")
p_clean_de <- plot_cleaning_impact("Germany")
p_clean_it <- plot_cleaning_impact("Italy")

# Display
print(p_clean_fr)
print(p_clean_de)
print(p_clean_it)

# 3. OVERWRITE GLOBAL DATA
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
results_nss_ecb <- list()

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
  results_nss_ecb[[country]] <- list(
    params = best_params,   # Les paramètres optimisés
    score = score,          # La performance
    data = df_ecb           # Les données nettoyées (IMPORTANT pour reproduire le fit)
  )
}
saveRDS(results_nss_ecb, file = "models/nss_ecb_results.rds")
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
results_nsswls <- list()

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
  results_nsswls[[country]] <- list(
    params = best_params,   # Les paramètres optimisés (WLS)
    score = score,          # La performance
    data = df               # Les données utilisées
  )
}
saveRDS(results_nsswls, file = "models/nss_wls_results.rds")
# Affichage des Résultats
final_scorecard_nsswls <- bind_rows(scorecards_nsswls)
print(final_scorecard_nsswls)

# 8. Modèle NSS WLS en Duration ------------------------------------------------

plots_nsswlsduration <- list()
scorecards_nsswlsduration <- list()
results_nss_duration <- list()

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
  results_nss_duration[[country]] <- list(
    params = best_params_dur, # Les paramètres optimisés sur la Duration
    score = score,            # La performance
    data = df                 # Les données (avec la colonne Duration incluse)
  )
}
saveRDS(results_nss_duration, file = "models/nss_duration_results.rds")
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
results_spline <- list()

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
  results_spline[[country]] <- list(
    model = best_spline,    # L'objet spline complet (nécessaire pour predict)
    score = score,          # La performance
    data = df               # Les données utilisées
  )
}
saveRDS(results_spline, file = "models/cubic_spline_results.rds")
# Affichage des Résultats
final_scorecard_spline <- bind_rows(scorecards_spline)
print(final_scorecard_spline)





# 10. SYNTHESIS & FINAL COMPARISON ---------------------------------------------

# 1. Aggregate all individual scorecards
# We use bind_rows to stack them vertically
global_scorecard <- bind_rows(
  bind_rows(scorecards_reglin) %>% mutate(Method = "Linear"),
  bind_rows(scorecards_nssols) %>% mutate(Method = "NSS OLS (Raw)"),
  bind_rows(scorecards_nss)    %>% mutate(Method = "NSS ECB"),
  bind_rows(scorecards_nsswls) %>% mutate(Method = "NSS WLS"),
  bind_rows(scorecards_nsswlsduration) %>% mutate(Method = "NSS Duration"),
  bind_rows(scorecards_spline) %>% mutate(Method = "Cubic Spline")
)

# 2. Reorder columns for readability
# We put 'Country' and 'Method' first, then the metrics
global_scorecard <- global_scorecard %>%
  # Extract Country from the Model Name if needed, or rely on row order
  # Assuming your 'Model' column looks like "Linear_France", we can clean it:
  separate(Model, into = c("Model_Type", "Country"), sep = "_", extra = "merge") %>%
  select(Country, Method, W_RMSE, Price_RMSE, Wiggliness, Tail_Stab_30_50, starts_with("Bias")) %>%
  arrange(Country, W_RMSE) # Sort by Country and then best performance (lowest RMSE)

# 3. Display the Final Leaderboard
print("--- GLOBAL MODEL LEADERBOARD ---")
print(global_scorecard)