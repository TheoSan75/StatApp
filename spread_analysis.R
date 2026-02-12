rm(list = ls())

# 1. Setup ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(gridExtra)
source("utils.R")
source("src/model_nss.R")

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
data_list <- split(df_all, df_all$Country)
data_list_clean <- lapply(data_list, filter_ecb_criteria)




tenors <- seq(1, 50, 0.5)
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