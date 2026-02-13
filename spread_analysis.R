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
nss_wls_results <- readRDS("models/nss_wls_results.rds")

grid_ttm <- seq(1, 30, length.out = 100)

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

# --- Calcul des courbes de taux (Yields) ---
# On utilise les paramètres stockés dans nss_wls_results
y_fr <- nss_func(grid_ttm, nss_wls_results[["France"]]$params)
y_de <- nss_func(grid_ttm, nss_wls_results[["Germany"]]$params)
y_it <- nss_func(grid_ttm, nss_wls_results[["Italy"]]$params)

# --- Calcul des Spreads vs Allemagne ---
# Le spread est la différence de rendement. 
df_spread <- data.frame(
  Maturity = grid_ttm,
  France = (y_fr - y_de), 
  Italy  = (y_it - y_de)
) %>%
  pivot_longer(cols = c("France", "Italy"), 
               names_to = "Country", 
               values_to = "Spread")

# --- Visualisation ---
p_spread <- ggplot(df_spread, aes(x = Maturity, y = Spread, color = Country)) +
  # Ligne de référence à 0 (le niveau de l'Allemagne)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_flags) +
  
  labs(title = "Sovereign Spreads vs Germany (NSS Model)",
       subtitle = "Ecart de rendement modélisé par rapport au Bund allemand",
       y = "Spread (bps)",
       x = "Maturité (Années)") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_spread)




# ==============================================================================
# 9. ANALYSE DE LA COURBE SWAP ET COMPARAISON (ASW / SPREADS)
# ==============================================================================

results_nsswls <- readRDS("models/nss_wls_results.rds")

# --- 9.1 Fonction de Fitting Spécifique aux Swaps (Calibration) ---
# Calibration parfaite requise pour la courbe Risk-Free

fit_nss_swap <- function(df) {
  
  maturities <- df$TTM
  yields <- df$Yield
  weights <- rep(1, length(yields)) 
  
  # Paramètres initiaux
  start_params <- c(b0=mean(yields), b1=-1, b2=0, b3=0, tau1=1, tau2=5)
  
  # Optimisation
  opt_res <- optim(
    par = start_params,
    fn = nss_objective,
    maturities = maturities,
    yields = yields,
    weights = weights,
    method = "L-BFGS-B",
    lower = c(0, -15, -15, -15, 0.1, 0.1),
    upper = c(20, 15, 15, 15, 10, 30),
    control = list(maxit = 5000, factr = 1e7)
  )
  
  message("Calibration Swap terminée. SSE: ", round(opt_res$value, 6))
  return(opt_res$par)
}

# --- 9.2 Nettoyage des données Swaps (CORRIGÉ) ---

clean_swap_data <- function(df) {
  convert_tenor <- function(t) {
    if (grepl("M", t)) {
      return(as.numeric(gsub("M", "", t)) / 12)
    } else if (grepl("Y", t)) {
      return(as.numeric(gsub("Y", "", t)))
    } else {
      return(NA)
    }
  }
  
  df_clean <- df %>%
    # CORRECTION ICI : Le Yield est dans la colonne 3, pas la 2 !
    select(TenorStr = 1, SwapRate = 3) %>% 
    mutate(
      TTM = sapply(TenorStr, convert_tenor),
      Yield = as.numeric(as.character(SwapRate)), # Conversion directe
      Country = "Swap (Risk Free)"
    ) %>%
    filter(!is.na(TTM), !is.na(Yield)) %>%
    arrange(TTM)
  
  # --- AUTO-SCALING (Sécurité) ---
  # Vérification que les données sont bien en % (ex: 3.0) et pas en décimales (0.03)
  avg_yield <- mean(df_clean$Yield, na.rm = TRUE)
  
  if (avg_yield < 0.5) {
    message("⚠️ Correction d'échelle Swap : Conversion Décimale -> Pourcentage (*100)")
    df_clean$Yield <- df_clean$Yield * 100
  } 
  # Note: On a supprimé la division par 100 car vos données brutes (2.11, 3.12) sont déjà correctes.
  
  return(df_clean)
}

df_swaps <- clean_swap_data(data_swaps)

# --- 9.3 Calibration du modèle ---

p_swap <- fit_nss_swap(df_swaps)
grid_ttm <- seq(0, 30, length.out = 100)

# Stockage
results_nsswls[["Swap"]] <- list(
  params = p_swap,
  data = df_swaps,
  preds = data.frame(TTM = grid_ttm, 
                     Yield = nss_func(grid_ttm, p_swap), 
                     Curve = "Swap")
)

# --- 9.4 Visualisation : Hiérarchie des Taux ---

get_curve_points <- function(country, label) {
  params <- results_nsswls[[country]]$params 
  # Rappel: Ordre (Time, Params)
  y <- nss_func(grid_ttm, params)
  data.frame(TTM = grid_ttm, Yield = y, Curve = label)
}

# On s'assure que la colonne Yield du Swap a le bon nom pour le bind_rows
swap_preds <- results_nsswls[["Swap"]]$preds %>% 
  rename(Yield = 2) # S'assure que la colonne 2 s'appelle 'Yield'

plot_data_comparison <- bind_rows(
  get_curve_points("Germany", "Allemagne"),
  get_curve_points("Italy", "Italie"),
  get_curve_points("France", "France"),
  swap_preds
)

p_superposition <- ggplot(plot_data_comparison, aes(x = TTM, y = Yield, color = Curve)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Allemagne" = "#BA0C2F", 
                                "Swap" = "black", 
                                "Italie" = "#008C45",
                                "France" = "#0055A4")) +
  labs(title = "Hiérarchie des Courbes de Taux (NSS WLS)",
       y = "Yield (%)", x = "Maturité (Années)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_superposition)

# --- 9.5 Analyse Asset Swap Spreads (ASW) ---

# Calcul de l'écart entre chaque obligation réelle et le taux swap interpolé
calculate_asw <- function(country_name) {
  # On reprend les données brutes filtrées ECB utilisées pour le fit
  bond_data <- results_nsswls[[country_name]]$data
  
  # Quel est le taux swap théorique pour la maturité exacte de chaque obligation ?
  swap_yields <- nss_func(bond_data$Maturity, p_swap)
  
  bond_data %>%
    mutate(
      SwapYield = swap_yields,
      ASW_Spread = Yield - SwapYield,
      Country = country_name
    )
}

asw_data <- bind_rows(
  calculate_asw("Germany"),
  calculate_asw("France"),
  calculate_asw("Italy")
)

p_asw <- ggplot(asw_data, aes(x = Maturity, y = ASW_Spread, color = Country)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Tendance lissée (LOESS) pour voir la structure par terme du spread
  geom_smooth(method = "loess", se = FALSE, size = 1, linetype="solid") +
  
  scale_color_manual(values = colors_flags) +
  scale_y_continuous(n.breaks = 10) + 
  labs(title = "Structure par Terme des Asset Swap Spreads (ASW)",
       y = "Spread vs Swap (bps)", x = "Maturité (Années)", size = "Volume") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_asw)

# --- 9.6 Analyse de Pente (2Y - 10Y) ---

calc_slope_2_10 <- function(params) {
  y2 <- nss_func(2, params)
  y10 <- nss_func(10, params)
  return((y10 - y2) * 100)
}

slope_df <- data.frame(
  Curve = c("Swap", "Italy", "Germany", "France"),
  Slope_2s10s = c(
    calc_slope_2_10(p_swap),
    calc_slope_2_10(results_nsswls[["Italy"]]$params),
    calc_slope_2_10(results_nsswls[["Germany"]]$params),
    calc_slope_2_10(results_nsswls[["France"]]$params)
  )
)

p_slope <- ggplot(slope_df, aes(x = reorder(Curve, Slope_2s10s), y = Slope_2s10s, fill = Curve)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = round(Slope_2s10s, 0)), vjust = -0.5, fontface="bold") +
  scale_fill_manual(values = c("Germany" = "#BA0C2F", 
                               "Swap" = "gray40", 
                               "Italy" = "#008C45", 
                               "France" = "#0055A4")) +
  labs(title = "Pente de la Courbe (2s10s)",
       y = "Pente (bps)", x = "") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_slope)