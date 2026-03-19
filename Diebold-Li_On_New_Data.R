rm(list = ls())
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(viridis)
library(tseries)
source("utils.R")

# --- 1. IMPORTATION ---
data_raw <- read_csv("r_data/NewData_Souverain_CrossTempo_Long_Clean.csv",
                     show_col_types = FALSE) |>
  mutate(Date = as.Date(Date), Yield = as.numeric(Yield)) |>
  filter(!is.na(Yield), Residual_Maturity > 0)

# --- 2. SÉLECTION DES 5 PAYS AVEC LE PLUS DE BONDS ---
top5_countries <- data_raw |>
  group_by(Country) |>
  summarise(N_bonds = n_distinct(ISIN), .groups = "drop") |>
  slice_max(N_bonds, n = 5) |>
  pull(Country)

message("Top 5 pays : ", paste(top5_countries, collapse = ", "))

# --- 3. FORMATAGE (identique à l'original) ---
df_final <- data_raw |>
  filter(Country %in% top5_countries) |>
  rename(Time_To_Maturity = Residual_Maturity) |>
  select(Date, Country, Time_To_Maturity, Yield) |>
  arrange(Date, Country, Time_To_Maturity)

# --- PHASE 1 : NETTOYAGE ECB ET CONSTRUCTION DE LA GRILLE ---

df_cleaned <- df_final |>
  filter(Time_To_Maturity >= 1,
         Time_To_Maturity <= 30)

# ── ECB PRE-PROCESSING ────────────────────────────────────────────────────────
# filter_ecb_criteria() attend une colonne "Maturity" et une colonne "Yield"
df_cleaned_ecb <- df_cleaned |>
  rename(Maturity = Time_To_Maturity) |>
  group_by(Date, Country) |>
  group_modify(~ filter_ecb_criteria(.x)) |>
  ungroup() |>
  rename(Time_To_Maturity = Maturity)

# Diagnostic : combien de points supprimés ?
n_before <- nrow(df_cleaned)
n_after  <- nrow(df_cleaned_ecb)
message("ECB filtering : ", n_before, " → ", n_after, " obs (",
        round(100 * (n_before - n_after) / n_before, 1), "% supprimés)")

# On remplace df_cleaned par la version nettoyée
df_cleaned <- df_cleaned_ecb
# ─────────────────────────────────────────────────────────────────────────────

target_maturities <- c(
  0.25, 0.5, 0.75, 1,
  1.25, 1.5, 1.75, 2,
  2.5, 3, 4, 5,
  6, 7, 8, 9, 10,
  15, 20, 30
)

interpolate_yields <- function(sub_df, targets) {
  sub_df <- sub_df |> arrange(Time_To_Maturity)
  if (nrow(sub_df) < 2) return(data.frame(Maturity = targets, Yield = NA))
  approx_res <- approx(
    x      = sub_df$Time_To_Maturity,
    y      = sub_df$Yield,
    xout   = targets,
    method = "linear",
    rule   = 2
  )
  return(data.frame(Maturity = approx_res$x, Yield = approx_res$y))
}

df_interpolated <- df_cleaned |>
  group_by(Date, Country) |>
  do(interpolate_yields(., target_maturities)) |>
  ungroup()

df_matrix_form <- df_interpolated |>
  pivot_wider(names_from = Maturity, values_from = Yield, names_prefix = "Mat_")

# --- PHASE 2 : CALIBRAGE DU LAMBDA ---
# Courbe moyenne sur les 5 pays combinés (identique à l'original)
df_avg_curve <- df_interpolated |>
  group_by(Maturity) |>
  summarise(Avg_Yield = mean(Yield, na.rm = TRUE))

# Remplacement du lambda unique par un lambda par pays
ns_rmse_function <- function(lambda, data) {
  tau <- data$Maturity
  x2  <- (1 - exp(-lambda * tau)) / (lambda * tau)
  x3  <- x2 - exp(-lambda * tau)
  model <- lm(data$Avg_Yield ~ x2 + x3)
  return(sqrt(mean(residuals(model)^2)))
}
# Puis courbe pour chaque pays (plus robuste)
df_avg_by_country <- df_interpolated |>
  group_by(Country, Maturity) |>
  summarise(Avg_Yield = mean(Yield, na.rm = TRUE), .groups = "drop")

lambda_by_country <- df_avg_by_country |>
  group_by(Country) |>
  summarise(
    Lambda = optimize(
      ns_rmse_function,
      interval = c(0.01, 3.0),
      data     = data.frame(Maturity  = Maturity,
                            Avg_Yield = Avg_Yield)
    )$minimum,
    .groups = "drop"
  )

message("\n── Lambda optimal par pays ──")
print(lambda_by_country)

# Lambda moyen pour les graphes globaux (factor loadings, etc.)
best_lambda <- mean(lambda_by_country$Lambda)
cat("Lambda moyen (référence) :", round(best_lambda, 4), "\n")

cat("Lambda optimal :", round(best_lambda, 4), "\n")
cat("Pic de courbure à :", round(1.7932 / best_lambda, 2), "ans\n")

# Factor loadings plot (Germany only)
plot_mats <- seq(0, 30, length.out = 100)
loadings_df <- data.frame(
  Maturity  = plot_mats,
  Level     = 1,
  Slope     = (1 - exp(-best_lambda * plot_mats)) / (best_lambda * plot_mats),
  Curvature = ((1 - exp(-best_lambda * plot_mats)) / (best_lambda * plot_mats)) - exp(-best_lambda * plot_mats)
) |>
  pivot_longer(cols = c("Level", "Slope", "Curvature"), names_to = "Factor", values_to = "Loading")

p_loadings <- ggplot(loadings_df, aes(x = Maturity, y = Loading, color = Factor, linetype = Factor)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 1.7932 / best_lambda, linetype = "dotted", color = "black") +
  annotate("text", x = (1.7932 / best_lambda) + 1, y = 0.1,
           label = paste0("Max Courbure\n~", round(1.7932 / best_lambda, 2), " ans"),
           hjust = 0, size = 3.5) +
  scale_color_manual(values = c("Curvature" = "#D55E00", "Level" = "#009E73", "Slope" = "#0072B2")) +
  labs(title = "Factor Loadings — Nelson-Siegel",
       subtitle = paste0("Lambda optimisé = ", round(best_lambda, 3)),
       x = "Maturité (Années)", y = "Poids (Loading)") +
  theme_minimal()
print(p_loadings)

# --- PHASE 3 : CALCUL DES SÉRIES TEMPORELLES DES BETAS ---

# Estimation betas avec lambda propre à chaque pays
estimate_betas <- function(sub_df, lambda) {
  tau <- sub_df$Maturity
  x2  <- (1 - exp(-lambda * tau)) / (lambda * tau)
  x3  <- x2 - exp(-lambda * tau)
  model  <- lm(sub_df$Yield ~ x2 + x3)
  coeffs <- coef(model)
  return(data.frame(
    Beta1_Level     = coeffs["(Intercept)"],
    Beta2_Slope     = coeffs["x2"],
    Beta3_Curvature = coeffs["x3"],
    R2              = summary(model)$r.squared
  ))
}

df_factors <- df_interpolated |>
  left_join(lambda_by_country, by = "Country") |>
  group_by(Date, Country, Lambda) |>
  do(estimate_betas(., lambda = unique(.$Lambda))) |>
  ungroup() |>
  select(-Lambda) |>
  arrange(Country, Date)

# ── Proxy empirique corrigé ───────────────────────────────────────────────────
# Pour chaque pays, on calcule numériquement les coefficients du butterfly
# qui isole β3 pour SON lambda, au lieu d'utiliser les maturités fixes de l'article

compute_butterfly_weights <- function(lambda) {
  peak  <- 1.7932 / lambda          # pic du loading β3
  mats  <- c(0.25, peak, 2 * peak)  # court, pic, symétrique long
  
  L2 <- function(tau) (1 - exp(-lambda * tau)) / (lambda * tau)
  L3 <- function(tau) L2(tau) - exp(-lambda * tau)
  
  tau_s <- mats[1]; tau_m <- mats[2]; tau_l <- mats[3]
  
  A_mat <- matrix(c(1, L2(tau_s), 1, L2(tau_l)), nrow = 2)
  b_vec <- c(-2, -2 * L2(tau_m))
  ac    <- solve(A_mat, b_vec)
  
  weights <- c(ac[1], 2, ac[2])
  names(weights) <- paste0("w_", round(mats, 2))
  
  attr(weights, "mats")    <- mats
  attr(weights, "coef_b3") <- sum(weights * sapply(mats, L3))
  attr(weights, "coef_b2") <- sum(weights * sapply(mats, L2))
  attr(weights, "coef_b1") <- sum(weights)
  return(weights)
}

butterfly_weights <- lambda_by_country |>
  rowwise() |>
  mutate(
    weights  = list(compute_butterfly_weights(Lambda)),
    mat_s    = attr(weights, "mats")[1],
    mat_m    = attr(weights, "mats")[2],   # = pic β3
    mat_l    = attr(weights, "mats")[3],   # = 2 x pic β3
    w_short  = weights[1],
    w_medium = weights[2],
    w_long   = weights[3],
    coef_b1  = attr(weights, "coef_b1"),
    coef_b2  = attr(weights, "coef_b2"),
    coef_b3  = attr(weights, "coef_b3")
  ) |>
  ungroup()

message("\n── Butterfly adaptatif par pays ──")
print(butterfly_weights |>
        select(Country, Lambda, mat_s, mat_m, mat_l, coef_b1, coef_b2, coef_b3))

# Proxy empirique avec les poids corrects par pays
df_empirical <- df_interpolated |>
  left_join(butterfly_weights |>
              select(Country, Lambda, mat_s, mat_m, mat_l, w_short, w_medium, w_long),
            by = "Country") |>
  group_by(Date, Country) |>
  summarise(
    Empirical_Level = approx(Maturity, Yield, xout = 30,   rule = 2)$y,
    Empirical_Slope = approx(Maturity, Yield, xout = 30,   rule = 2)$y -
      approx(Maturity, Yield, xout = 0.25, rule = 2)$y,
    # Butterfly avec les maturités propres au pays
    y_s = approx(Maturity, Yield, xout = first(mat_s), rule = 2)$y,
    y_m = approx(Maturity, Yield, xout = first(mat_m), rule = 2)$y,
    y_l = approx(Maturity, Yield, xout = first(mat_l), rule = 2)$y,
    w_s = first(w_short),
    w_m = first(w_medium),
    w_l = first(w_long),
    .groups = "drop"
  ) |>
  mutate(
    Empirical_Curvature = w_s * y_s + w_m * y_m + w_l * y_l
  ) |>
  select(Date, Country, Empirical_Level, Empirical_Slope, Empirical_Curvature)

df_comparison <- left_join(df_factors, df_empirical, by = c("Date", "Country"))

# Plot facteurs — Germany only
df_plot_factors <- df_comparison |>
  left_join(butterfly_weights |> select(Country, coef_b3), by = "Country") |>
  pivot_longer(
    cols      = c(starts_with("Beta"), starts_with("Empirical")),
    names_to  = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Factor_Type = case_when(
      grepl("Level",     Metric) ~ "1. Niveau (Level)",
      grepl("Slope",     Metric) ~ "2. Pente (Slope)",
      grepl("Curvature", Metric) ~ "3. Courbure (Curvature)"
    ),
    Source = ifelse(grepl("Beta", Metric), "Estimé (Modèle)", "Empirique (Proxy)")
  ) |>
  mutate(Value = ifelse(Metric == "Beta2_Slope",     -Value,          Value)) |>
  mutate(Value = ifelse(Metric == "Beta3_Curvature",  Value * coef_b3, Value))

germany_name <- top5_countries[grepl("DE|Germany|Allemagne", top5_countries, ignore.case = TRUE)][1]

p_germany <- df_plot_factors |>
  filter(Country == germany_name) |>
  ggplot(aes(x = Date, y = Value, color = Source, linetype = Source)) +
  geom_line(size = 0.8) +
  facet_wrap(~Factor_Type, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Estimé (Modèle)" = "#BA0C2F", "Empirique (Proxy)" = "black")) +
  scale_linetype_manual(values = c("Estimé (Modèle)" = "solid", "Empirique (Proxy)" = "dashed")) +
  labs(title    = paste0("Dynamique des Facteurs - ", germany_name),
       subtitle = "Modèle Nelson-Siegel vs Proxies Empiriques",
       x = "Temps", y = "Valeur (%)",
       caption  = "Beta2 inversé pour lecture 'Long - Court'") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_germany)

# Corrélations par pays
correlations <- df_comparison |>
  group_by(Country) |>
  summarise(
    Corr_Level     = cor(Beta1_Level,     Empirical_Level,     use = "complete.obs"),
    Corr_Slope     = cor(Beta2_Slope,     Empirical_Slope,     use = "complete.obs"),
    Corr_Curvature = cor(Beta3_Curvature, Empirical_Curvature, use = "complete.obs")
  )
print(correlations)

# --- PHASE 4 : PERFORMANCE IN-SAMPLE ---

ns_yield_calc <- function(lambda, beta1, beta2, beta3, maturity) {
  x2 <- (1 - exp(-lambda * maturity)) / (lambda * maturity)
  x3 <- x2 - exp(-lambda * maturity)
  return(beta1 + beta2 * x2 + beta3 * x3)
}

df_performance <- df_interpolated |>
  left_join(df_factors, by = c("Date", "Country")) |>
  left_join(lambda_by_country, by = "Country") |>   # <-- ajouter le lambda par pays
  mutate(
    Fitted   = ns_yield_calc(Lambda, Beta1_Level, Beta2_Slope, Beta3_Curvature, Maturity),
    Residual = Yield - Fitted
  ) |>
  select(-Lambda)

# Métriques par pays + moyenne des 5
rmse_by_country <- df_performance |>
  group_by(Country) |>
  summarise(
    Mean_Error = round(mean(Residual, na.rm = TRUE), 4),
    RMSE       = round(sqrt(mean(Residual^2, na.rm = TRUE)), 4),
    .groups    = "drop"
  )

rmse_avg <- rmse_by_country |>
  summarise(
    Country    = "MOYENNE (5 pays)",
    Mean_Error = round(mean(Mean_Error), 4),
    RMSE       = round(mean(RMSE), 4)
  )

rmse_final <- bind_rows(rmse_by_country, rmse_avg)
message("\n--- Performance Globale (RMSE en %) ---")
print(rmse_final)

# Snapshot Germany only
last_date <- max(df_performance$Date)

df_points_last <- df_performance |> filter(Date == last_date, Country == germany_name)

mats_fine   <- seq(0.25, 30, by = 0.1)
betas_de    <- df_factors |> filter(Date == last_date, Country == germany_name)
lambda_de   <- lambda_by_country$Lambda[lambda_by_country$Country == germany_name]
yields_fine <- ns_yield_calc(lambda_de, betas_de$Beta1_Level,   # <-- lambda DE
                             betas_de$Beta2_Slope, betas_de$Beta3_Curvature, mats_fine)
df_curve_last <- data.frame(Maturity = mats_fine, Yield = yields_fine)

p_snapshot <- ggplot() +
  geom_line(data  = df_curve_last,
            aes(x = Maturity, y = Yield), color = "#BA0C2F", size = 1.2) +
  geom_point(data = df_points_last,
             aes(x = Maturity, y = Yield), color = "#BA0C2F", size = 3, alpha = 0.6) +
  labs(title    = paste0("Ajustement NS — ", germany_name,
                         " (", format(last_date, "%d/%m/%Y"), ")"),
       subtitle = "Points = Taux Réels | Ligne = Courbe Nelson-Siegel",
       x = "Maturité (Années)", y = "Rendement (%)") +
  theme_minimal()
print(p_snapshot)

# Heatmap résidus — Germany only
p_residuals <- df_performance |>
  filter(Country == germany_name) |>
  ggplot(aes(x = Date, y = as.factor(Maturity), fill = Residual)) +
  geom_tile() +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, name = "Erreur (%)") +
  labs(title    = paste0("Résidus (Heatmap) — ", germany_name),
       subtitle = "Rouge = Modèle trop bas | Bleu = Modèle trop haut",
       x = "Temps", y = "Maturité") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_residuals)

# --- PHASE 5 : MODÉLISATION DYNAMIQUE ---

df_factors    <- df_factors |> arrange(Country, Date)
df_long_factors <- df_factors |>
  pivot_longer(cols = starts_with("Beta"), names_to = "Factor", values_to = "Value")

# Tests stationnarité — tous pays, affichage complet
run_stationarity_tests <- function(series) {
  if (length(series) < 10) return(data.frame(ADF_p_value = NA, KPSS_p_value = NA))
  adf_res  <- try(adf.test(series, alternative = "stationary"), silent = TRUE)
  kpss_res <- try(kpss.test(series, null = "Level"), silent = TRUE)
  data.frame(
    ADF_p_value  = if (inherits(adf_res,  "try-error")) NA else adf_res$p.value,
    KPSS_p_value = if (inherits(kpss_res, "try-error")) NA else kpss_res$p.value
  )
}

stationarity_results <- df_long_factors |>
  group_by(Country, Factor) |>
  reframe(run_stationarity_tests(Value)) |>
  mutate(
    Verdict_ADF  = ifelse(ADF_p_value  < 0.05, "Stationnaire", "Non-Stationnaire"),
    Verdict_KPSS = ifelse(KPSS_p_value > 0.05, "Stationnaire", "Non-Stationnaire"),
    Synthese     = case_when(
      is.na(ADF_p_value) ~ "Erreur",
      Verdict_ADF == "Stationnaire"     & Verdict_KPSS == "Stationnaire"     ~ "SÛR: Stationnaire",
      Verdict_ADF == "Non-Stationnaire" & Verdict_KPSS == "Non-Stationnaire" ~ "SÛR: Non-Stationnaire",
      TRUE ~ "Ambigu"
    )
  )
print(stationarity_results)

# ACF — Germany only
calc_acf_df <- function(series, lags = 20) {
  acf_res <- acf(series, plot = FALSE, lag.max = lags)
  data.frame(Lag = as.numeric(acf_res$lag),
             ACF = as.numeric(acf_res$acf)) |>
    filter(Lag > 0)
}

df_acf_de <- df_long_factors |>
  filter(Country == germany_name) |>
  group_by(Country, Factor) |>
  reframe(calc_acf_df(Value))

p_acf <- ggplot(df_acf_de, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "#BA0C2F", width = 0.7) +
  facet_wrap(~Factor, ncol = 1) +
  geom_hline(yintercept = c(-1.96 / sqrt(nrow(df_factors) / length(top5_countries)),
                            1.96 / sqrt(nrow(df_factors) / length(top5_countries))),
             linetype = "dashed", color = "blue") +
  labs(title    = paste0("ACF des Facteurs — ", germany_name),
       subtitle = "Persistance élevée → modèles AR justifiés",
       x = "Retard (Jours)", y = "Autocorrélation") +
  theme_minimal()
print(p_acf)

# AR(1) — tous pays
fit_ar1 <- function(series) {
  model <- arima(series, order = c(1, 0, 0))
  list(ar_coef   = as.numeric(model$coef["ar1"]),
       intercept = as.numeric(model$coef["intercept"]),
       residuals = as.numeric(residuals(model)))
}

models_results <- df_long_factors |>
  group_by(Country, Factor) |>
  summarise(Model_Out = list(fit_ar1(Value)), .groups = "drop") |>
  mutate(
    Phi       = sapply(Model_Out, function(x) x$ar_coef),
    Intercept = sapply(Model_Out, function(x) x$intercept),
    Residuals = lapply(Model_Out, function(x) x$residuals)
  )

message("\n--- Coefficients AR(1) par pays ---")
print(models_results |> select(Country, Factor, Phi, Intercept))

# ACF résidus — Germany only
df_residuals_acf_de <- models_results |>
  filter(Country == germany_name) |>
  select(Country, Factor, Residuals) |>
  unnest(Residuals) |>
  group_by(Country, Factor) |>
  reframe(calc_acf_df(Residuals))

p_acf_res <- ggplot(df_residuals_acf_de, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "#BA0C2F", width = 0.7) +
  facet_wrap(~Factor, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = c(-1.96 / sqrt(nrow(df_factors) / length(top5_countries)),
                            1.96 / sqrt(nrow(df_factors) / length(top5_countries))),
             linetype = "dashed", color = "blue") +
  labs(title    = paste0("ACF Résidus AR(1) — ", germany_name),
       subtitle = "Absence d'autocorrélation = Modèle validé",
       x = "Retard (Jours)", y = "Autocorrélation Résiduelle") +
  theme_minimal()
print(p_acf_res)