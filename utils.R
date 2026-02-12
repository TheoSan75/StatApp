################################################################################
######    Fonction Calcul de Wiggliness
################################################################################
# Fonction pour calculer la "Wiggliness" (Rugosité) d'une courbe
calc_wiggliness <- function(predict_func, min_mat = 0.25, max_mat = 50, step = 0.25) {
  # On commence à 3 mois (0.25) pour éviter la division par zéro du NSS à t=0
  grid_mat <- seq(min_mat, max_mat, by = step)
  grid_yields <- predict_func(grid_mat)
  # Sécurité supplémentaire : suppression des NA éventuels
  if(any(is.na(grid_yields))) {
    return(NA) 
  }
  second_deriv <- diff(diff(grid_yields))
  # Métrique
  wiggliness_score <- sum(second_deriv^2) * 1000 
  return(wiggliness_score)
}


################################################################################
######    Fonction Pre-Processing ECB compliant
################################################################################
filter_ecb_criteria <- function(df) {
  # 1. Filtre Maturité Résiduelle
  df_filtered <- df %>%
    filter(Maturity >= 0.25 & Maturity <= 30)
  clean_data <- df_filtered
  
  # 2. Filtrage Itératif (2 passages comme recommandé par la BCE)
  for(i in 1:2) {
    # Calcul dynamique du nombre de buckets (cible : ~15 points par bucket)
    n_obs <- nrow(clean_data)
    target_n <- 15 
    n_buckets <- max(1, floor(n_obs / target_n))
    clean_data <- clean_data %>%
      # On crée les groupes par quantiles (densité égale)
      mutate(Bucket = ggplot2::cut_number(Maturity, n = n_buckets)) %>%
      group_by(Bucket) %>%
      mutate(
        mu = mean(Yield, na.rm=TRUE),
        sigma = sd(Yield, na.rm=TRUE),
        Lower = mu - 2*sigma,
        Upper = mu + 2*sigma
      ) %>%
      ungroup() %>%
      # On ne garde que ce qui est dans le tunnel
      filter(Yield >= Lower & Yield <= Upper) %>%
      select(-Bucket, -mu, -sigma, -Lower, -Upper) # Nettoyage des colonnes temporaires
  }
  return(clean_data)
}

################################################################################
######    Fonction Pre-Processing classique
################################################################################
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
    # Nettoyage des NA
    filter(!is.na(Yield), !is.na(Maturity), !is.na(Cpn), Maturity > 0) %>%
    select(Country, MaturityDate, Maturity, Yield, Amt_Out, Cpn) 
  
  return(df_clean)
}

################################################################################
######    Fonction de calcul NSS
################################################################################

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

################################################################################
######    Fonction de calcul de la Macaulay Duration
################################################################################

# (formule fermée)
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

################################################################################
######    Fonction Calcul des métriques de performance
################################################################################
evaluate_model_performance <- function(model_name, data, predict_func, col_weight = "Amt_Out") {
  # 1. Prédictions
  preds <- predict_func(data$Maturity)
  residuals <- data$Yield - preds  # <--- Le vecteur est ici !
  
  # 2. Poids
  weights <- data[[col_weight]]
  if(is.null(weights) | all(is.na(weights))) { weights <- rep(1, nrow(data)) }
  weights[is.na(weights)] <- 0 
  
  # 3. Métriques Standard
  w_rmse <- sqrt(sum(weights * residuals^2, na.rm = TRUE) / sum(weights, na.rm = TRUE))
  wiggliness <- calc_wiggliness(predict_func)
  
  pred_30 <- predict_func(30)
  pred_50 <- predict_func(50)
  tail_diff <- abs(pred_50 - pred_30)
  
  # 4. Métrique Price Error (CORRECTION ICI)
  # On recalcule la duration localement pour être sûr
  durations <- calculate_duration(data$Yield, data$Cpn, data$Maturity)
  
  # On multiplie le vecteur 'residuals' par le vecteur 'durations'
  # (Et non pas data$Residual qui n'existe pas)
  price_rmse <- sqrt(mean((residuals * durations)^2, na.rm = TRUE))
  
  # 5. Bias par Bucket
  data_res <- data %>%
    mutate(
      Residual = residuals,
      Bucket = cut(Maturity, 
                   breaks = c(0, 2, 5, 10, 30, 100), 
                   labels = c("0-2Y", "2-5Y", "5-10Y", "10-30Y", "30Y+"))
    )
  
  bucket_bias <- data_res %>%
    group_by(Bucket) %>%
    summarise(Bias = mean(Residual, na.rm = TRUE), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = Bucket, values_from = Bias, names_prefix = "Bias_")
  
  # 6. Assemblage
  scorecard <- tibble(
    Model = model_name,
    W_RMSE = round(w_rmse, 4),
    Price_RMSE = round(price_rmse, 4),
    Wiggliness = round(wiggliness, 4),
    Tail_Stab_30_50 = round(tail_diff, 4)
  ) %>%
    bind_cols(bucket_bias)
  
  return(scorecard)
}
