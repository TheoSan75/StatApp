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