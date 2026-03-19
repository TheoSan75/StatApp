rm(list = ls())

# ==============================================================================
# CORPORATE BOND YIELD CURVE FITTING
# Goal: Show that sovereign curve models (NSS, Splines) are NOT appropriate
# when applied issuer-by-issuer on sparse corporate bond data
# ==============================================================================

# 1. Setup ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(progress)
source("utils.R")

# 2. Load Corporate Data -------------------------------------------------------
corpo <- read_csv("r_data/NewData_Corpo_clean.csv", show_col_types = FALSE) |>
  mutate(across(starts_with("YLD_"), ~ as.numeric(gsub(",", ".", .x))))

yld_cols <- grep("^YLD_", names(corpo), value = TRUE)

corpo_long <- corpo |>
  select(-Maturity) |>
  pivot_longer(cols = all_of(yld_cols),
               names_to  = "date_label",
               values_to = "Yield") |>
  filter(!is.na(Yield)) |>
  mutate(date_label = gsub("YLD_", "", date_label)) |>
  rename(
    Maturity = Residual_Maturity,
    Amt_Out  = `Amt Out`,
    Ticker   = Ticker        # keep as Ticker, not Country
  ) |>
  mutate(Cpn = as.numeric(gsub(",", ".", Cpn)))

message("Loaded ", nrow(corpo_long), " rows | ",
        length(unique(corpo_long$Ticker)), " tickers | ",
        length(unique(corpo_long$date_label)), " dates")

# 3. Helper: safe model fitters ------------------------------------------------

fit_nss_ols_safe <- function(df) {
  tryCatch({
    ssr_cost <- function(params) {
      if (params[5] <= 0 || params[6] <= 0) return(1e9)
      preds <- nss_func(df$Maturity, params)
      sum((df$Yield - preds)^2)
    }
    start_params <- c(b0 = mean(df$Yield), b1 = 0, b2 = 0, b3 = 0, tau1 = 2, tau2 = 10)
    opt <- optim(par = start_params, fn = ssr_cost, control = list(maxit = 5000))
    opt$par
  }, error = function(e) NULL)
}

fit_nss_ecb_safe <- function(df) {
  tryCatch({
    start_params <- c(b0 = mean(df$Yield), b1 = 0, b2 = 0, b3 = 0, tau1 = 1, tau2 = 5)
    opt <- optim(
      par        = start_params,
      fn         = nss_objective,
      maturities = df$Maturity,
      yields     = df$Yield,
      weights    = df$Amt_Out,
      method     = "L-BFGS-B",
      lower      = c(0, -15, -15, -15, 0.1, 0.1),
      upper      = c(15,  15,  15,  15, 10,  30),
      control    = list(maxit = 2000)
    )
    opt$par
  }, error = function(e) NULL)
}

fit_nss_wls_safe <- function(df) {
  tryCatch({
    weights <- df$Amt_Out
    weights[is.na(weights)] <- 0
    if (sum(weights) == 0) weights <- rep(1, nrow(df))
    weights <- weights / sum(weights)
    
    ssr_cost_wls <- function(params) {
      if (params[5] <= 0.1 || params[6] <= 0.5) return(1e9)
      preds     <- nss_func(df$Maturity, params)
      tail_rate <- nss_func(30, params)
      penalty   <- if (tail_rate < -2 || tail_rate > 15) 1000 else 0
      sum(weights * (df$Yield - preds)^2) + penalty
    }
    start_params <- c(b0 = mean(df$Yield), b1 = 0, b2 = 0, b3 = 0, tau1 = 1, tau2 = 5)
    opt <- optim(par = start_params, fn = ssr_cost_wls,
                 method = "Nelder-Mead", control = list(maxit = 10000))
    opt$par
  }, error = function(e) NULL)
}

fit_spline_safe <- function(df) {
  tryCatch({
    if (length(unique(df$Maturity)) < 4) return(NULL)
    weights <- df$Amt_Out
    weights[is.na(weights)] <- 0
    if (sum(weights) == 0) weights <- rep(1, nrow(df))
    
    best_score <- Inf
    best_model <- NULL
    for (s in seq(0.1, 1.5, by = 0.05)) {
      fit <- tryCatch(
        smooth.spline(x = df$Maturity, y = df$Yield, w = weights, spar = s, cv = TRUE),
        error = function(e) NULL
      )
      if (!is.null(fit) && fit$cv.crit < best_score) {
        best_score <- fit$cv.crit
        best_model <- fit
      }
    }
    best_model
  }, error = function(e) NULL)
}

# 4. Main Loop: issuer x date --------------------------------------------------
all_scorecards <- list()
fit_summary    <- list()
tickers <- unique(corpo_long$Ticker)
dates   <- unique(corpo_long$date_label)

pb <- progress_bar$new(
  format = "[:bar] :current/:total tickers | :ticker | ETA: :eta",
  total  = length(tickers),
  clear  = FALSE
)

for (ticker in tickers) {
  pb$tick(tokens = list(ticker = formatC(ticker, width = 8, flag = "-")))
  for (d in dates) {
    
    df <- corpo_long |>
      filter(Ticker == ticker, date_label == d) |>
      filter(!is.na(Maturity), !is.na(Yield)) |>
      rename(Country = Ticker)   # evaluate_model_performance doesn't need this
    # but keeps df self-contained; Amt_Out/Yield/Maturity/Cpn present
    
    n <- nrow(df)
    
    fit_summary[[paste(ticker, d, sep = "_")]] <- data.frame(
      Ticker = ticker, Date = d, N_bonds = n
    )
    
    if (n < 2) next
    
    tag <- paste0(ticker, " | ", d)
    
    # ── NSS OLS ────────────────────────────────────────────────────────────
    params_ols <- fit_nss_ols_safe(df)
    if (!is.null(params_ols)) {
      score <- tryCatch(
        evaluate_model_performance(
          model_name   = paste("NSS_OLS", tag),
          data         = df,
          predict_func = function(m) nss_func(m, params_ols),
          col_weight   = "Amt_Out"
        ),
        error = function(e) NULL
      )
      if (!is.null(score)) {
        all_scorecards[[paste("NSS_OLS", tag)]] <- cbind(
          as.data.frame(score),
          Method = "NSS OLS", Ticker = ticker, Date = d, N_bonds = n
        )
      }
    }
    
    # ── NSS ECB ────────────────────────────────────────────────────────────
    params_ecb <- fit_nss_ecb_safe(df)
    if (!is.null(params_ecb)) {
      score <- tryCatch(
        evaluate_model_performance(
          model_name   = paste("NSS_ECB", tag),
          data         = df,
          predict_func = function(m) nss_func(m, params_ecb),
          col_weight   = "Amt_Out"
        ),
        error = function(e) NULL
      )
      if (!is.null(score)) {
        all_scorecards[[paste("NSS_ECB", tag)]] <- cbind(
          as.data.frame(score),
          Method = "NSS ECB", Ticker = ticker, Date = d, N_bonds = n
        )
      }
    }
    
    # ── NSS WLS ────────────────────────────────────────────────────────────
    params_wls <- fit_nss_wls_safe(df)
    if (!is.null(params_wls)) {
      score <- tryCatch(
        evaluate_model_performance(
          model_name   = paste("NSS_WLS", tag),
          data         = df,
          predict_func = function(m) nss_func(m, params_wls),
          col_weight   = "Amt_Out"
        ),
        error = function(e) NULL
      )
      if (!is.null(score)) {
        all_scorecards[[paste("NSS_WLS", tag)]] <- cbind(
          as.data.frame(score),
          Method = "NSS WLS", Ticker = ticker, Date = d, N_bonds = n
        )
      }
    }
    
    # ── Cubic Spline ───────────────────────────────────────────────────────
    spline_model <- fit_spline_safe(df)
    if (!is.null(spline_model)) {
      score <- tryCatch(
        evaluate_model_performance(
          model_name   = paste("Spline", tag),
          data         = df,
          predict_func = function(m) predict(spline_model, x = m)$y,
          col_weight   = "Amt_Out"
        ),
        error = function(e) NULL
      )
      if (!is.null(score)) {
        all_scorecards[[paste("Spline", tag)]] <- cbind(
          as.data.frame(score),
          Method = "Cubic Spline", Ticker = ticker, Date = d, N_bonds = n
        )
      }
    }
    
  }  # end for d
}  # end for ticker

message("Loop done. Scorecards collected: ", length(all_scorecards))

# 5. Aggregate Results ---------------------------------------------------------

bond_counts <- bind_rows(fit_summary) |>
  group_by(Ticker) |>
  summarise(Avg_N_bonds = mean(N_bonds), Min_N = min(N_bonds), Max_N = max(N_bonds),
            .groups = "drop")

full_scorecard <- bind_rows(all_scorecards)
message("full_scorecard: ", nrow(full_scorecard), " rows")

# Average per issuer x method (across dates)
avg_by_issuer <- full_scorecard |>
  group_by(Ticker, Method) |>
  summarise(
    Avg_W_RMSE     = mean(W_RMSE,         na.rm = TRUE),
    Avg_Price_RMSE = mean(Price_RMSE,      na.rm = TRUE),
    Avg_Wiggliness = mean(Wiggliness,      na.rm = TRUE),
    Avg_Tail_Stab  = mean(Tail_Stab_30_50, na.rm = TRUE),
    N_fits         = n(),
    .groups = "drop"
  ) |>
  arrange(Ticker, Method)

# Average per method only (global headline)
avg_global <- full_scorecard |>
  group_by(Method) |>
  summarise(
    Avg_W_RMSE     = mean(W_RMSE,         na.rm = TRUE),
    Avg_Price_RMSE = mean(Price_RMSE,      na.rm = TRUE),
    Avg_Wiggliness = mean(Wiggliness,      na.rm = TRUE),
    Avg_Tail_Stab  = mean(Tail_Stab_30_50, na.rm = TRUE),
    N_fits         = n(),
    .groups = "drop"
  ) |>
  arrange(Avg_W_RMSE)

# 6. Plot 1: Sparsity distribution ---------------------------------------------

# Top 5 issuers by avg bond count
top5 <- bond_counts |>
  slice_max(Avg_N_bonds, n = 5)

p_sparsity <- ggplot(bond_counts, aes(x = Avg_N_bonds)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 1, fill = "#378ADD", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 6, linetype = "dashed", color = "#E24B4A", linewidth = 0.9) +
  annotate("text", x = 6.3, y = Inf, vjust = 2,
           label = "Minimum NSS (6 params)",
           color = "#E24B4A", size = 3.2, hjust = 0, fontface = "italic") +
  ggrepel::geom_text_repel(
    data = top5,
    aes(x = Avg_N_bonds, y = 0, label = Ticker),
    nudge_y       = 0.05,
    size          = 3.2,
    color         = "#185FA5",
    fontface      = "bold",
    segment.color = "#B5D4F4",
    segment.size  = 0.4,
    direction     = "x",
    min.segment.length = 0
  ) +
  scale_x_continuous(breaks = seq(0, max(bond_counts$Max_N), by = 10)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Distribution du nombre moyen d'obligations par entreprise",
    subtitle = "La plupart des émetteurs ont moins d'obligations que de paramètres NSS (6)",
    x        = "Nombre moyen d'obligations par date d'observation",
    y        = "Part des entreprises"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_sparsity)

# 7. Plot 2: RMSE distribution per method -------------------------------------

p_rmse <- ggplot(full_scorecard, aes(x = Method, y = W_RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_y_log10() +
  labs(
    title    = "Distribution of Weighted RMSE by Model -- Corporate Bonds",
    subtitle = "High variance and extreme values reveal model instability from sparse data",
    x = NULL, y = "W_RMSE (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_rmse)

# 8. Print scorecards ----------------------------------------------------------

message("\n========================================")
message("SCORECARD: Average per Issuer x Method")
message("========================================")
print(avg_by_issuer, n = Inf)

message("\n========================================")
message("SCORECARD: Global Average per Method")
message("========================================")
print(avg_global)

# 9. Save outputs --------------------------------------------------------------
saveRDS(full_scorecard,  "models/corpo_full_scorecard.rds")
write.csv(avg_by_issuer, "r_data/corpo_avg_by_issuer.csv", row.names = FALSE)
write.csv(avg_global,    "r_data/corpo_avg_global.csv",    row.names = FALSE)

message("Done.")