# ==============================================================================
# MODÉLISATION HIÉRARCHIQUE DU SPREAD CORPORATE — RIDGE PÉNALISÉE ADAPTATIVE
# ==============================================================================
#
# Méthodologie (5 étapes) :
#
#   Étape 0 : NS3 à τ fixé par groupe (réduction 6→3 paramètres libres)
#   Étape 1 : Calibration de τ_groupe sur émetteurs riches (≥ 8 bonds)
#   Étape 2 : A priori (μ_group, Σ_group) par cellule rating × secteur
#   Étape 3 : Ridge pénalisée avec λ_i = λ_0 / N_i,
#             λ_0 calibré par LOO-CV sur émetteurs modérés (6–9 bonds)
#   Étape 4 : AR(1) sur les séries temporelles des paramètres (dynamique légère)
#   Étape 5 : Métriques complètes in-sample et out-of-sample,
#             comparaison avec Ren (2008)
#
# Métriques calculées (in-sample ET out-of-sample, par classe de sparsité) :
#   — RMSE yield, MAE yield, RMSE prix (Duration × résidu yield)
#   — Biais par bucket de maturité (0–2Y, 2–5Y, 5–10Y, 10–30Y)
#   — Wiggliness (rugosité de la courbe)
#   — Stabilité temporelle des paramètres (σ inter-dates)
#   — Analogues souverains : R², corrélation empirique/modèle des facteurs NS3
#
# Fonctions réutilisées depuis utils.R :
#   nss_func(), calculate_duration(), filter_ecb_criteria(),
#   evaluate_model_performance(), nss_objective(), calc_wiggliness()
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
library(tseries)        # ADF/KPSS pour Section 4
source("utils.R")

set.seed(42)

# ── PALETTES (cohérence graphique avec le reste du projet) ───────────────────
rating_colors <- c(
  "AA"  = "#1B4F8A", "A" = "#2E86C1", "BBB" = "#F39C12",
  "HY"  = "#C0392B", "NR" = "#7F8C8D"
)
sparsity_colors <- c(
  "Très sparse (< 3)"  = "#E74C3C",
  "Sparse (3–5)"       = "#E67E22",
  "Modéré (6–9)"       = "#3498DB",
  "Riche (≥ 10)"       = "#27AE60"
)

# ==============================================================================
# SECTION 0 — CHARGEMENT, PRÉ-TRAITEMENT ET CONSTRUCTION DU SPREAD
# ==============================================================================

# ── 0.1 Données corporate ─────────────────────────────────────────────────────
corpo_raw <- read_csv("r_data/NewData_Corpo_clean.csv", show_col_types = FALSE) |>
  mutate(across(starts_with("YLD_"), ~ suppressWarnings(as.numeric(gsub(",", ".", .x)))))

yld_cols <- grep("^YLD_", names(corpo_raw), value = TRUE)

corpo_long <- corpo_raw |>
  pivot_longer(cols = all_of(yld_cols), names_to = "Date_Label", values_to = "Yield_Corp") |>
  mutate(
    Date     = dmy(gsub("^YLD_", "", Date_Label)),
    Maturity = as.numeric(Residual_Maturity),
    Amt_Out  = suppressWarnings(as.numeric(`Amt Out`)),
    Cpn      = suppressWarnings(as.numeric(Cpn)),
    Yield_Corp = as.numeric(Yield_Corp)
  ) |>
  filter(!is.na(Yield_Corp), !is.na(Maturity), Maturity > 0.25, Maturity <= 30) |>
  rename(
    Rating     = `BBG Composite`,
    Sector_L1  = `BICS Level 1`,
    Sector_L2  = `BICS Level 2`,
    Country_Risk = `Cntry of Risk`
  ) |>
  mutate(
    Rating_Bucket = case_when(
      Rating %in% c("AAA","AA+","AA","AA-")              ~ "AA",
      Rating %in% c("A+","A","A-")                       ~ "A",
      Rating %in% c("BBB+","BBB","BBB-")                 ~ "BBB",
      Rating %in% c("BB+","BB","BB-","B+","B","B-",
                    "CCC","CC","C")                       ~ "HY",
      TRUE                                                ~ "NR"
    ),
    # Classe de sparsité (calculée sur snapshot dernière date plus bas)
    Group_Key = paste(Rating_Bucket, Sector_L1, sep = " | ")
  ) |>
  select(ISIN, `Issuer Name`, Ticker, Rating, Rating_Bucket, Sector_L1, Sector_L2,
         Country_Risk, Date, Maturity, Yield_Corp, Amt_Out, Cpn, Group_Key)

# ── 0.2 Courbe swap risk-free (identique à Ren_2008__On_NewCorpo.R) ──────────
swap_raw <- read_csv("r_data/StatApp_Data1_Swap.csv", show_col_types = FALSE)

clean_swap <- function(df) {
  cvt <- function(t) {
    if (grepl("M", t)) return(as.numeric(gsub("M","",t)) / 12)
    if (grepl("Y", t)) return(as.numeric(gsub("Y","",t)))
    NA
  }
  out <- df |>
    select(TenorStr = 1, SwapRate = 3) |>
    mutate(Maturity = sapply(TenorStr, cvt),
           Yield    = suppressWarnings(as.numeric(as.character(SwapRate)))) |>
    filter(!is.na(Maturity), !is.na(Yield))
  if (mean(out$Yield, na.rm = TRUE) < 0.5) out$Yield <- out$Yield * 100
  out
}

df_swap <- clean_swap(swap_raw)

swap_opt <- optim(
  par    = c(b0=4, b1=-0.5, b2=10, b3=-5, tau1=2.5, tau2=0.7),
  fn     = nss_objective,
  maturities = df_swap$Maturity,
  yields     = df_swap$Yield,
  weights    = rep(1, nrow(df_swap)),
  method = "L-BFGS-B",
  lower  = c(0, -15, -15, -15, 0.1, 0.1),
  upper  = c(20,  15,  15,  15, 15,  30),
  control = list(maxit = 5000, factr = 1e6)
)
p_swap <- swap_opt$par
message("Calibration swap : RMSE = ",
        round(sqrt(mean((df_swap$Yield - nss_func(df_swap$Maturity, p_swap))^2)), 4), "%")

# ── 0.3 Spread observé ────────────────────────────────────────────────────────
corpo_long <- corpo_long |>
  mutate(
    Yield_RF = nss_func(Maturity, p_swap),
    Spread   = Yield_Corp - Yield_RF
  )

# ── 0.4 Classe de sparsité (basée sur la dernière date, utilisée tout au long) ─
last_date  <- max(corpo_long$Date)
dates_all  <- sort(unique(corpo_long$Date))
n_dates    <- length(dates_all)

# Création d'un mapping STRICT : 1 Ticker = 1 Rating = 1 Secteur = 1 Sparsity
# (On garde le nom sparsity_map pour que la suite du script fonctionne)
sparsity_map <- corpo_long |>
  group_by(Ticker) |>
  summarise(
    # On calcule le nombre d'obligations présentes sur la dernière date
    N_bonds_snap = n_distinct(ISIN[Date == last_date]),
    
    # On attribue le Rating et le Group_Key MAJORITAIRES à l'émetteur
    Rating_Bucket = names(sort(table(Rating_Bucket), decreasing = TRUE))[1],
    Group_Key     = names(sort(table(Group_Key), decreasing = TRUE))[1],
    .groups = "drop"
  ) |>
  mutate(
    Sparsity_Class = case_when(
      N_bonds_snap < 3  ~ "Très sparse (< 3)",
      N_bonds_snap < 6  ~ "Sparse (3–5)",
      N_bonds_snap < 10 ~ "Modéré (6–9)",
      TRUE              ~ "Riche (≥ 10)"
    )
  )

# On nettoie corpo_long en écrasant les attributs potentiellement multiples
# par les attributs uniques (majoritaires) du mapping
corpo_long <- corpo_long |>
  select(-Rating_Bucket, -Group_Key) |> # Suppression des colonnes multiples
  left_join(sparsity_map, by = "Ticker")

message("Dataset nettoyé : ", nrow(corpo_long), " obs | ",
        length(unique(corpo_long$ISIN)), " ISINs | ",
        length(unique(corpo_long$Ticker)), " tickers | ",
        n_dates, " dates")

# ==============================================================================
# SECTION 1 — ÉTAPE 0+1 : NS3 ET CALIBRATION DE τ PAR GROUPE
# ==============================================================================
# NS3 (Nelson-Siegel pur) : Spread(τ) = b0 + b1·L(τ;τ1) + b2·C(τ;τ1)
# τ1 est fixé par groupe rating×secteur, estimé sur les émetteurs riches (≥8 bonds)

# ── Fonctions NS3 ─────────────────────────────────────────────────────────────
ns3_loadings <- function(maturity, tau1) {
  # Sécurité pour τ=0 (évite division par zéro)
  tau1 <- pmax(tau1, 0.01)
  x <- maturity / tau1
  L <- (1 - exp(-x)) / x          # Loading Niveau/Pente
  C <- L - exp(-x)                 # Loading Courbure
  cbind(L = L, C = C)
}

ns3_spread <- function(maturity, params, tau1) {
  # params = c(b0, b1, b2)
  ld <- ns3_loadings(maturity, tau1)
  params[1] + params[2] * ld[, "L"] + params[3] * ld[, "C"]
}

# ── Calibration de τ1 pour un sous-dataset donné ─────────────────────────────
calibrate_tau_group <- function(df_group) {
  # Minimise le RMSE sur le spread en ne laissant libre que τ1
  # (b0, b1, b2 sont estimés par OLS conditionnellement à τ1)
  cost_tau <- function(tau1) {
    if (tau1 <= 0.05 || tau1 > 30) return(1e9)
    ld   <- ns3_loadings(df_group$Maturity, tau1)
    X    <- cbind(1, ld[, "L"], ld[, "C"])
    w    <- pmax(df_group$Amt_Out, 0, na.rm = TRUE)
    if (sum(w, na.rm = TRUE) == 0) w <- rep(1, nrow(df_group))
    W    <- diag(w / sum(w))
    # OLS pondéré : β = (X'WX)^{-1} X'Wy
    XtWX <- t(X) %*% W %*% X
    if (abs(det(XtWX)) < 1e-12) return(1e9)
    beta <- solve(XtWX) %*% t(X) %*% W %*% df_group$Spread
    fitted <- X %*% beta
    sqrt(mean(w * (df_group$Spread - fitted)^2) / mean(w))
  }
  opt <- optimize(cost_tau, interval = c(0.1, 15))
  opt$minimum
}

# ── Application : τ par groupe, estimé sur émetteurs riches ──────────────────
rich_tickers <- sparsity_map |> filter(N_bonds_snap >= 8) |> pull(Ticker)

# On utilise toutes les dates pour maximiser l'information sur les émetteurs riches
tau_by_group <- corpo_long |>
  filter(Ticker %in% rich_tickers, !is.na(Spread), abs(Spread) < 10) |>
  group_by(Group_Key) |>
  filter(n() >= 10) |>   # Sécurité : au moins 10 obs dans le groupe
  summarise(
    tau1_group = calibrate_tau_group(cur_data()),
    N_obs_group = n(),
    N_tickers_group = n_distinct(Ticker),
    .groups = "drop"
  )

# Fallback global pour les groupes sans émetteurs riches
tau_global <- corpo_long |>
  filter(Ticker %in% rich_tickers, !is.na(Spread), abs(Spread) < 10) |>
  (\(d) calibrate_tau_group(d))()

message("\n── τ calibrés par groupe ──")
print(tau_by_group)
message("τ global (fallback) : ", round(tau_global, 4))

# Merge τ dans corpo_long
corpo_long <- corpo_long |>
  left_join(tau_by_group |> select(Group_Key, tau1_group), by = "Group_Key") |>
  mutate(tau1_group = ifelse(is.na(tau1_group), tau_global, tau1_group))

# ==============================================================================
# SECTION 2 — ÉTAPE 2 : A PRIORI PAR GROUPE (μ_group, Σ_group)
# ==============================================================================
# Estimé sur les émetteurs riches de chaque groupe.
# Pour les groupes sans émetteurs riches : fallback sur la moyenne globale par rating.

estimate_ns3_ols <- function(df_emetteur, tau1) {
  # OLS pondéré par Amt_Out, retourne c(b0, b1, b2)
  ld <- ns3_loadings(df_emetteur$Maturity, tau1)
  X  <- cbind(1, ld[, "L"], ld[, "C"])
  w  <- pmax(df_emetteur$Amt_Out, 0, na.rm = TRUE)
  if (sum(w, na.rm = TRUE) == 0) w <- rep(1, nrow(df_emetteur))
  w  <- w / sum(w)
  W  <- diag(w)
  XtWX <- t(X) %*% W %*% X
  if (abs(det(XtWX)) < 1e-10) return(rep(NA, 3))
  beta <- as.numeric(solve(XtWX) %*% t(X) %*% W %*% df_emetteur$Spread)
  beta
}

# ── Estimation des paramètres NS3 pour chaque émetteur riche × date ──────────
params_rich <- corpo_long |>
  filter(Ticker %in% rich_tickers, !is.na(Spread), abs(Spread) < 10) |>
  group_by(Ticker, Group_Key, Rating_Bucket, Date) |>
  filter(n() >= 3) |>
  summarise(
    tau1 = first(tau1_group),
    {
      p <- estimate_ns3_ols(cur_data(), first(tau1_group))
      tibble(b0 = p[1], b1 = p[2], b2 = p[3])
    },
    N_obs = n(),
    .groups = "drop"
  ) |>
  filter(!is.na(b0))

# ── Construction de μ_group et Σ_group ───────────────────────────────────────
prior_by_group <- params_rich |>
  group_by(Group_Key) |>
  summarise(
    mu_b0   = mean(b0, na.rm = TRUE),
    mu_b1   = mean(b1, na.rm = TRUE),
    mu_b2   = mean(b2, na.rm = TRUE),
    # Variances et covariances (diagonale de Σ + terme croisé b0-b1)
    var_b0  = var(b0, na.rm = TRUE),
    var_b1  = var(b1, na.rm = TRUE),
    var_b2  = var(b2, na.rm = TRUE),
    cov_b01 = cov(b0, b1, use = "complete.obs"),
    cov_b02 = cov(b0, b2, use = "complete.obs"),
    cov_b12 = cov(b1, b2, use = "complete.obs"),
    N_rich  = n(),
    .groups = "drop"
  ) |>
  # Remplacement des NA (groupes avec 1 seul émetteur riche : pas de variance empirique)
  mutate(across(starts_with("var_"), ~ ifelse(is.na(.x) | .x == 0, 0.1, .x)),
         across(starts_with("cov_"), ~ ifelse(is.na(.x), 0, .x)))

# Fallback : prior par rating bucket (pour groupes sans émetteurs riches)
prior_by_rating <- params_rich |>
  group_by(Rating_Bucket) |>
  summarise(
    mu_b0 = mean(b0, na.rm = TRUE), mu_b1 = mean(b1, na.rm = TRUE),
    mu_b2 = mean(b2, na.rm = TRUE), var_b0 = var(b0, na.rm = TRUE),
    var_b1 = var(b1, na.rm = TRUE), var_b2 = var(b2, na.rm = TRUE),
    cov_b01 = cov(b0, b1, use="complete.obs"),
    cov_b02 = cov(b0, b2, use="complete.obs"),
    cov_b12 = cov(b1, b2, use="complete.obs"),
    .groups = "drop"
  ) |>
  mutate(across(starts_with("var_"), ~ ifelse(is.na(.x) | .x == 0, 0.1, .x)),
         across(starts_with("cov_"), ~ ifelse(is.na(.x), 0, .x)))

message("\n── A priori par groupe (μ_group) ──")
print(prior_by_group |> select(Group_Key, mu_b0, mu_b1, mu_b2, N_rich))

# ==============================================================================
# SECTION 3 — ÉTAPE 3 : RIDGE PÉNALISÉE ADAPTATIVE
# ==============================================================================
# Minimise : Σ_j w_j (s_j - ŝ(τ_j; θ))² + λ_i (θ - μ_g)' Σ_g⁻¹ (θ - μ_g)
# avec λ_i = λ_0 / N_i  (shrinkage adaptatif)
#
# λ_0 est calibré par LOO-CV sur les émetteurs modérés (6–9 bonds).

# ── Fonction de récupération du prior pour un émetteur ───────────────────────
get_prior <- function(group_key, rating_bucket) {
  p <- prior_by_group |> filter(Group_Key == group_key)
  if (nrow(p) == 0) {
    p <- prior_by_rating |> filter(Rating_Bucket == rating_bucket)
  }
  if (nrow(p) == 0) {
    # Fallback ultime : prior neutre
    return(list(mu = c(0.5, -0.1, 0),
                Sigma_inv = diag(c(10, 10, 10))))
  }
  mu <- c(p$mu_b0[1], p$mu_b1[1], p$mu_b2[1])
  Sigma <- matrix(
    c(p$var_b0[1],  p$cov_b01[1], p$cov_b02[1],
      p$cov_b01[1], p$var_b1[1],  p$cov_b12[1],
      p$cov_b02[1], p$cov_b12[1], p$var_b2[1]),
    nrow = 3
  )
  # Régularisation de Σ pour garantir l'inversibilité (ridge sur Σ)
  Sigma <- Sigma + diag(1e-4, 3)
  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) diag(c(10, 10, 10)))
  list(mu = mu, Sigma_inv = Sigma_inv)
}

# ── Estimateur Ridge NS3 ──────────────────────────────────────────────────────
fit_ns3_ridge <- function(df_emetteur, tau1, lambda0, prior) {
  # Données
  mat  <- df_emetteur$Maturity
  spr  <- df_emetteur$Spread
  w    <- pmax(df_emetteur$Amt_Out, 0, na.rm = TRUE)
  if (sum(w, na.rm = TRUE) == 0) w <- rep(1, nrow(df_emetteur))
  w    <- w / sum(w)
  N    <- nrow(df_emetteur)
  
  ld   <- ns3_loadings(mat, tau1)
  X    <- cbind(1, ld[, "L"], ld[, "C"])
  mu   <- prior$mu
  Si   <- prior$Sigma_inv
  
  # λ adaptatif
  lam  <- lambda0 / N
  
  # Solution analytique Ridge pénalisée vers μ_g :
  # θ̂ = (X'WX + λ·Σ⁻¹)⁻¹ (X'Ws + λ·Σ⁻¹·μ)
  W_mat <- diag(w)
  A     <- t(X) %*% W_mat %*% X + lam * Si
  b_vec <- t(X) %*% W_mat %*% spr + lam * Si %*% mu
  
  theta <- tryCatch(
    as.numeric(solve(A) %*% b_vec),
    error = function(e) mu   # Fallback sur le prior si singulier
  )
  theta
}

# ── Calibration de λ_0 par LOO-CV sur émetteurs modérés (6–9 bonds) ──────────
moderate_tickers <- sparsity_map |> filter(N_bonds_snap >= 6, N_bonds_snap <= 9) |> pull(Ticker)

# Grille de λ_0 à tester
lambda_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50)

loo_cv_lambda <- function(lambda0, df_moderate) {
  errors <- c()
  tickers_mod <- unique(df_moderate$Ticker)
  
  for (tk in tickers_mod) {
    df_tk <- df_moderate |> filter(Ticker == tk, Date == last_date, !is.na(Spread))
    if (nrow(df_tk) < 3) next
    
    tau1   <- df_tk$tau1_group[1]
    gk     <- df_tk$Group_Key[1]
    rb     <- df_tk$Rating_Bucket[1]
    prior  <- get_prior(gk, rb)
    
    # LOO : laisser une obligation de côté, estimer sur les autres, prédire sur celle laissée
    for (j in seq_len(nrow(df_tk))) {
      train <- df_tk[-j, ]
      test  <- df_tk[j, ]
      if (nrow(train) < 2) next
      
      theta <- tryCatch(
        fit_ns3_ridge(train, tau1, lambda0, prior),
        error = function(e) prior$mu
      )
      pred  <- ns3_spread(test$Maturity, theta, tau1)
      errors <- c(errors, (test$Spread - pred)^2)
    }
  }
  sqrt(mean(errors, na.rm = TRUE))
}

message("\n── Calibration LOO-CV de λ_0 ──")
df_moderate_loo <- corpo_long |>
  filter(Ticker %in% moderate_tickers, !is.na(Spread), abs(Spread) < 10)

cv_results <- sapply(lambda_grid, function(l) {
  loo_cv_lambda(l, df_moderate_loo)
})
cv_df       <- data.frame(Lambda0 = lambda_grid, LOO_RMSE = cv_results)
best_lambda0 <- lambda_grid[which.min(cv_results)]

message("Résultats LOO-CV :")
print(cv_df)
message("λ_0 optimal : ", best_lambda0)

p_loocv <- ggplot(cv_df, aes(x = log10(Lambda0), y = LOO_RMSE)) +
  geom_line(color = "#2E86C1", linewidth = 1.2) +
  geom_point(size = 3, color = "#2E86C1") +
  geom_point(data = cv_df |> filter(Lambda0 == best_lambda0),
             aes(x = log10(Lambda0), y = LOO_RMSE),
             color = "#E74C3C", size = 5, shape = 18) +
  annotate("text", x = log10(best_lambda0) + 0.3, y = min(cv_results) + 0.002,
           label = paste0("λ_0 optimal = ", best_lambda0),
           color = "#E74C3C", fontface = "bold", hjust = 0) +
  labs(
    title    = "Calibration de λ_0 par LOO-CV (émetteurs modérés, 6–9 bonds)",
    subtitle = "Diamant rouge = λ_0 retenu | Axe x en log₁₀",
    x = "log₁₀(λ_0)", y = "RMSE LOO-CV (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
print(p_loocv)

# ── Estimation de tous les émetteurs × dates ──────────────────────────────────
message("\n── Estimation Ridge NS3 pour tous les émetteurs × dates ──")

all_estimates <- corpo_long |>
  filter(!is.na(Spread), abs(Spread) < 10) |>
  group_by(Ticker, Rating_Bucket, Group_Key, Date, Sparsity_Class) |>
  filter(n() >= 1) |>
  summarise(
    tau1   = first(tau1_group),
    {
      df_sub <- cur_data()
      gk     <- first(Group_Key)
      rb     <- first(Rating_Bucket)
      prior  <- get_prior(gk, rb)
      theta  <- tryCatch(
        fit_ns3_ridge(df_sub, first(tau1_group), best_lambda0, prior),
        error = function(e) get_prior(gk, rb)$mu
      )
      tibble(b0 = theta[1], b1 = theta[2], b2 = theta[3])
    },
    N_obs = n(),
    .groups = "drop"
  ) |>
  filter(!is.na(b0))

message("Estimations produites : ", nrow(all_estimates), " (ticker × date)")

# ==============================================================================
# SECTION 4 — ÉTAPE 4 : COMPOSANTE DYNAMIQUE AR(1)
# ==============================================================================
# On modélise l'évolution temporelle de (b0, b1, b2) par AR(1) par émetteur.
# Reproduit la structure de Diebold-Li.R pour cohérence avec le reste du projet.

df_factors_long <- all_estimates |>
  filter(!is.na(b0)) |>
  pivot_longer(cols = c(b0, b1, b2), names_to = "Factor", values_to = "Value") |>
  arrange(Ticker, Factor, Date)

# ── Tests de stationnarité (ADF + KPSS) ───────────────────────────────────────
run_stationarity_tests <- function(series) {
  if (length(series) < 8) return(data.frame(ADF_p = NA, KPSS_p = NA))
  adf_r  <- try(adf.test(series, alternative = "stationary"),  silent = TRUE)
  kpss_r <- try(kpss.test(series, null = "Level"),             silent = TRUE)
  data.frame(
    ADF_p  = if (inherits(adf_r,  "try-error")) NA else adf_r$p.value,
    KPSS_p = if (inherits(kpss_r, "try-error")) NA else kpss_r$p.value
  )
}

stationarity <- df_factors_long |>
  group_by(Ticker, Factor) |>
  reframe(run_stationarity_tests(Value)) |>
  mutate(
    Verdict_ADF  = ifelse(ADF_p  < 0.05, "Stationnaire", "Non-Stat."),
    Verdict_KPSS = ifelse(KPSS_p > 0.05, "Stationnaire", "Non-Stat.")
  )

message("\n── Tests de stationnarité (résumé) ──")
print(stationarity |> count(Factor, Verdict_ADF, Verdict_KPSS))

# ── Estimation AR(1) par émetteur × facteur ───────────────────────────────────
fit_ar1_safe <- function(series) {
  if (length(series) < 4) return(list(phi = NA, intercept = mean(series, na.rm=TRUE),
                                      resid = rep(NA, length(series))))
  m <- tryCatch(arima(series, order = c(1,0,0)), error = function(e) NULL)
  if (is.null(m)) return(list(phi = 0, intercept = mean(series, na.rm=TRUE),
                              resid = rep(0, length(series))))
  list(phi = as.numeric(m$coef["ar1"]),
       intercept = as.numeric(m$coef["intercept"]),
       resid = as.numeric(residuals(m)))
}

ar1_models <- df_factors_long |>
  group_by(Ticker, Factor) |>
  summarise(
    ar1_out = list(fit_ar1_safe(Value)),
    .groups = "drop"
  ) |>
  mutate(
    Phi       = sapply(ar1_out, \(x) x$phi),
    Intercept = sapply(ar1_out, \(x) x$intercept)
  )

message("\n── Coefficients AR(1) moyens par facteur ──")
print(ar1_models |> group_by(Factor) |>
        summarise(Phi_Mean = round(mean(Phi, na.rm=TRUE), 3),
                  Phi_SD   = round(sd(Phi,   na.rm=TRUE), 3)))

# ── ACF des facteurs (subset émetteurs riches pour lisibilité) ────────────────
calc_acf_df <- function(series, lags = 12) {
  acf_res <- acf(series, plot = FALSE, lag.max = lags)
  data.frame(Lag = as.numeric(acf_res$lag), ACF = as.numeric(acf_res$acf)) |>
    filter(Lag > 0)
}

acf_rich <- df_factors_long |>
  filter(Ticker %in% rich_tickers[1:min(6, length(rich_tickers))]) |>
  group_by(Ticker, Factor) |>
  reframe(calc_acf_df(Value)) |>
  mutate(Factor_Label = recode(Factor, b0="β0 (Niveau)", b1="β1 (Pente)", b2="β2 (Courbure)"))

p_acf <- ggplot(acf_rich, aes(x = Lag, y = ACF, fill = Ticker)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(n_dates),
             linetype = "dashed", color = "blue") +
  facet_wrap(~Factor_Label, ncol = 1, scales = "free_y") +
  labs(title = "ACF des Facteurs NS3 (émetteurs riches)",
       subtitle = "Persistance → AR(1) justifié | Bandes bleues = seuils 5%",
       x = "Retard", y = "Autocorrélation") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
print(p_acf)

# ==============================================================================
# SECTION 5 — ÉTAPE 5 : MÉTRIQUES COMPLÈTES
# ==============================================================================
# Pour chaque émetteur × date, on calcule les résidus in-sample.
# Out-of-sample : on utilise les prévisions AR(1) pour les 2 dernières dates.

# ── 5.1 Reconstruction des résidus in-sample ──────────────────────────────────
corpo_with_fit <- corpo_long |>
  filter(!is.na(Spread), abs(Spread) < 10) |>
  left_join(all_estimates |> select(Ticker, Date, b0, b1, b2, tau1, N_obs),
            by = c("Ticker", "Date")) |>
  filter(!is.na(b0)) |>
  mutate(
    Spread_Fitted  = ns3_spread(Maturity, cbind(b0, b1, b2)[1,], tau1),
    # Vectorisé proprement
    Spread_Fitted  = mapply(function(m, b0, b1, b2, t1) {
      ns3_spread(m, c(b0, b1, b2), t1)
    }, Maturity, b0, b1, b2, tau1),
    Yield_Fitted   = Yield_RF + Spread_Fitted,
    Resid_Yield    = Yield_Corp - Yield_Fitted,
    Resid_Spread   = Spread - Spread_Fitted,
    Duration_Mac   = calculate_duration(Yield_Corp, Cpn, Maturity),
    Resid_Price    = -Resid_Yield * Duration_Mac   # ΔP ≈ -D × Δy
  )

# ── 5.2 Prévisions out-of-sample (AR(1) sur les 2 dernières dates) ───────────
n_oos     <- 2
cut_date  <- dates_all[n_dates - n_oos]
oos_dates <- dates_all[(n_dates - n_oos + 1):n_dates]

# Paramètres prévus pour les dates OOS (prévision 1-step ahead AR(1))
preds_oos <- all_estimates |>
  filter(Date <= cut_date) |>
  left_join(ar1_models |> select(Ticker, Factor, Phi, Intercept), by = "Ticker") |>
  # On prend la dernière valeur connue et on applique AR(1)
  group_by(Ticker, Rating_Bucket, Group_Key, Sparsity_Class, tau1) |>
  summarise(
    b0_last = last(b0[order(Date)]),
    b1_last = last(b1[order(Date)]),
    b2_last = last(b2[order(Date)]),
    .groups = "drop"
  ) |>
  left_join(ar1_models |> filter(Factor == "b0") |> select(Ticker, Phi_b0=Phi, Int_b0=Intercept), by="Ticker") |>
  left_join(ar1_models |> filter(Factor == "b1") |> select(Ticker, Phi_b1=Phi, Int_b1=Intercept), by="Ticker") |>
  left_join(ar1_models |> filter(Factor == "b2") |> select(Ticker, Phi_b2=Phi, Int_b2=Intercept), by="Ticker") |>
  mutate(
    across(starts_with("Phi_"), ~ ifelse(is.na(.x), 0.8, .x)),
    across(starts_with("Int_"),  ~ ifelse(is.na(.x), 0,   .x)),
    b0_pred = Int_b0 + Phi_b0 * (b0_last - Int_b0),
    b1_pred = Int_b1 + Phi_b1 * (b1_last - Int_b1),
    b2_pred = Int_b2 + Phi_b2 * (b2_last - Int_b2)
  )

# On applique les prévisions aux dates OOS
corpo_oos <- corpo_long |>
  filter(Date %in% oos_dates, !is.na(Spread), abs(Spread) < 10) |>
  left_join(preds_oos |> select(Ticker, tau1, b0_pred, b1_pred, b2_pred),
            by = "Ticker") |>
  filter(!is.na(b0_pred)) |>
  mutate(
    Spread_Pred   = mapply(function(m, b0, b1, b2, t1) ns3_spread(m, c(b0,b1,b2), t1),
                           Maturity, b0_pred, b1_pred, b2_pred, tau1),
    Yield_Pred    = Yield_RF + Spread_Pred,
    Resid_Yield   = Yield_Corp - Yield_Pred,
    Resid_Spread  = Spread - Spread_Pred,
    Duration_Mac  = calculate_duration(Yield_Corp, Cpn, Maturity),
    Resid_Price   = -Resid_Yield * Duration_Mac
  )

# ── 5.3 Fonction de calcul des métriques complètes ───────────────────────────
compute_metrics <- function(df, sample_label) {
  df |>
    mutate(
      Mat_Bucket = cut(Maturity, breaks = c(0, 2, 5, 10, 30),
                       labels = c("0–2Y", "2–5Y", "5–10Y", "10–30Y"))
    ) |>
    group_by(Sparsity_Class, Mat_Bucket) |>
    summarise(
      Sample      = sample_label,
      N_obs       = n(),
      # ── Métriques yield ──
      RMSE_Yield  = round(sqrt(mean(Resid_Yield^2,  na.rm=TRUE)), 5),
      MAE_Yield   = round(mean(abs(Resid_Yield),    na.rm=TRUE), 5),
      Bias_Yield  = round(mean(Resid_Yield,          na.rm=TRUE), 5),
      # ── Métriques prix (ΔP = -D × Δy) ──
      RMSE_Price  = round(sqrt(mean(Resid_Price^2,  na.rm=TRUE)), 5),
      MAE_Price   = round(mean(abs(Resid_Price),    na.rm=TRUE), 5),
      # ── Métriques spread ──
      RMSE_Spread = round(sqrt(mean(Resid_Spread^2, na.rm=TRUE)), 5),
      MAE_Spread  = round(mean(abs(Resid_Spread),   na.rm=TRUE), 5),
      .groups = "drop"
    )
}

metrics_is  <- compute_metrics(corpo_with_fit, "In-Sample")
metrics_oos <- compute_metrics(corpo_oos,       "Out-of-Sample")
metrics_all <- bind_rows(metrics_is, metrics_oos)

message("\n── Métriques in-sample par classe de sparsité × bucket de maturité ──")
print(as.data.frame(metrics_is))
message("\n── Métriques out-of-sample ──")
print(as.data.frame(metrics_oos))

# ── 5.4 Métriques agrégées (résumé global — analogues souverains) ─────────────
# Analogue du tableau RMSE global de Diebold-Li.R (Phase 4)
metrics_global <- bind_rows(
  corpo_with_fit |>
    group_by(Sparsity_Class) |>
    summarise(
      Sample = "In-Sample", N_obs = n(),
      RMSE_Yield  = round(sqrt(mean(Resid_Yield^2,  na.rm=TRUE)), 5),
      MAE_Yield   = round(mean(abs(Resid_Yield),    na.rm=TRUE), 5),
      Bias_Yield  = round(mean(Resid_Yield,          na.rm=TRUE), 5),
      RMSE_Price  = round(sqrt(mean(Resid_Price^2,  na.rm=TRUE)), 5),
      RMSE_Spread = round(sqrt(mean(Resid_Spread^2, na.rm=TRUE)), 5),
      .groups = "drop"
    ),
  corpo_oos |>
    group_by(Sparsity_Class) |>
    summarise(
      Sample = "Out-of-Sample", N_obs = n(),
      RMSE_Yield  = round(sqrt(mean(Resid_Yield^2,  na.rm=TRUE)), 5),
      MAE_Yield   = round(mean(abs(Resid_Yield),    na.rm=TRUE), 5),
      Bias_Yield  = round(mean(Resid_Yield,          na.rm=TRUE), 5),
      RMSE_Price  = round(sqrt(mean(Resid_Price^2,  na.rm=TRUE)), 5),
      RMSE_Spread = round(sqrt(mean(Resid_Spread^2, na.rm=TRUE)), 5),
      .groups = "drop"
    )
) |> arrange(Sample, Sparsity_Class)

# ── Stabilité temporelle des paramètres (analogue corrélations Diebold-Li) ───
stability_params <- all_estimates |>
  group_by(Ticker, Sparsity_Class) |>
  summarise(
    SD_b0 = sd(b0, na.rm=TRUE),
    SD_b1 = sd(b1, na.rm=TRUE),
    SD_b2 = sd(b2, na.rm=TRUE),
    .groups = "drop"
  ) |>
  group_by(Sparsity_Class) |>
  summarise(across(starts_with("SD_"), ~ round(mean(.x, na.rm=TRUE), 4)))

# ── R² par émetteur × date ────────────────────────────────────────────────────
r2_by_emetteur <- corpo_with_fit |>
  group_by(Ticker, Date, Sparsity_Class, Rating_Bucket) |>
  summarise(
    SS_res = sum(Resid_Yield^2, na.rm=TRUE),
    SS_tot = sum((Yield_Corp - mean(Yield_Corp, na.rm=TRUE))^2, na.rm=TRUE),
    R2     = ifelse(SS_tot > 0, 1 - SS_res/SS_tot, NA),
    N_obs  = n(),
    .groups = "drop"
  )

message("\n── Résumé global des métriques ──")
print(as.data.frame(metrics_global))
message("\n── Stabilité temporelle des paramètres ──")
print(stability_params)
message("\n── R² médian par classe de sparsité ──")
print(r2_by_emetteur |> group_by(Sparsity_Class) |>
        summarise(R2_Median = round(median(R2, na.rm=TRUE), 4),
                  R2_Mean   = round(mean(R2,   na.rm=TRUE), 4)))

# ==============================================================================
# SECTION 6 — VISUALISATIONS POUR LE RAPPORT
# ==============================================================================

# ── VIZ 1 : Courbes de spread pour des émetteurs représentatifs ───────────────
# Un par classe de sparsité → illustre le comportement du shrinkage

grid_mat <- seq(0.5, 30, by = 0.25)

# Sélection d'un représentant par classe
reps <- c(
  "Riche (≥ 10)"      = rich_tickers[1],
  "Modéré (6–9)"      = (sparsity_map |> filter(Sparsity_Class=="Modéré (6–9)") |> pull(Ticker))[1],
  "Sparse (3–5)"      = (sparsity_map |> filter(Sparsity_Class=="Sparse (3–5)") |> pull(Ticker))[1],
  "Très sparse (< 3)" = (sparsity_map |> filter(Sparsity_Class=="Très sparse (< 3)") |> pull(Ticker))[1]
)
reps <- reps[!is.na(reps)]

df_curves_viz <- bind_rows(lapply(names(reps), function(cls) {
  tk     <- reps[[cls]]
  est_tk <- all_estimates |> filter(Ticker == tk, Date == last_date)
  if (nrow(est_tk) == 0) return(NULL)
  
  prior  <- get_prior(est_tk$Group_Key[1], est_tk$Rating_Bucket[1])
  
  # Courbe estimée (Ridge)
  spread_curve <- ns3_spread(grid_mat, c(est_tk$b0[1], est_tk$b1[1], est_tk$b2[1]),
                             est_tk$tau1[1])
  # Courbe a priori (groupe)
  spread_prior <- ns3_spread(grid_mat, prior$mu, est_tk$tau1[1])
  
  bind_rows(
    data.frame(Maturity=grid_mat, Spread=spread_curve, Type="Estimée (Ridge)",
               Ticker=tk, Class=cls),
    data.frame(Maturity=grid_mat, Spread=spread_prior, Type="A priori (groupe)",
               Ticker=tk, Class=cls)
  )
}))

# Points observés
df_obs_viz <- corpo_long |>
  filter(Ticker %in% reps, Date == last_date, !is.na(Spread)) |>
  left_join(data.frame(Ticker=reps, Class=names(reps)), by="Ticker")

p_viz1 <- ggplot() +
  geom_line(data = df_curves_viz,
            aes(x=Maturity, y=Spread, color=Type, linetype=Type), linewidth=1) +
  geom_point(data = df_obs_viz,
             aes(x=Maturity, y=Spread), color="black", size=2.5, alpha=0.7) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  facet_wrap(~Class, scales="free_y", ncol=2) +
  scale_color_manual(values=c("Estimée (Ridge)"="#E74C3C","A priori (groupe)"="#95A5A6"),
                     name=NULL) +
  scale_linetype_manual(values=c("Estimée (Ridge)"="solid","A priori (groupe)"="dashed"),
                        name=NULL) +
  scale_x_continuous(breaks=seq(0,30,5)) +
  labs(
    title    = "Viz 1 — Courbes de spread NS3 : Estimée vs A priori (par classe de sparsité)",
    subtitle = paste0("Points noirs = obligations observées | Date : ",
                      format(last_date, "%d/%m/%Y")),
    x = "Maturité (Années)", y = "Spread vs Swap (%)"
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
print(p_viz1)

# ── VIZ 2 : Heatmap des résidus in-sample (analogue Diebold-Li heatmap) ───────
# Pour les 10 émetteurs avec le plus d'observations

top_tickers_heatmap <- corpo_with_fit |>
  count(Ticker) |> slice_max(n, n=10) |> pull(Ticker)

p_viz2 <- corpo_with_fit |>
  filter(Ticker %in% top_tickers_heatmap) |>
  mutate(Mat_Bucket = cut(Maturity, breaks=c(0,2,5,10,30),
                          labels=c("0–2Y","2–5Y","5–10Y","10–30Y"))) |>
  group_by(Ticker, Date, Mat_Bucket) |>
  summarise(Resid_Mean = mean(Resid_Yield, na.rm=TRUE), .groups="drop") |>
  ggplot(aes(x=Date, y=Ticker, fill=Resid_Mean)) +
  geom_tile(color="white", linewidth=0.3) +
  facet_wrap(~Mat_Bucket, ncol=2) +
  scale_fill_gradient2(low="#B2182B", mid="white", high="#2166AC",
                       midpoint=0, name="Résidu yield (%)") +
  scale_x_date(date_labels="%d/%m", date_breaks="3 weeks") +
  labs(
    title    = "Viz 2 — Heatmap des résidus in-sample (top 10 émetteurs)",
    subtitle = "Rouge = modèle trop bas | Bleu = modèle trop haut",
    x=NULL, y=NULL
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=8))
print(p_viz2)

# ── VIZ 3 : RMSE in vs out-of-sample par classe de sparsité (barplot comparatif) ─
p_viz3 <- metrics_global |>
  filter(!is.na(Sparsity_Class)) |>
  mutate(Sparsity_Class = factor(Sparsity_Class,
                                 levels=c("Très sparse (< 3)","Sparse (3–5)",
                                          "Modéré (6–9)","Riche (≥ 10)"))) |>
  ggplot(aes(x=Sparsity_Class, y=RMSE_Yield, fill=Sample)) +
  geom_col(position="dodge", alpha=0.85) +
  geom_text(aes(label=round(RMSE_Yield*100,2)), position=position_dodge(width=0.9),
            vjust=-0.4, size=3, fontface="bold") +
  scale_fill_manual(values=c("In-Sample"="#2E86C1","Out-of-Sample"="#E74C3C"), name=NULL) +
  labs(
    title    = "Viz 3 — RMSE yield (%) : In-Sample vs Out-of-Sample par classe de sparsité",
    subtitle = "Le shrinkage réduit l'écart IS/OOS pour les émetteurs sparse",
    x=NULL, y="RMSE Yield (%)"
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=20,hjust=1),
        legend.position="bottom")
print(p_viz3)

# ── VIZ 4 : Dynamique temporelle des facteurs NS3 (émetteurs riches) ──────────
# Analogue direct des plots de Diebold-Li Phase 3

df_factors_viz <- all_estimates |>
  filter(Ticker %in% rich_tickers[1:4]) |>
  select(Ticker, Date, b0, b1, b2, Sparsity_Class) |>
  pivot_longer(cols=c(b0,b1,b2), names_to="Factor", values_to="Value") |>
  mutate(Factor_Label = recode(Factor,
                               b0 = "β0 — Niveau du spread",
                               b1 = "β1 — Pente",
                               b2 = "β2 — Courbure"))

p_viz4 <- ggplot(df_factors_viz, aes(x=Date, y=Value, color=Ticker)) +
  geom_line(linewidth=0.8) +
  geom_point(size=2) +
  facet_wrap(~Factor_Label, ncol=1, scales="free_y") +
  scale_x_date(date_labels="%d/%m", date_breaks="2 weeks") +
  labs(
    title    = "Viz 4 — Dynamique temporelle des facteurs NS3 (émetteurs riches)",
    subtitle = "Analogues des facteurs Niveau/Pente/Courbure de Diebold-Li",
    x=NULL, y="Valeur du paramètre"
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.position="bottom")
print(p_viz4)

# ── VIZ 5 : Comparaison Ridge vs Ren (2008) — RMSE yield par rating bucket ───
# On charge les résultats Ren si disponibles, sinon on le note

ren_results_path <- "models/ren_2008_perf.rds"

if (file.exists(ren_results_path)) {
  ren_perf <- readRDS(ren_results_path)
  
  compare_df <- bind_rows(
    metrics_global |>
      filter(Sample == "In-Sample") |>
      left_join(sparsity_map |> select(Ticker, Sparsity_Class) |> distinct(),
                by="Sparsity_Class") |>
      mutate(Method = "Ridge NS3 (Hiérarchique)") |>
      select(Method, Sparsity_Class, RMSE_Yield),
    ren_perf |> mutate(Method = "Ren (2008)", Sparsity_Class = Rating_Bucket)
  )
  
  p_viz5 <- ggplot(compare_df, aes(x=Sparsity_Class, y=RMSE_Yield, fill=Method)) +
    geom_col(position="dodge", alpha=0.85) +
    scale_fill_manual(values=c("Ridge NS3 (Hiérarchique)"="#2E86C1","Ren (2008)"="#E74C3C"),
                      name=NULL) +
    labs(title="Viz 5 — Ridge NS3 vs Ren (2008) : RMSE yield in-sample",
         x=NULL, y="RMSE Yield (%)") +
    theme_minimal() +
    theme(plot.title=element_text(face="bold"), legend.position="bottom")
  print(p_viz5)
} else {
  message("Viz 5 : fichier Ren non trouvé — générer d'abord Ren_2008__On_NewCorpo.R")
}

# ── VIZ 6 : Distribution des R² par classe de sparsité ───────────────────────
p_viz6 <- r2_by_emetteur |>
  filter(!is.na(R2), R2 >= 0) |>
  mutate(Sparsity_Class = factor(Sparsity_Class,
                                 levels=c("Très sparse (< 3)","Sparse (3–5)",
                                          "Modéré (6–9)","Riche (≥ 10)"))) |>
  ggplot(aes(x=Sparsity_Class, y=R2, fill=Sparsity_Class)) +
  geom_boxplot(alpha=0.8, outlier.alpha=0.3, outlier.size=1) +
  geom_hline(yintercept=0.9, linetype="dashed", color="gray40") +
  scale_fill_manual(values=sparsity_colors, guide="none") +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title    = "Viz 6 — Distribution du R² in-sample par classe de sparsité",
    subtitle = "Ligne tiretée = seuil R² = 0.90",
    x=NULL, y="R²"
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=15,hjust=1))
print(p_viz6)

# ── VIZ 7 : Distribution du shrinkage appliqué (||θ̂ - μ|| / ||θ̂_OLS - μ||) ──
# Montre à quel point chaque émetteur est "tiré" vers l'a priori

shrinkage_measure <- all_estimates |>
  filter(Date == last_date) |>
  left_join(
    corpo_long |>
      filter(Date == last_date, !is.na(Spread)) |>
      group_by(Ticker, Group_Key, Rating_Bucket) |>
      summarise(.groups="drop"),
    by = c("Ticker","Group_Key","Rating_Bucket")
  ) |>
  rowwise() |>
  mutate(
    prior    = list(get_prior(Group_Key, Rating_Bucket)),
    mu_b0    = prior$mu[1],
    mu_b1    = prior$mu[2],
    mu_b2    = prior$mu[3],
    Dist_To_Prior = sqrt((b0-mu_b0)^2 + (b1-mu_b1)^2 + (b2-mu_b2)^2)
  ) |>
  ungroup()
# Le left_join inutile a été supprimé ici, car Sparsity_Class est déjà dans all_estimates

p_viz7 <- shrinkage_measure |>
  filter(!is.na(Sparsity_Class)) |>
  mutate(Sparsity_Class = factor(Sparsity_Class,
                                 levels=c("Très sparse (< 3)","Sparse (3–5)",
                                          "Modéré (6–9)","Riche (≥ 10)"))) |>
  ggplot(aes(x=Sparsity_Class, y=Dist_To_Prior, fill=Sparsity_Class)) +
  geom_boxplot(alpha=0.8, outlier.alpha=0.4) +
  scale_fill_manual(values=sparsity_colors, guide="none") +
  labs(
    title    = "Viz 7 — Distance ||θ̂ - μ_groupe|| par classe de sparsité",
    subtitle = "Distance faible = forte attraction vers l'a priori (shrinkage effectif)",
    x=NULL, y="Distance euclidienne au prior"
  ) +
  theme_minimal() +
  theme(plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=15,hjust=1))
print(p_viz7)

# ==============================================================================
# SECTION 7 — TABLEAU DE SYNTHÈSE FINAL
# ==============================================================================

message("\n", paste(rep("=",65), collapse=""))
message("TABLEAU DE SYNTHÈSE — MODÈLE HIÉRARCHIQUE RIDGE NS3")
message(paste(rep("=",65), collapse=""))

message("\n[1] MÉTRIQUES GLOBALES (in-sample)")
print(metrics_global |> filter(Sample=="In-Sample") |>
        select(Sparsity_Class, N_obs, RMSE_Yield, MAE_Yield, Bias_Yield,
               RMSE_Price, RMSE_Spread))

message("\n[2] MÉTRIQUES GLOBALES (out-of-sample)")
print(metrics_global |> filter(Sample=="Out-of-Sample") |>
        select(Sparsity_Class, N_obs, RMSE_Yield, MAE_Yield, Bias_Yield,
               RMSE_Price, RMSE_Spread))

message("\n[3] STABILITÉ TEMPORELLE DES PARAMÈTRES (σ inter-dates)")
print(stability_params)

message("\n[4] R² MÉDIAN PAR CLASSE DE SPARSITÉ")
print(r2_by_emetteur |> group_by(Sparsity_Class) |>
        summarise(R2_Median=round(median(R2,na.rm=TRUE),4),
                  R2_P10=round(quantile(R2,0.1,na.rm=TRUE),4),
                  R2_P90=round(quantile(R2,0.9,na.rm=TRUE),4)))

message("\n[5] λ_0 RETENU : ", best_lambda0)
message("[6] τ GLOBAL (fallback) : ", round(tau_global, 4))
message(paste(rep("=",65), collapse=""))