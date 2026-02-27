rm(list = ls())
library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(viridis)
library(knitr)
source("utils.R")

# --- 1. IMPORTATION ---

# On lit le fichier France puis Allemagne
df_france <- read_excel("r_data/souverain_cross_temporel_france.xlsx") %>% 
  mutate(Country = "France")
df_germany <- read_excel("r_data/souverain_cross_temporel_allemagne.xlsx") %>% 
  mutate(Country = "Germany")


# --- 2. FUSION (MERGE) ---
# On colle les deux l'un en dessous de l'autre
df_global <- bind_rows(df_france, df_germany)

# --- 3. FORMATAGE FINAL (Sécurité) ---
df_final <- df_global %>%
  mutate(
    # On s'assure que R comprend bien que ce sont des dates
    Date = as.Date(Date_observation),
    Maturity_Date = as.Date(Maturity),
    Yield = as.numeric(Yield),
    Time_To_Maturity = as.numeric(Time_To_Maturity)
  ) %>%
  select(Date, Country, Time_To_Maturity, Yield) %>%
  arrange(Date, Country, Time_To_Maturity)


# --- 3. EDA ---

# --- A. CRÉATION DES BUCKETS (CATÉGORIES) ---
df_eda <- df_final %>%
  mutate(
    Maturity_Bucket = cut(Time_To_Maturity, 
                          breaks = c(0, 2, 5, 10, 30, 100),
                          labels = c("Court-Terme (<2Y)", 
                                     "Moyen-Terme (2-5Y)", 
                                     "Intermédiaire (5-10Y)", 
                                     "Long-Terme (10-30Y)", 
                                     "Ultra-Long (>30Y)"),
                          include.lowest = TRUE)
  )

# --- B. LE TABLEAU SYNTHÉTIQUE ---
summary_table <- df_eda %>%
  group_by(Country, Maturity_Bucket) %>%
  summarise(
    N_Obs = n(),
    Yield_Mean = round(mean(Yield, na.rm = TRUE), 3),
    Yield_Min = round(min(Yield, na.rm = TRUE), 3),
    Yield_Max = round(max(Yield, na.rm = TRUE), 3),
    Volatilite_SD = round(sd(Yield, na.rm = TRUE), 3),
    .groups = 'drop'
  )

# Affichage propre dans la console
print(as.data.frame(summary_table))


# --- 1. Graphique pour la FRANCE ---
p_france <- df_eda %>%
  filter(Country == "France") %>%
  ggplot(aes(x = Time_To_Maturity, y = Yield, group = Date, color = as.numeric(Date))) + # On force en numeric pour le gradient
  geom_line(alpha = 0.25) +
  
  # C'est ICI que la magie opère pour la légende :
  scale_color_viridis_c(
    option = "magma", 
    name = "Année",            # Titre de la légende plus propre
    labels = date_labeller     # On applique notre formateur
  ) +
  
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  theme_minimal() +
  labs(
    title = "Dynamique des Taux - FRANCE",
    x = "Maturité (Années)",
    y = "Yield (%)"
  )

# --- 2. Graphique pour l'ALLEMAGNE ---
p_germany <- df_eda %>%
  filter(Country == "Germany") %>%
  ggplot(aes(x = Time_To_Maturity, y = Yield, group = Date, color = as.numeric(Date))) +
  geom_line(alpha = 0.25) +
  
  # Même correction ici
  scale_color_viridis_c(
    option = "magma", 
    name = "Année", 
    labels = date_labeller
  ) +
  
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  theme_minimal() +
  labs(
    title = "Dynamique des Taux - ALLEMAGNE",
    x = "Maturité (Années)",
    y = "Yield (%)"
  )

# --- AFFICHAGE ---
print(p_france)
print(p_germany)

p3 <- ggplot(df_eda, aes(x = Maturity_Bucket, y = Yield, fill = Country)) +
  geom_boxplot(outlier.alpha = 0.2, alpha = 0.7) +
  scale_fill_manual(values = c("France" = "#0055A4", "Germany" = "#BA0C2F")) + # Couleurs officielles
  theme_minimal() +
  labs(
    title = "Boxplot des Taux par Maturité",
    x = "Segment",
    y = "Yield (%)"
  ) +
  theme(legend.position = "top")

print(p3)




# --- PHASE 1 : NETTOYAGE ET CONSTRUCTION DE LA GRILLE (Fidèle à Diebold-Li) ---

# 1. Nettoyage des données (Data Cleaning) selon Diebold-Li 
# On enlève les maturités très courtes (< 1 mois) qui sont souvent illiquides ou bruitées.
# Diebold-Li enlèvent les bills < 1 mois.
df_cleaned <- df_final %>%
  filter(Time_To_Maturity >= 1/12) # On garde uniquement ce qui a plus d'un mois

# 2. Définition des maturités cibles EXACTES de Diebold-Li 
# Converties en années. J'ajoute 15, 20 et 30 pour la France/Allemagne.
target_maturities <- c(
  0.25, 0.5, 0.75, 1,    # 3, 6, 9, 12 mois
  1.25, 1.5, 1.75, 2,    # 15, 18, 21, 24 mois
  2.5, 3, 4, 5,          # 30, 36, 48, 60 mois
  6, 7, 8, 9, 10,        # 72, 84, 96, 108, 120 mois
  15, 20, 30             # Ajout spécifique pour nos données (Long terme)
)

# 3. Fonction d'interpolation linéaire
# Méthode citée explicitement dans le papier: 
# "we linearly interpolate nearby maturities to pool into fixed maturities"
interpolate_yields <- function(sub_df, targets) {
  
  # On trie par maturité pour être sûr
  sub_df <- sub_df %>% arrange(Time_To_Maturity)
  
  # Sécurité : il faut au moins 2 points pour interpoler
  if(nrow(sub_df) < 2) return(data.frame(Maturity = targets, Yield = NA))
  
  # Interpolation
  approx_res <- approx(
    x = sub_df$Time_To_Maturity, 
    y = sub_df$Yield, 
    xout = targets,
    method = "linear",
    rule = 2 # Extrapolation plate si on dépasse les bornes (ex: on n'a pas de 30 ans ce jour-là)
  )
  
  return(data.frame(Maturity = approx_res$x, Yield = approx_res$y))
}

# 4. Application à l'ensemble du dataset
df_interpolated <- df_cleaned %>%
  group_by(Date, Country) %>%
  do(interpolate_yields(., target_maturities)) %>%
  ungroup()

# 5. Transformation en format "Wide" (Matrice Date x Maturité)
df_matrix_form <- df_interpolated %>%
  pivot_wider(names_from = Maturity, values_from = Yield, names_prefix = "Mat_")

# Affichage des maturités pour vérifier qu'elles correspondent à l'article
print(names(df_matrix_form))

# --- VISUALISATION RAPIDE ---
# La Forme de la Courbe (Snapshots)
# On sélectionne quelques dates clés pour ne pas surcharger le graphique
# Par exemple : Première date, une date au milieu, et la dernière date
dates_dispo <- unique(df_interpolated$Date)
selected_dates <- dates_dispo[c(1, floor(length(dates_dispo)/2), length(dates_dispo))]

p_curves <- df_interpolated %>%
  filter(Date %in% selected_dates) %>%
  ggplot(aes(x = Maturity, y = Yield, group = Date, color = as.factor(Date))) +
  geom_line(size = 1) +
  geom_point(size = 2) + # Ajoute les points pour bien voir la "grille" Diebold-Li
  facet_wrap(~Country, ncol = 1) +
  
  labs(
    title = "Formes de la Courbe des Taux (Snapshots)",
    subtitle = "Comparaison de la structure par terme à différentes dates",
    x = "Maturité (Années)",
    y = "Rendement (%)",
    color = "Date"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = target_maturities) + # Affiche bien nos maturités fixes
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

print(p_curves)

# --- VISUALISATION DE L'INTERPOLATION LINEAIRE ---
# --- DONNÉES DU 01/12/2025 ---
# Point 1 : OAT 25/11/2028 (Ligne 36 FR)
x1 <- 2.9842574
y1 <- 2.410
label1 <- "Obs 1: ~2.98 ans"

# Point 2 : OAT 25/02/2029 (Ligne 14 FR)
x2 <- 3.2361396
y2 <- 2.459
label2 <- "Obs 2: ~3.24 ans"

# Cible Diebold-Li
target_x <- 3.00

# --- CALCUL DE L'INTERPOLATION ---
slope <- (y2 - y1) / (x2 - x1)
target_y <- y1 + slope * (target_x - x1)

# Création d'un petit dataframe pour ggplot
df_demo <- data.frame(
  x = c(x1, x2, target_x),
  y = c(y1, y2, target_y),
  type = c("Observé", "Observé", "Interpolé")
)

# --- CRÉATION DU GRAPHIQUE PÉDAGOGIQUE ---
p_demo <- ggplot() +
  
  # 1. La ligne d'interpolation (segment entre les deux points réels)
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
               linetype = "dashed", color = "black", alpha = 0.6, size = 1) +
  
  # 2. Les lignes de projection (pointillés rouges vers les axes)
  geom_segment(aes(x = target_x, xend = target_x, y = 2.38, yend = target_y), 
               linetype = "dotted", color = "#D55E00", size = 0.8) +
  geom_segment(aes(x = 2.9, xend = target_x, y = target_y, yend = target_y), 
               linetype = "dotted", color = "#D55E00", size = 0.8) +
  
  # 3. Les Points (Bleus pour réels, Rouge Losange pour cible)
  geom_point(data = df_demo %>% filter(type == "Observé"), 
             aes(x = x, y = y), color = "#0072B2", size = 5) +
  geom_point(data = df_demo %>% filter(type == "Interpolé"), 
             aes(x = x, y = y), color = "#D55E00", size = 6, shape = 18) +
  
  # 4. Annotations Textuelles (Labels)
  # Label Point 1
  annotate("text", x = x1, y = y1 - 0.005, label = paste0(y1, "%"), 
           color = "#0072B2", vjust = 1, fontface = "bold") +
  # Label Point 2
  annotate("text", x = x2, y = y2 + 0.005, label = paste0(y2, "%"), 
           color = "#0072B2", vjust = 0, fontface = "bold") +
  # Label Cible
  annotate("text", x = target_x, y = target_y + 0.008, 
           label = paste0("Cible 3 ans\n", round(target_y, 3), "%"), 
           color = "#D55E00", fontface = "bold", vjust = 0) +
  
  # 5. Flèche explicative
  annotate("text", x = 3.1, y = 2.40, label = "Hypothèse de linéarité", fontface = "italic", color = "gray30") +
  geom_curve(aes(x = 3.08, y = 2.405, xend = 3.05, yend = 2.42), 
             arrow = arrow(length = unit(0.03, "npc")), curvature = 0.2, color = "gray30") +
  
  # 6. Mise en forme
  scale_x_continuous(limits = c(2.9, 3.3), breaks = seq(2.9, 3.3, 0.1)) +
  scale_y_continuous(limits = c(2.38, 2.48)) +
  labs(
    title = "Principe de l'Interpolation Linéaire (Diebold-Li)",
    subtitle = "Exemple concret : Estimation du taux 3 ans (France, 01/12/2025)",
    x = "Maturité (Années)",
    y = "Rendement (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_demo)


# --- PHASE 2 : CALIBRAGE DU LAMBDA ---

# 1. Calcul de la Courbe Moyenne (Yield Moyen par Maturité)
# On prend l'ensemble des données interpolées pour avoir une "forme typique"
df_avg_curve <- df_interpolated %>%
  group_by(Maturity) %>%
  summarise(Avg_Yield = mean(Yield, na.rm = TRUE))

# 2. Fonction de Coût (RMSE)
# Cette fonction prend un lambda, estime les betas par OLS, et renvoie l'erreur
ns_rmse_function <- function(lambda, data) {
  
  # Calcul des "regressors" (les charges factorielles) pour ce lambda
  tau <- data$Maturity
  
  # Loading Pente (Slope) : (1 - exp(-lambda*tau)) / (lambda*tau)
  x2 <- (1 - exp(-lambda * tau)) / (lambda * tau)
  
  # Loading Courbure (Curvature) : x2 - exp(-lambda*tau)
  x3 <- x2 - exp(-lambda * tau)
  
  # Régression linéaire : Yield ~ 1 + Slope_Factor + Curvature_Factor
  # Le "1" (l'intercept) correspond au Beta1 (Niveau)
  model <- lm(data$Avg_Yield ~ x2 + x3)
  
  # On retourne la racine de l'erreur quadratique moyenne (RMSE)
  return(sqrt(mean(residuals(model)^2)))
}

# 3. Optimisation
# On cherche le lambda entre 0.01 et 2 qui minimise l'erreur
opt_result <- optimize(ns_rmse_function, interval = c(0.01, 3.0), data = df_avg_curve)
best_lambda <- opt_result$minimum

cat("Le Lambda optimal pour nos données est :", round(best_lambda, 4), "\n")
cat("Cela place le pic de courbure à une maturité de :", round(1.7932 / best_lambda, 2), "ans.\n")

# --- GRAPHIQUE 1 : LES FACTOR LOADINGS (Indispensable) ---
# On trace la forme des coefficients pour le lambda optimal

# Création d'une séquence de maturités pour que le tracé soit lisse
plot_mats <- seq(0, 30, length.out = 100)

# Calcul des loadings
loadings_df <- data.frame(
  Maturity = plot_mats,
  Level = 1, # Beta 1 est toujours à 1
  Slope = (1 - exp(-best_lambda * plot_mats)) / (best_lambda * plot_mats),
  Curvature = ((1 - exp(-best_lambda * plot_mats)) / (best_lambda * plot_mats)) - exp(-best_lambda * plot_mats)
) %>%
  pivot_longer(cols = c("Level", "Slope", "Curvature"), names_to = "Factor", values_to = "Loading")

p_loadings <- ggplot(loadings_df, aes(x = Maturity, y = Loading, color = Factor, linetype = Factor)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = 1.7932 / best_lambda, linetype = "dotted", color = "black") +
  annotate("text", x = (1.7932 / best_lambda) + 1, y = 0.1, 
           label = paste0("Max Courbure\n~", round(1.7932 / best_lambda, 2), " ans"), 
           hjust = 0, size = 3.5) +
  scale_color_manual(values = c("Curvature" = "#D55E00", "Level" = "#009E73", "Slope" = "#0072B2")) +
  labs(
    title = "Factor Loadings (Charges Factorielles) de Nelson-Siegel",
    subtitle = paste0("Lambda optimisé = ", round(best_lambda, 3)),
    x = "Maturité (Années)",
    y = "Poids (Loading)"
  ) +
  theme_minimal()

print(p_loadings)


# --- GRAPHIQUE 2 : SENSIBILITÉ DU RMSE (Illustratif / Pédagogique) ---
# Ce graph montre pourquoi on a choisi CE lambda (le creux de la vague)

lambda_seq <- seq(0.1, 2.0, by = 0.05)
rmse_seq <- sapply(lambda_seq, ns_rmse_function, data = df_avg_curve)

p_optim <- data.frame(Lambda = lambda_seq, RMSE = rmse_seq) %>%
  ggplot(aes(x = Lambda, y = RMSE)) +
  geom_line(color = "darkgrey", size = 1) +
  geom_point(data = data.frame(Lambda = best_lambda, RMSE = opt_result$objective), 
             color = "red", size = 4) +
  labs(
    title = "Optimisation du Paramètre Lambda",
    subtitle = "Recherche de la valeur minimisant l'erreur sur la courbe moyenne",
    x = "Valeur de Lambda",
    y = "Erreur (RMSE en %)"
  ) +
  theme_minimal()

print(p_optim)

# --- GRAPHIQUE 3 : Plot de la Courbe moyenne ---
# --- ÉTAPE 1 : Préparation des données pour le graphique ---

# 1. La Courbe Empirique (Moyenne des données réelles)
# (Déjà calculé dans le bloc précédent, mais je le remets pour être autonome)
df_avg_curve <- df_interpolated %>%
  group_by(Maturity) %>%
  summarise(Avg_Yield = mean(Yield, na.rm = TRUE)) %>%
  mutate(Type = "Données Moyennes (Empirique)")

# 2. La Courbe Théorique (Nelson-Siegel avec le Lambda optimal)
# On réutilise les bêtas obtenus par la régression sur la moyenne
# Rappel de la régression faite dans l'optimisation :
x2 <- (1 - exp(-best_lambda * df_avg_curve$Maturity)) / (best_lambda * df_avg_curve$Maturity)
x3 <- x2 - exp(-best_lambda * df_avg_curve$Maturity)
model_avg <- lm(df_avg_curve$Avg_Yield ~ x2 + x3)

# On calcule les points de la courbe ajustée (Fitted)
df_fitted_curve <- data.frame(
  Maturity = df_avg_curve$Maturity,
  Avg_Yield = predict(model_avg),
  Type = "Modèle Nelson-Siegel (Ajusté)"
)

# 3. Fusion pour affichage
df_plot_calibration <- rbind(df_avg_curve, df_fitted_curve)


# --- ÉTAPE 2 : Tracé du Graphique (Diebold-Li Fig. 4) ---

p_calibration <- ggplot() +
  
  # 1. La Courbe Nelson-Siegel (Le Modèle)
  # On trace une ligne continue rouge/orange comme dans l'article
  geom_line(data = df_fitted_curve, 
            aes(x = Maturity, y = Avg_Yield, color = "Modèle Nelson-Siegel"), 
            size = 1.2) +
  
  # 2. Les Points (Les Données Réelles)
  # On affiche uniquement les points, sans les relier, pour bien montrer la "grille"
  geom_point(data = df_avg_curve, 
             aes(x = Maturity, y = Avg_Yield, color = "Données Réelles (Moyenne)"), 
             size = 3, shape = 16, alpha=0.75) +
  
  # 3. Habillage (Style Académique)
  scale_color_manual(values = c("Modèle Nelson-Siegel" = "#D55E00", 
                                "Données Réelles (Moyenne)" = "black")) +
  
  labs(
    title = "Ajustement du Modèle Nelson-Siegel sur la Courbe Moyenne",
    x = "Maturité (Années)",
    y = "Rendement (%)",
    color = NULL # Enlève le titre de la légende pour plus de propreté
  ) +
  
  theme_minimal() +
  theme(
    legend.position = c(0.75, 0.2), # Légende flottante en bas à droite
    legend.background = element_rect(fill = "white", color = "grey90"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank() # Enlève les petites grilles pour plus de clarté
  )

# Courbe NS plot en linéaire par morceau (car 1 image de la fonction par maturité de la grille).
# Particulièrement visible sur les grandes maturités car faible densité de points.
print(p_calibration)



# --- PHASE 3 : Calcul des Séries Temp des Betas ---
# --- ÉTAPE 1 : Estimation Date par Date (Identique à avant) ---

estimate_betas <- function(sub_df, lambda) {
  tau <- sub_df$Maturity
  
  # Calcul des charges factorielles (Loadings)
  # Level = 1
  # Slope = (1 - exp(-lambda * tau)) / (lambda * tau)
  x2 <- (1 - exp(-lambda * tau)) / (lambda * tau)
  # Curvature = Slope - exp(-lambda * tau)
  x3 <- x2 - exp(-lambda * tau)
  
  # Régression OLS : Yield ~ 1 + Slope + Curvature
  model <- lm(sub_df$Yield ~ x2 + x3)
  coeffs <- coef(model)
  
  return(data.frame(
    Beta1_Level = coeffs["(Intercept)"],
    Beta2_Slope = coeffs["x2"],
    Beta3_Curvature = coeffs["x3"],
    R2 = summary(model)$r.squared
  ))
}

# Application de l'estimation
df_factors <- df_interpolated %>%
  group_by(Date, Country) %>%
  do(estimate_betas(., lambda = best_lambda)) %>%
  ungroup() %>%
  arrange(Country, Date)

# --- ÉTAPE 2 : Proxies Empiriques CORRIGÉS (Adaptation 30 ans) ---

# On utilise Mat_30 (le vrai long terme de nos données) au lieu de Mat_10.
# (Car Diebold-Li n'avaient que des données jusqu'à 10Y)
df_empirical <- df_matrix_form %>%
  mutate(
    # Le Niveau est comparé au taux à 30 ans
    # Car Beta1 est l'asymptote (taux infini), donc le taux le plus long est le meilleur proxy.
    Empirical_Level = Mat_30,                        
    
    # La Pente est alors définie comme (Taux Long - Taux Court)
    # Pour être cohérent, on prend aussi 30 ans comme borne longue.
    # Note : Diebold-Li utilisent (10Y - 3M), ici on adapte à nos données.
    Empirical_Slope = Mat_30 - Mat_0.25,             
    
    # La Courbure reste le "Papillon" autour du moyen terme (2 ans)
    Empirical_Curvature = 2 * Mat_2 - (Mat_0.25 + Mat_30) 
  ) %>%
  select(Date, Country, Empirical_Level, Empirical_Slope, Empirical_Curvature)

# Fusion
df_comparison <- left_join(df_factors, df_empirical, by = c("Date", "Country"))

# --- ÉTAPE 3 : Visualisation ---

df_plot_factors <- df_comparison %>%
  pivot_longer(
    cols = c(starts_with("Beta"), starts_with("Empirical")),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Factor_Type = case_when(
      grepl("Level", Metric) ~ "1. Niveau (Level)",
      grepl("Slope", Metric) ~ "2. Pente (Slope)",
      grepl("Curvature", Metric) ~ "3. Courbure (Curvature)"
    ),
    Source = ifelse(grepl("Beta", Metric), "Estimé (Modèle)", "Empirique (Proxy)")
  ) %>%
  # Inversion du Beta2 pour superposition graphique (car Beta2 = -(Pente))
  mutate(Value = ifelse(Metric == "Beta2_Slope", -Value, Value))

# Le Graphique France
p_france <- df_plot_factors %>%
  filter(Country == "France") %>%
  ggplot(aes(x = Date, y = Value, color = Source, linetype = Source)) +
  geom_line(size = 0.8) +
  facet_wrap(~Factor_Type, ncol = 1, scales = "free_y") + # Facet wrap vertical
  
  scale_color_manual(values = c("Estimé (Modèle)" = "#0055A4", "Empirique (Proxy)" = "black")) + # Bleu France
  scale_linetype_manual(values = c("Estimé (Modèle)" = "solid", "Empirique (Proxy)" = "dashed")) +
  
  labs(
    title = "Dynamique des Facteurs - FRANCE",
    subtitle = "Modèle Nelson-Siegel vs Proxies Empiriques (Ajustés 30 ans)",
    x = "Temps",
    y = "Valeur (%)",
    caption = "Beta2 inversé pour lecture 'Long - Court'"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# La Graphique Allemagne
p_germany <- df_plot_factors %>%
  filter(Country == "Germany") %>%
  ggplot(aes(x = Date, y = Value, color = Source, linetype = Source)) +
  geom_line(size = 0.8) +
  facet_wrap(~Factor_Type, ncol = 1, scales = "free_y") +
  
  scale_color_manual(values = c("Estimé (Modèle)" = "#BA0C2F", "Empirique (Proxy)" = "black")) + # Rouge Allemand (approx)
  scale_linetype_manual(values = c("Estimé (Modèle)" = "solid", "Empirique (Proxy)" = "dashed")) +
  
  labs(
    title = "Dynamique des Facteurs - ALLEMAGNE",
    subtitle = "Modèle Nelson-Siegel vs Proxies Empiriques (Ajustés 30 ans)",
    x = "Temps",
    y = "Valeur (%)",
    caption = "Beta2 inversé pour lecture 'Long - Court'"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Affichage et Sauvegarde
# ATTENTION: On observe un décalage entre les TS empiriques et théoriques.
# C'est normal, c'est dû à la différence entre le beta "infini" et celui à 30Y.
# Puisque les 3 coeff sont estimés conjointement, le shift s'observe sur les 3 courbes.
# En revanche, les corrélations empiriques/théoriques (sauf pour la slope Allemande) montrent que
# la time series suit très convenablement les observations.
print(p_france)
print(p_germany)

# --- VÉRIFICATION ---
# Les corrélations doivent rester très élevées (> 0.9)
correlations <- df_comparison %>%
  group_by(Country) %>%
  summarise(
    Corr_Level = cor(Beta1_Level, Empirical_Level),
    Corr_Slope = cor(Beta2_Slope, Empirical_Slope), 
    Corr_Curvature = cor(Beta3_Curvature, Empirical_Curvature)
  )

print(correlations)


# --- PHASE 4 : Analyse de la performance "In-Sample" ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- ÉTAPE 1 : CALCUL DES RÉSIDUS ---

# Fonction pour recalculer le taux théorique à partir des betas
ns_yield_calc <- function(lambda, beta1, beta2, beta3, maturity) {
  x2 <- (1 - exp(-lambda * maturity)) / (lambda * maturity)
  x3 <- x2 - exp(-lambda * maturity)
  return(beta1 + beta2 * x2 + beta3 * x3)
}

# Calcul des erreurs pour l'ensemble des données
df_performance <- df_interpolated %>%
  left_join(df_factors, by = c("Date", "Country")) %>%
  mutate(
    Fitted = ns_yield_calc(best_lambda, Beta1_Level, Beta2_Slope, Beta3_Curvature, Maturity),
    Residual = Yield - Fitted # C'est ici qu'on fait "Réalité - Modèle"
  )

# --- ÉTAPE 2 : STATISTIQUES (RMSE) ---
# Le "Bulletin de notes" du modèle
rmse_stats <- df_performance %>%
  group_by(Country) %>%
  summarise(
    Mean_Error = mean(Residual),      # Biais moyen
    RMSE = sqrt(mean(Residual^2)),    # Erreur type (doit être faible, ex: < 0.10)
    .groups = 'drop'
  )

print("--- Performance Globale (RMSE en %) ---")
print(rmse_stats)


# --- ÉTAPE 3 : GRAPHIQUE 'SNAPSHOT' (UNE SEULE DATE) ---

# On sélectionne uniquement la DERNIÈRE date disponible
last_date <- max(df_performance$Date)

# A. Les Points Réels pour cette date
df_points_last <- df_performance %>%
  filter(Date == last_date)

# B. La Courbe Lisse Modélisée pour cette date (Grille fine 0 à 30 ans)
df_curve_last <- data.frame()
for(ctry in unique(df_factors$Country)) {
  # Récupérer les betas de ce pays à la dernière date
  betas <- df_factors %>% filter(Date == last_date, Country == ctry)
  
  if(nrow(betas) > 0) {
    # Grille fine pour un tracé lisse
    mats_fine <- seq(0.25, 30, by = 0.1)
    yields_fine <- ns_yield_calc(best_lambda, betas$Beta1_Level, betas$Beta2_Slope, betas$Beta3_Curvature, mats_fine)
    
    df_curve_last <- rbind(df_curve_last, data.frame(Date = last_date, Country = ctry, Maturity = mats_fine, Yield = yields_fine))
  }
}

# Le Graphique
p_snapshot_unique <- ggplot() +
  # La Ligne (Le Modèle)
  geom_line(data = df_curve_last, aes(x = Maturity, y = Yield, color = Country), size = 1.2) +
  # Les Points (La Réalité)
  geom_point(data = df_points_last, aes(x = Maturity, y = Yield, color = Country), size = 3, alpha = 0.6) +
  
  scale_color_manual(values = c("France" = "#0055A4", "Germany" = "#BA0C2F")) +
  
  labs(
    title = paste("Qualité de l'Ajustement (Snapshot au", format(last_date, "%d/%m/%Y"), ")"),
    subtitle = "Points = Taux Réels | Ligne = Courbe Nelson-Siegel",
    x = "Maturité (Années)",
    y = "Rendement (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_snapshot_unique)


# --- ÉTAPE 4 : HEATMAP DES RÉSIDUS ---

p_residuals <- ggplot(df_performance, aes(x = Date, y = as.factor(Maturity), fill = Residual)) +
  geom_tile() +
  facet_wrap(~Country, ncol = 1) +
  
  # Gradient : Rouge = Modèle trop bas, Bleu = Modèle trop haut
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, name = "Erreur (%)") +
  
  labs(
    title = "Analyse des Résidus (Heatmap)",
    subtitle = "Visualisation des erreurs de modélisation dans le temps et par maturité",
    x = "Temps",
    y = "Maturité"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Lecture: Chaque groupe de 5 colonnes verticales = 1 semaine.
# Chaque carré est la corrélation entre les taux fixes interpolés en phase 1 
# et le modèle NS fit en phase 3.
print(p_residuals)


# --- PHASE 5 : Modélisation dynamique ---
# On s'assure que les données sont triées chronologiquement
df_factors <- df_factors %>% arrange(Country, Date)

# --- PRÉPARATION : FORMAT LONG ---
# C'est beaucoup plus simple de travailler avec une colonne "Factor" et une colonne "Value"
df_long_factors <- df_factors %>%
  pivot_longer(cols = starts_with("Beta"), names_to = "Factor", values_to = "Value")

# --- ÉTAPE 0 : TESTS DE STATIONNARITÉ (ADF & KPSS) ---

# 1. Installation et Chargement du package indispensable
if(!require(tseries)) install.packages("tseries")
library(tseries)

# 2. On redéfinit la fonction de test (sans le nettoyage excessif car vos données sont propres)
run_stationarity_tests <- function(series) {
  
  # Sécurité de base
  if(length(series) < 10) return(data.frame(ADF_p_value = NA, KPSS_p_value = NA))
  
  # Test ADF (Augmented Dickey-Fuller)
  # H0: Non-Stationnaire. Si p < 0.05 -> Stationnaire.
  # On enlève le suppressWarnings pour voir s'il y a un souci
  adf_res <- try(adf.test(series, alternative = "stationary"), silent = TRUE)
  pval_adf <- if(inherits(adf_res, "try-error")) NA else adf_res$p.value
  
  # Test KPSS
  # H0: Stationnaire. Si p < 0.05 -> Non-Stationnaire.
  kpss_res <- try(kpss.test(series, null = "Level"), silent = TRUE)
  pval_kpss <- if(inherits(kpss_res, "try-error")) NA else kpss_res$p.value
  
  return(data.frame(
    ADF_p_value = pval_adf,
    KPSS_p_value = pval_kpss
  ))
}

# 3. Exécution sur vos données format Long
stationarity_results <- df_long_factors %>%
  group_by(Country, Factor) %>%
  reframe(run_stationarity_tests(Value)) %>%
  mutate(
    Verdict_ADF = ifelse(ADF_p_value < 0.05, "Stationnaire", "Non-Stationnaire"),
    Verdict_KPSS = ifelse(KPSS_p_value > 0.05, "Stationnaire", "Non-Stationnaire"),
    Synthese = case_when(
      is.na(ADF_p_value) ~ "Erreur de calcul",
      Verdict_ADF == "Stationnaire" & Verdict_KPSS == "Stationnaire" ~ "SÛR: Stationnaire",
      Verdict_ADF == "Non-Stationnaire" & Verdict_KPSS == "Non-Stationnaire" ~ "SÛR: Non-Stationnaire",
      TRUE ~ "Ambigu (Mixte)"
    )
  )

print("--- RÉSULTATS FINAUX DES TESTS ---")
print(stationarity_results)
# --- ÉTAPE 1 : ANALYSE DE L'AUTOCORRÉLATION (ACF) DES FACTEURS ---

# Fonction robuste pour calculer l'ACF et retourner un DataFrame
calc_acf_df <- function(series, lags = 20) {
  acf_res <- acf(series, plot = FALSE, lag.max = lags)
  data.frame(
    Lag = as.numeric(acf_res$lag),
    ACF = as.numeric(acf_res$acf)
  ) %>% filter(Lag > 0) # On enlève le lag 0 (toujours 1)
}

# Calcul de l'ACF pour chaque facteur de chaque pays
df_acf_factors <- df_long_factors %>%
  group_by(Country, Factor) %>%
  reframe(calc_acf_df(Value)) # reframe remplace summarise pour les sorties multiples

# GRAPHIQUE 1 : Autocorrélogrammes des Facteurs
p_acf_factors <- ggplot(df_acf_factors, aes(x = Lag, y = ACF, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~Factor, ncol = 1) +
  # Seuil de significativité approximatif (Bartlett)
  geom_hline(yintercept = c(-1.96/sqrt(nrow(df_factors)/2), 1.96/sqrt(nrow(df_factors)/2)), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Mémoire des Facteurs (ACF Brute)",
    subtitle = "Une persistance élevée justifie l'utilisation de modèles auto-régressifs",
    x = "Retard (Jours)",
    y = "Autocorrélation"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("France" = "#0055A4", "Germany" = "#BA0C2F"))

print(p_acf_factors)


# --- ÉTAPE 2 : ESTIMATION DES MODÈLES AR(1) ---
# Modèle : Beta(t) = c + phi * Beta(t-1) + epsilon(t)

fit_ar1 <- function(series) {
  # Estimation du modèle ARIMA(1,0,0) = AR(1)
  model <- arima(series, order = c(1, 0, 0))
  # On capture les résultats importants
  list(
    ar_coef = as.numeric(model$coef["ar1"]), 
    intercept = as.numeric(model$coef["intercept"]),
    residuals = as.numeric(residuals(model))
  )
}

# Application du modèle à chaque série
models_results <- df_long_factors %>%
  group_by(Country, Factor) %>%
  summarise(
    Model_Out = list(fit_ar1(Value)), 
    .groups = "drop"
  ) %>%
  mutate(
    Phi = sapply(Model_Out, function(x) x$ar_coef),
    Intercept = sapply(Model_Out, function(x) x$intercept),
    Residuals = lapply(Model_Out, function(x) x$residuals)
  )

# Tableau des coefficients pour le rapport
print("--- Coefficients AR(1) Estimés ---")
print(models_results %>% select(Country, Factor, Phi, Intercept))


# --- ÉTAPE 3 : DIAGNOSTIC (ACF DES RÉSIDUS) ---
# Vérification que l'AR(1) a bien capturé toute la dynamique

# On extrait et calcule l'ACF des résidus
df_residuals_acf <- models_results %>%
  select(Country, Factor, Residuals) %>%
  unnest(Residuals) %>%
  group_by(Country, Factor) %>%
  reframe(calc_acf_df(Residuals))

# GRAPHIQUE 2 : Autocorrélogrammes des Résidus
p_acf_residuals <- ggplot(df_residuals_acf, aes(x = Lag, y = ACF, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~Factor, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = c(-1.96/sqrt(nrow(df_factors)/2), 1.96/sqrt(nrow(df_factors)/2)), 
             linetype = "dashed", color = "blue") +
  labs(
    title = "Qualité des Modèles AR(1) : ACF des Résidus",
    subtitle = "Absence d'autocorrélation significative = Modèle validé",
    x = "Retard (Jours)",
    y = "Autocorrélation Résiduelle"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("France" = "#0055A4", "Germany" = "#BA0C2F"))

print(p_acf_residuals)











################################################################################
################################################################################
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################################################################################
################################################################################

# DO NOT RUN YET (NOT ENOUGH DATA FOR CONSISTENT INSIGHTS)
# --- PHASE 6 : Prédiction et Comparaison avec Random Walk ---
# 1. Préparation des données de test
n_days_test <- 5
all_dates <- unique(df_factors$Date) %>% sort()
cut_off_date <- all_dates[length(all_dates) - n_days_test]

# On passe df_factors en format long dès le début pour faciliter les calculs
df_factors_long <- df_factors %>%
  pivot_longer(cols = starts_with("Beta"), names_to = "Factor", values_to = "Value")

df_train_long <- df_factors_long %>% filter(Date <= cut_off_date)
df_test_long  <- df_factors_long %>% filter(Date > cut_off_date)

# 2. Prévision des facteurs (AR1 vs RW)
df_forecast_factors <- df_test_long %>%
  left_join(models_results %>% select(Country, Factor, Phi, Intercept), 
            by = c("Country", "Factor")) %>%
  group_by(Country, Factor) %>%
  arrange(Date) %>%
  mutate(
    # Valeur du jour précédent pour la prévision
    # Pour le premier jour du test, on va chercher la dernière valeur du train
    Beta_Prev = lag(Value),
    # Remplissage de la première valeur NA du lag avec la fin du train
    Beta_Prev = ifelse(is.na(Beta_Prev), 
                       df_train_long$Value[df_train_long$Country == first(Country) & 
                                             df_train_long$Factor == first(Factor) & 
                                             df_train_long$Date == cut_off_date], 
                       Beta_Prev),
    
    # Prévision Diebold-Li (AR1)
    Beta_Pred_DL = Intercept + Phi * (Beta_Prev - Intercept),
    # Prévision Naïve (Random Walk)
    Beta_Pred_RW = Beta_Prev
  ) %>%
  ungroup()

# 3. Reconstruction des courbes de taux (Yields)
# On repasse en format large pour avoir les 3 bêtas sur la même ligne
df_preds_wide <- df_forecast_factors %>%
  select(Date, Country, Factor, Beta_Pred_DL, Beta_Pred_RW) %>%
  pivot_wider(names_from = Factor, values_from = c(Beta_Pred_DL, Beta_Pred_RW))

df_final_forecast <- df_interpolated %>%
  filter(Date > cut_off_date) %>%
  left_join(df_preds_wide, by = c("Date", "Country")) %>%
  mutate(
    # Reconstruction via Nelson-Siegel
    Yield_Pred_DL = ns_yield_calc(best_lambda, 
                                  Beta_Pred_DL_Beta1_Level, 
                                  Beta_Pred_DL_Beta2_Slope, 
                                  Beta_Pred_DL_Beta3_Curvature, Maturity),
    Yield_Pred_RW = ns_yield_calc(best_lambda, 
                                  Beta_Pred_RW_Beta1_Level, 
                                  Beta_Pred_RW_Beta2_Slope, 
                                  Beta_Pred_RW_Beta3_Curvature, Maturity)
  )

# 4. Calcul des performances
forecast_perf <- df_final_forecast %>%
  group_by(Country) %>%
  summarise(
    RMSE_DieboldLi = sqrt(mean((Yield - Yield_Pred_DL)^2, na.rm = TRUE)),
    RMSE_RandomWalk = sqrt(mean((Yield - Yield_Pred_RW)^2, na.rm = TRUE))
  )

print("--- PERFORMANCES DE PRÉVISION (OUT-OF-SAMPLE) ---")
print(forecast_perf)

# 5. Visualisation (Exemple Taux 10 ans)
p_forecast_compare <- df_final_forecast %>%
  filter(Maturity == 10) %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = Yield, color = "Réel"), linewidth = 1.2) +
  geom_line(aes(y = Yield_Pred_DL, color = "Diebold-Li (AR1)"), linetype = "dashed") +
  geom_line(aes(y = Yield_Pred_RW, color = "Random Walk"), linetype = "dotted") +
  facet_wrap(~Country, scales = "free_y") +
  labs(
    title = "Prévision du Taux à 10 ans (Hors-Échantillon)",
    subtitle = "Comparaison des modèles sur les 5 derniers jours",
    y = "Yield (%)", x = "Date", color = "Modèle"
  ) +
  theme_minimal()

print(p_forecast_compare)
