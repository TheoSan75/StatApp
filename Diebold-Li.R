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