rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)

# 1. Chargement ----------------------------------------------------------------
data <- read_csv("r_data/NewData_Souverain_CrossTempo_Long_Clean.csv",
                 show_col_types = FALSE) |>
  mutate(Date = as.Date(Date), Yield = as.numeric(Yield)) |>
  filter(!is.na(Yield))

message(nrow(data), " lignes | ",
        length(unique(data$Country)), " pays | ",
        length(unique(data$ISIN)), " obligations")

# 2. Nombre de bonds observés par pays x date ----------------------------------
bonds_per_date <- data |>
  group_by(Country, Date) |>
  summarise(N_bonds = n_distinct(ISIN), .groups = "drop")

# 3. Statistiques résumées par pays -------------------------------------------
summary_by_country <- bonds_per_date |>
  group_by(Country) |>
  summarise(
    N_dates     = n(),
    Avg_bonds   = round(mean(N_bonds), 1),
    Min_bonds   = min(N_bonds),
    Max_bonds   = max(N_bonds),
    .groups = "drop"
  ) |>
  arrange(desc(Avg_bonds))

message("\n── Résumé par pays ──")
print(summary_by_country)

# 4. Plot par pays -------------------------------------------------------------
countries <- unique(bonds_per_date$Country)

for (country in countries) {
  df_plot <- bonds_per_date |> filter(Country == country)
  
  # Seuil minimum recommandé pour Diebold-Li
  threshold <- 6
  
  p <- ggplot(df_plot, aes(x = Date, y = N_bonds)) +
    geom_line(color = "#378ADD", linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 0.8, color = "#185FA5", alpha = 0.5) +
    # Zone sous le seuil en rouge
    geom_hline(yintercept = threshold, linetype = "dashed",
               color = "#E24B4A", linewidth = 0.8) +
    annotate("text",
             x      = min(df_plot$Date),
             y      = threshold + 0.3,
             label  = paste0("Seuil min. Diebold-Li (", threshold, " bonds)"),
             color  = "#E24B4A", size = 3, hjust = 0, fontface = "italic") +
    # Aire sous le seuil
    geom_ribbon(aes(ymin = pmin(N_bonds, threshold), ymax = threshold),
                fill = "#E24B4A", alpha = 0.15) +
    scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      title    = paste0(country, " — Nombre de bonds observés par date"),
      subtitle = paste0("Moy: ", round(mean(df_plot$N_bonds), 1),
                        " | Min: ", min(df_plot$N_bonds),
                        " | Max: ", max(df_plot$N_bonds)),
      x = NULL,
      y = "Nombre de bonds observés"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
  
  print(p)
}

# 5. Vue globale : tous les pays superposés -----------------------------------
p_global <- ggplot(bonds_per_date, aes(x = Date, y = N_bonds, color = Country)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(
    title    = "Nombre d'obligations observées par date — Tous pays",
    subtitle = "Ligne noire = seuil minimum pour NSS (6 points)",
    x = NULL, y = "Nombre d'obligations observées", color = "Pays"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_global)

# 6. Sauvegarder en PDF -------------------------------------------------------
pdf("r_data/bonds_per_date_by_country.pdf", width = 10, height = 5)
for (country in countries) {
  df_plot <- bonds_per_date |> filter(Country == country)
  p <- ggplot(df_plot, aes(x = Date, y = N_bonds)) +
    geom_line(color = "#378ADD", linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 0.8, color = "#185FA5", alpha = 0.5) +
    geom_hline(yintercept = 6, linetype = "dashed", color = "#E24B4A", linewidth = 0.8) +
    geom_ribbon(aes(ymin = pmin(N_bonds, 6), ymax = 6),
                fill = "#E24B4A", alpha = 0.15) +
    scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      title    = paste0(country, " — Nombre de bonds observés par date"),
      subtitle = paste0("Moy: ", round(mean(df_plot$N_bonds), 1),
                        " | Min: ", min(df_plot$N_bonds),
                        " | Max: ", max(df_plot$N_bonds)),
      x = NULL, y = "Nombre de bonds observés"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )
  print(p)
}
print(p_global)
dev.off()

message("Done. PDF sauvegardé : r_data/bonds_per_date_by_country.pdf")