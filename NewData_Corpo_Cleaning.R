library(readxl)
library(dplyr)
library(tidyr)

# ── CONFIG ──────────────────────────────────────────────────────────────────

file_path <- "r_data/NewData_Corpo.xlsx"

# Feuille 1 = Info des bonds
info_sheet  <- "16032026"
price_dates <- c(
  "16032026",
  "13032026",
  "06032026",
  "27022026",
  "20022026",
  "13022026",
  "06022026",
  "30012026",
  "23012026",
  "16012026",
  "09012026",
  "02012026"
)

# Map date label → sheet name in your workbook (edit right-hand side if needed)
sheet_names <- setNames(price_dates, price_dates)

# ── BOND INFO (Feuille 1) ─────────────────────────────────────────────────

bond_info <- read_excel(file_path, sheet = info_sheet) |>
  rename_with(trimws)  # strip any accidental whitespace in headers

# The 16/03/2026 yield is already embedded in Sheet 1 as YLD_CNV_MID_16032026
# Rename it to a standard format so it merges cleanly with the other dates
bond_info <- bond_info |>
  rename(YLD_16032026 = YLD_CNV_MID_16032026)

# ── ADD RESIDUAL_MATURITY COL ────────────────────────────────────────────────

bond_info <- bond_info |>
  mutate(
    Maturity = as.Date(Maturity, format = "%d/%m/%Y"),
    Residual_Maturity = as.numeric(Maturity - Sys.Date()) / 365.25
  )

# ── READ & STACK WEEKLY PRICE SHEETS ────────────────────────────────────────

read_price_sheet <- function(sheet, date_label) {
  df <- read_excel(file_path, sheet = sheet) |>
    rename_with(trimws)
  
  # Each weekly sheet has: date header in row 1, then ID_BB / YLD / ISIN columns
  # read_excel picks up the column names; normalise them
  df <- df |>
    select(ISIN, yield = starts_with("YLD")) |>
    filter(!is.na(ISIN)) |>
    mutate(
      date  = date_label,
      yield = as.numeric(yield)
    )
  df
}

# Skip 16/03/2026 — already in Sheet 1; process the rest
other_dates <- price_dates[price_dates != "16032026"]

price_long <- bind_rows(
  lapply(other_dates, function(d) {
    tryCatch(
      read_price_sheet(sheet_names[[d]], d),
      error = function(e) {
        warning("Could not read sheet for ", d, ": ", e$message)
        NULL
      }
    )
  })
)

# ── PIVOT YIELDS WIDE (one column per date) ──────────────────────────────────

# Clean date label → valid column name: "13/03/2026" → "YLD_13032026"
price_wide <- price_long |>
  mutate(col = paste0("YLD_", gsub("/", "", date))) |>
  select(ISIN, col, yield) |>
  pivot_wider(names_from = col, values_from = yield)

# ── JOIN EVERYTHING ──────────────────────────────────────────────────────────

combined <- bond_info |>
  left_join(price_wide, by = "ISIN")

# Reorder yield columns chronologically (oldest → newest)
yld_cols <- sort(grep("^YLD_", names(combined), value = TRUE))  # lexicographic sort works for DDMMYYYY
non_yld  <- setdiff(names(combined), yld_cols)
combined <- combined |> select(all_of(non_yld), all_of(yld_cols))

# ── EXPORT ──────────────────────────────────────────────────────────────────

output_path <- "r_data/NewData_Corpo_clean.csv"
write.csv(combined, output_path, row.names = FALSE)
message("Done! Saved to: ", output_path)
message(nrow(combined), " bonds × ", ncol(combined), " columns")

# ── CHECK ──────────────────────────────────────────────────────────────────

clean_corpo <- read.csv2("r_data/NewData_Corpo_clean.csv", sep=",")