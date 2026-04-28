import os
import pandas as pd
import numpy as np


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "raw", "NewData_Corpo_clean.csv")

    df = pd.read_csv(file_path)

    print("=" * 100)
    print("RAW DATASET OVERVIEW")
    print("=" * 100)
    print("Shape:", df.shape)
    print()
    print("Columns:")
    print(df.columns.tolist())
    print()
    print("First 5 rows:")
    print(df.head())
    print()

    # Identify yield columns automatically
    yield_cols = [col for col in df.columns if str(col).startswith("YLD_")]
    static_cols = [col for col in df.columns if col not in yield_cols]

    print("=" * 100)
    print("COLUMN SPLIT")
    print("=" * 100)
    print(f"Number of static columns: {len(static_cols)}")
    print(f"Number of yield columns: {len(yield_cols)}")
    print()
    print("Static columns:")
    print(static_cols)
    print()
    print("First 10 yield columns:")
    print(yield_cols[:10])
    print()

    # Basic checks on likely useful fields
    candidate_cols = [
        "Residual_Maturity",
        "Cntry of Risk",
        "BICS Level",
        "BICS Level 1",
        "RTG_MOODY",
        "RTG_SP",
        "RTG_FITCH",
        "BB Composite",
        "Ticker",
        "Issuer Name",
        "ISIN"
    ]

    print("=" * 100)
    print("USEFUL COLUMN PREVIEW")
    print("=" * 100)
    for col in candidate_cols:
        if col in df.columns:
            print(f"\nColumn: {col}")
            print(df[col].dropna().astype(str).unique()[:20])

    # Melt wide -> long
    df_long = df.melt(
        id_vars=static_cols,
        value_vars=yield_cols,
        var_name="Date_Column",
        value_name="Yield"
    )

    # Extract date from columns such as YLD_02012026
    df_long["Date"] = (
        df_long["Date_Column"]
        .str.replace("YLD_", "", regex=False)
    )

    df_long["Date"] = pd.to_datetime(
        df_long["Date"],
        format="%d%m%Y",
        errors="coerce"
    )

    # Basic numeric cleaning
    if "Residual_Maturity" in df_long.columns:
        df_long["Residual_Maturity"] = pd.to_numeric(df_long["Residual_Maturity"], errors="coerce")

    df_long["Yield"] = pd.to_numeric(df_long["Yield"], errors="coerce")

    # Keep only usable rows
    df_long = df_long.dropna(subset=["Date", "Yield"])

    if "Residual_Maturity" in df_long.columns:
        df_long = df_long.dropna(subset=["Residual_Maturity"])
        df_long = df_long[df_long["Residual_Maturity"] > 0]

    print()
    print("=" * 100)
    print("LONG DATASET OVERVIEW")
    print("=" * 100)
    print("Shape:", df_long.shape)
    print()
    print("First 10 rows:")
    print(df_long.head(10))
    print()

    if "Date" in df_long.columns:
        print("Min date:", df_long["Date"].min())
        print("Max date:", df_long["Date"].max())
        print()

    if "Residual_Maturity" in df_long.columns:
        print("Residual maturity summary:")
        print(df_long["Residual_Maturity"].describe())
        print()

    print("Yield summary:")
    print(df_long["Yield"].describe())
    print()

    # Save processed long dataset for later use
    output_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_long.to_csv(output_path, index=False)

    print(f"Saved long dataset to: {output_path}")


if __name__ == "__main__":
    main()