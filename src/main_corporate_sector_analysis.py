import os
import numpy as np
import pandas as pd

from main_corporate_rf_comparison import run_corporate_model_comparison


def run_sector_analysis(file_path, selected_date="2026-01-02", n_splits=5, min_obs=60):
    df = pd.read_csv(file_path)

    # Nettoyage de base
    df = df[(df["Residual_Maturity"] > 0) & (df["Residual_Maturity"] <= 50)]
    df = df[(df["Yield"] > 0) & (df["Yield"] <= 10)]
    df = df[df["Date"] == selected_date].copy()

    # Comptage par secteur
    sector_counts = df["BICS Level 1"].value_counts()

    print("=" * 90)
    print("SECTOR DISTRIBUTION")
    print("=" * 90)
    print(sector_counts)
    print()

    results = []

    for sector, count in sector_counts.items():
        if count < min_obs:
            continue

        print("=" * 90)
        print(f"SECTOR: {sector} | n_obs={count}")
        print("=" * 90)

        df_sector = df[df["BICS Level 1"] == sector].copy()

        # Sauvegarde temporaire pour réutiliser la pipeline existante
        temp_path = "temp_sector.csv"
        df_sector.to_csv(temp_path, index=False)

        res_df = run_corporate_model_comparison(
            file_path=temp_path,
            selected_date=selected_date,
            n_splits=n_splits
        )

        results.append({
            "sector": sector,
            "n_obs": count,
            "ridge_rmse": res_df["ridge_test_rmse"].mean(),
            "ridge_gp_rmse": res_df["ridge_gp_test_rmse"].mean(),
            "rf_rmse": res_df["rf_test_rmse"].mean(),
        })

    results_df = pd.DataFrame(results)

    print("=" * 90)
    print("SECTOR SUMMARY")
    print("=" * 90)
    print(results_df.sort_values("ridge_gp_rmse"))

    return results_df


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")

    run_sector_analysis(
        file_path=file_path,
        selected_date="2026-01-02",
        n_splits=5,
        min_obs=60
    )


if __name__ == "__main__":
    main()