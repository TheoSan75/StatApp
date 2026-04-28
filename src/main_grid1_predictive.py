import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_loader import load_grid1_sheet
from preprocessing import (
    filter_static_nominal_bonds,
    compute_residual_maturity_from_reference_date,
    prepare_static_curve_data,
)
from mean_functions import fit_nelson_siegel, nelson_siegel_yield
from gp_model import fit_gp_residuals, predict_gp_residuals


def make_blocked_folds_by_maturity(tau, n_splits=5):
    """
    Create blocked folds by sorting maturities and splitting them into contiguous blocks.
    Returns a list of (train_idx, test_idx).
    """
    tau = np.asarray(tau)
    sorted_idx = np.argsort(tau)
    blocks = np.array_split(sorted_idx, n_splits)

    folds = []
    for k in range(n_splits):
        test_idx = blocks[k]
        train_idx = np.concatenate([blocks[j] for j in range(n_splits) if j != k])
        folds.append((train_idx, test_idx))

    return folds


def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


def run_predictive_cv_for_country(file_path, country_code, reference_date, n_splits=5, show_fold_plot=False):
    sheet_map = {
        "FR": "France",
        "DE": "Allemagne",
        "IT": "Italie",
    }

    df = load_grid1_sheet(file_path, sheet_map[country_code])
    df = filter_static_nominal_bonds(df, country_code=country_code)
    df = compute_residual_maturity_from_reference_date(df, reference_date=reference_date)
    X, y, df_curve = prepare_static_curve_data(df)

    tau = X.flatten()
    folds = make_blocked_folds_by_maturity(tau, n_splits=n_splits)

    fold_results = []

    print("=" * 90)
    print(f"Country: {country_code} | Reference date: {reference_date} | Blocked CV with {n_splits} folds")
    print("=" * 90)
    print(f"Number of observations: {len(df_curve)}")
    print()

    for fold_id, (train_idx, test_idx) in enumerate(folds, start=1):
        X_train = X[train_idx]
        y_train = y[train_idx]
        X_test = X[test_idx]
        y_test = y[test_idx]

        tau_train = X_train.flatten()
        tau_test = X_test.flatten()

        # Fit NSS on train only
        fit_result = fit_nelson_siegel(tau_train, y_train)

        y_nss_train = fit_result["y_fitted"]
        residuals_train = y_train - y_nss_train

        # Predict NSS on test
        y_nss_test = nelson_siegel_yield(
            tau_test,
            fit_result["beta0"],
            fit_result["beta1"],
            fit_result["beta2"],
            fit_result["lambda"]
        )

        nss_test_rmse = rmse(y_test, y_nss_test)

        # Fit GP on NSS residuals (train only)
        gp = fit_gp_residuals(X_train, residuals_train)

        gp_residual_test_mean, gp_residual_test_std = predict_gp_residuals(gp, X_test)
        y_nss_gp_test = y_nss_test + gp_residual_test_mean

        nss_gp_test_rmse = rmse(y_test, y_nss_gp_test)

        fold_results.append({
            "fold": fold_id,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "tau_test_min": float(np.min(tau_test)),
            "tau_test_max": float(np.max(tau_test)),
            "nss_test_rmse": nss_test_rmse,
            "nss_gp_test_rmse": nss_gp_test_rmse,
            "beta0": fit_result["beta0"],
            "beta1": fit_result["beta1"],
            "beta2": fit_result["beta2"],
            "lambda": fit_result["lambda"],
            "gp_kernel": str(gp.kernel_),
        })

        print(f"Fold {fold_id}")
        print(f"  Test maturity block: [{np.min(tau_test):.2f}, {np.max(tau_test):.2f}]")
        print(f"  NSS test RMSE:     {nss_test_rmse:.6f}")
        print(f"  NSS+GP test RMSE:  {nss_gp_test_rmse:.6f}")
        print()

        if show_fold_plot:
            tau_grid = np.linspace(tau.min(), tau.max(), 300)
            X_grid = tau_grid.reshape(-1, 1)

            y_nss_grid = nelson_siegel_yield(
                tau_grid,
                fit_result["beta0"],
                fit_result["beta1"],
                fit_result["beta2"],
                fit_result["lambda"]
            )
            gp_grid_mean, gp_grid_std = predict_gp_residuals(gp, X_grid)
            y_nss_gp_grid = y_nss_grid + gp_grid_mean

            plt.figure(figsize=(10, 5))
            plt.scatter(X_train.flatten(), y_train, label="Train points")
            plt.scatter(X_test.flatten(), y_test, label="Test points")
            plt.plot(tau_grid, y_nss_grid, label="NSS fit")
            plt.plot(tau_grid, y_nss_gp_grid, label="NSS + GP fit")
            plt.fill_between(
                tau_grid,
                y_nss_gp_grid - 1.96 * gp_grid_std,
                y_nss_gp_grid + 1.96 * gp_grid_std,
                alpha=0.2,
                label="95% CI"
            )
            plt.xlabel("Residual Maturity (years)")
            plt.ylabel("Yield (%)")
            plt.title(f"{country_code} | Fold {fold_id}")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()

    results_df = pd.DataFrame(fold_results)

    print("-" * 90)
    print("AVERAGE RESULTS")
    print("-" * 90)
    print(f"Mean NSS test RMSE:    {results_df['nss_test_rmse'].mean():.6f}")
    print(f"Mean NSS+GP test RMSE: {results_df['nss_gp_test_rmse'].mean():.6f}")
    print()

    return results_df


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "raw", "grid1_uzd01ebt.xlsx")

    reference_date = "2025-01-02"
    countries = ["FR", "DE", "IT"]
    n_splits = 5

    all_results = []

    for country_code in countries:
        results_df = run_predictive_cv_for_country(
            file_path=file_path,
            country_code=country_code,
            reference_date=reference_date,
            n_splits=n_splits,
            show_fold_plot=False
        )
        results_df["country"] = country_code
        all_results.append(results_df)

    summary_df = pd.concat(all_results, ignore_index=True)

    print("\n" + "=" * 90)
    print("GLOBAL SUMMARY")
    print("=" * 90)

    grouped = summary_df.groupby("country")[["nss_test_rmse", "nss_gp_test_rmse"]].mean()
    print(grouped)


if __name__ == "__main__":
    main()