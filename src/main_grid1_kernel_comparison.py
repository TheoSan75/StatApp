import os
import numpy as np
import pandas as pd

from data_loader import load_grid1_sheet
from preprocessing import (
    filter_static_nominal_bonds,
    compute_residual_maturity_from_reference_date,
    prepare_static_curve_data,
)
from mean_functions import fit_nelson_siegel, nelson_siegel_yield
from gp_model import predict_gp_residuals, fit_gp_residuals_with_kernel


def make_blocked_folds_by_maturity(tau, n_splits=5):
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


def run_kernel_comparison(file_path, country_code, reference_date, kernels, n_splits=5):
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

    results = []

    print("=" * 80)
    print(f"Country: {country_code}")
    print("=" * 80)

    for kernel_type in kernels:
        rmses = []

        for train_idx, test_idx in folds:
            X_train, y_train = X[train_idx], y[train_idx]
            X_test, y_test = X[test_idx], y[test_idx]

            tau_train = X_train.flatten()
            tau_test = X_test.flatten()

            # NSS
            fit_result = fit_nelson_siegel(tau_train, y_train)

            y_nss_train = fit_result["y_fitted"]
            residuals_train = y_train - y_nss_train

            y_nss_test = nelson_siegel_yield(
                tau_test,
                fit_result["beta0"],
                fit_result["beta1"],
                fit_result["beta2"],
                fit_result["lambda"]
            )

            # GP
            gp = fit_gp_residuals_with_kernel(X_train, residuals_train, kernel_type=kernel_type)

            gp_mean_test, _ = predict_gp_residuals(gp, X_test)
            y_pred = y_nss_test + gp_mean_test

            rmses.append(rmse(y_test, y_pred))

        avg_rmse = np.mean(rmses)

        print(f"{kernel_type} | mean RMSE: {avg_rmse:.6f}")

        results.append({
            "country": country_code,
            "kernel": kernel_type,
            "rmse": avg_rmse
        })

    print()

    return pd.DataFrame(results)


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "raw", "grid1_uzd01ebt.xlsx")

    reference_date = "2025-01-02"
    countries = ["FR", "DE", "IT"]

    kernels = ["rbf", "matern32", "matern52"]

    all_results = []

    for country in countries:
        df_res = run_kernel_comparison(
            file_path,
            country,
            reference_date,
            kernels
        )
        all_results.append(df_res)

    final_df = pd.concat(all_results, ignore_index=True)

    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(final_df.pivot(index="country", columns="kernel", values="rmse"))


if __name__ == "__main__":
    main()