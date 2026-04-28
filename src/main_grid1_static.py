import os
import numpy as np
import matplotlib.pyplot as plt

from data_loader import load_grid1_sheet
from preprocessing import (
    filter_static_nominal_bonds,
    compute_residual_maturity_from_reference_date,
    prepare_static_curve_data,
)
from mean_functions import fit_nelson_siegel, nelson_siegel_yield
from gp_model import fit_gp_residuals, predict_gp_residuals


def plot_curve_with_nss(X, y, tau_grid, y_nss_grid, country_code, reference_date):
    plt.figure(figsize=(9, 5))
    plt.scatter(X.flatten(), y, label="Observed yields")
    plt.plot(tau_grid, y_nss_grid, label="Nelson-Siegel fit")
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"Static Sovereign Curve - {country_code} - {reference_date}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_residuals(tau, residuals, country_code, reference_date):
    plt.figure(figsize=(9, 4))
    plt.scatter(tau, residuals, label="NSS residuals")
    plt.axhline(0.0, linestyle="--", label="Zero line")
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Residual (%)")
    plt.title(f"NSS Residuals - {country_code} - {reference_date}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_nss_gp_fit(X, y, tau_grid, y_nss_grid, y_total_grid, std_gp, country_code, reference_date):
    plt.figure(figsize=(10, 5))
    plt.scatter(X.flatten(), y, label="Observed yields")
    plt.plot(tau_grid, y_nss_grid, label="Nelson-Siegel fit")
    plt.plot(tau_grid, y_total_grid, label="NSS + GP fit")

    plt.fill_between(
        tau_grid,
        y_total_grid - 1.96 * std_gp,
        y_total_grid + 1.96 * std_gp,
        alpha=0.2,
        label="95% CI (GP part)"
    )

    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"Sovereign Curve - NSS + GP - {country_code} - {reference_date}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_static_pipeline_for_country(file_path, country_code, reference_date, show_plots=True):
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

    # NSS fit
    fit_result = fit_nelson_siegel(tau, y)
    y_fitted = fit_result["y_fitted"]
    residuals = y - y_fitted
    nss_rmse = np.sqrt(np.mean((y - y_fitted) ** 2))

    # GP on residuals
    gp = fit_gp_residuals(X, residuals)

    # Prediction grid
    tau_grid = np.linspace(tau.min(), tau.max(), 300)
    X_grid = tau_grid.reshape(-1, 1)

    y_nss_grid = nelson_siegel_yield(
        tau_grid,
        fit_result["beta0"],
        fit_result["beta1"],
        fit_result["beta2"],
        fit_result["lambda"]
    )

    gp_mean_grid, gp_std_grid = predict_gp_residuals(gp, X_grid)
    y_total_grid = y_nss_grid + gp_mean_grid

    # In-sample fitted values
    gp_mean_train, _ = predict_gp_residuals(gp, X)
    gp_residual_rmse = np.sqrt(np.mean((residuals - gp_mean_train) ** 2))
    y_nss_gp_train = y_fitted + gp_mean_train
    nss_gp_rmse = np.sqrt(np.mean((y - y_nss_gp_train) ** 2))

    print("=" * 80)
    print(f"Country: {country_code} | Reference date: {reference_date}")
    print("=" * 80)
    print(df_curve[[
        "Issuer Name",
        "Series",
        "Maturity",
        "Residual_Maturity",
        "Mid Yield to Convention"
    ]].head())
    print()
    print(f"Number of observations: {len(df_curve)}")
    print()
    print("Nelson-Siegel parameters:")
    print(f"beta0  = {fit_result['beta0']:.6f}")
    print(f"beta1  = {fit_result['beta1']:.6f}")
    print(f"beta2  = {fit_result['beta2']:.6f}")
    print(f"lambda = {fit_result['lambda']:.6f}")
    print(f"Optimization success: {fit_result['success']}")
    print(f"Loss (SSE): {fit_result['loss']:.6f}")
    print(f"NSS RMSE: {nss_rmse:.6f}")
    print()
    print("Fitted GP kernel:")
    print(gp.kernel_)
    print(f"NSS + GP RMSE: {nss_gp_rmse:.6f}")
    print(f"GP residual RMSE: {gp_residual_rmse:.6f}")
    print()

    if show_plots:
        plot_curve_with_nss(X, y, tau_grid, y_nss_grid, country_code, reference_date)
        plot_residuals(tau, residuals, country_code, reference_date)
        plot_nss_gp_fit(
            X,
            y,
            tau_grid,
            y_nss_grid,
            y_total_grid,
            gp_std_grid,
            country_code,
            reference_date
        )

    return {
        "country": country_code,
        "reference_date": reference_date,
        "n_obs": len(df_curve),
        "beta0": fit_result["beta0"],
        "beta1": fit_result["beta1"],
        "beta2": fit_result["beta2"],
        "lambda": fit_result["lambda"],
        "nss_rmse": nss_rmse,
        "nss_gp_rmse": nss_gp_rmse,
        "gp_kernel": str(gp.kernel_),
        "gp_residual_rmse": gp_residual_rmse,
    }


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "raw", "grid1_uzd01ebt.xlsx")

    reference_date = "2025-01-02"
    countries = ["FR", "DE", "IT"]

    results = []

    for country_code in countries:
        result = run_static_pipeline_for_country(
            file_path=file_path,
            country_code=country_code,
            reference_date=reference_date,
            show_plots=True
        )
        results.append(result)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for res in results:
        print(
            f"{res['country']} | "
            f"n_obs={res['n_obs']} | "
            f"NSS RMSE={res['nss_rmse']:.6f} | "
            f"NSS+GP RMSE={res['nss_gp_rmse']:.6f}"
        )


if __name__ == "__main__":
    main()