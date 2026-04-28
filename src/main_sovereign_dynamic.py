import os
import numpy as np
import matplotlib.pyplot as plt

from data_loader import load_sovereign_data
from preprocessing import (
    prepare_curve_data,
    get_available_dates_for_country,
    filter_country_date_min_obs,
    merge_isin_series,
    filter_nominal_only,
    normalize_isin_column,
)
from data_loader import load_sovereign_data, load_isin_mapping

from mean_functions import fit_nelson_siegel, nelson_siegel_yield, fit_nelson_siegel_grid_lambda, fit_nelson_siegel_multi_start

def plot_dynamic_nss_fit(X, y, tau_grid, y_nss_grid, country_code, date_str):
    plt.figure(figsize=(9, 5))
    plt.scatter(X.flatten(), y, label="Observed yields")
    plt.plot(tau_grid, y_nss_grid, label="Nelson-Siegel fit")
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"Sovereign Dynamic Curve - {country_code} - {date_str}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_dynamic_nss_for_date(df, country_code, date_str, show_plot=True):
    df_curve = filter_country_date_min_obs(
        df=df,
        country=country_code,
        date=date_str,
        min_obs=15
    )

    if len(df_curve) == 0:
        print(f"{country_code} | {date_str} | skipped (not enough observations)")
        return None

    X, y, df_curve = prepare_curve_data(df_curve)
    tau = X.flatten()

    print("Longest maturities:")
    print(
    df_curve[["ISIN", "Maturity", "Residual_Maturity", "Yield"]]
    .sort_values("Residual_Maturity", ascending=False)
    .head(10)
)
    print()
    fit_result = fit_nelson_siegel_grid_lambda(tau, y)
    y_fitted = fit_result["y_fitted"]
    rmse = np.sqrt(np.mean((y - y_fitted) ** 2))

    # Truncated fit to diagnose long-end influence
    fit_result_35, rmse_35, df_fit_35 = fit_nss_with_optional_tau_cap(df_curve, tau_cap=35.0)
    fit_result_40, rmse_40, df_fit_40 = fit_nss_with_optional_tau_cap(df_curve, tau_cap=40.0)

    tau_grid = np.linspace(tau.min(), tau.max(), 300)
    y_nss_grid = nelson_siegel_yield(
        tau_grid,
        fit_result["beta0"],
        fit_result["beta1"],
        fit_result["beta2"],
        fit_result["lambda"]
    )
    tau_grid_35 = np.linspace(tau.min(), tau.max(), 300)
    y_nss_grid_35 = nelson_siegel_yield(
    tau_grid_35,
    fit_result_35["beta0"],
    fit_result_35["beta1"],
    fit_result_35["beta2"],
    fit_result_35["lambda"]
)

    print("=" * 80)
    print(f"Country: {country_code} | Date: {date_str}")
    print("=" * 80)
    print(df_curve.head())
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
    print(f"NSS RMSE: {rmse:.6f}")
    print("Truncated-curve diagnostics:")
    print(f"NSS RMSE (tau <= 35): {rmse_35:.6f} | n_obs={len(df_fit_35)}")
    print(f"NSS RMSE (tau <= 40): {rmse_40:.6f} | n_obs={len(df_fit_40)}")
    print()
    print()

    if show_plot:
        plot_dynamic_nss_fit(X, y, tau_grid, y_nss_grid, country_code, date_str)
        plot_dynamic_nss_comparison(
            X,
            y,
            tau_grid,
            y_nss_grid,
            tau_grid_35,
            y_nss_grid_35,
            country_code,
            date_str,
            tau_cap=35
        )

    return {
        "country": country_code,
        "date": date_str,
        "n_obs": len(df_curve),
        "beta0": fit_result["beta0"],
        "beta1": fit_result["beta1"],
        "beta2": fit_result["beta2"],
        "lambda": fit_result["lambda"],
        "rmse": rmse,
    }

def fit_nss_with_optional_tau_cap(df_curve, tau_cap=None):
    """
    Fit NSS on full curve or truncated curve up to tau_cap.
    """
    df_fit = df_curve.copy()

    if tau_cap is not None:
        df_fit = df_fit[df_fit["Residual_Maturity"] <= tau_cap].copy()

    X_fit, y_fit, df_fit = prepare_curve_data(df_fit)
    tau_fit = X_fit.flatten()

    fit_result = fit_nelson_siegel_grid_lambda(tau_fit, y_fit)
    rmse = np.sqrt(np.mean((y_fit - fit_result["y_fitted"]) ** 2))

    return fit_result, rmse, df_fit

def plot_dynamic_nss_comparison(X, y, tau_grid_full, y_nss_full, tau_grid_trunc, y_nss_trunc, country_code, date_str, tau_cap):
    plt.figure(figsize=(10, 5))
    plt.scatter(X.flatten(), y, label="Observed yields")
    plt.plot(tau_grid_full, y_nss_full, label="NSS full fit")
    plt.plot(tau_grid_trunc, y_nss_trunc, label=f"NSS fit (tau <= {tau_cap})")
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"NSS Comparison - {country_code} - {date_str}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(
        base_dir,
        "data",
        "raw",
        "NewData_Souverain_CrossTempo_Long_Clean.csv"
    )

    country_code = "FR"

    df = load_sovereign_data(file_path)

    mapping_path = os.path.join(base_dir, "data", "raw", "souverains2(1).xlsx")
    df_mapping = load_isin_mapping(mapping_path)

    # Normalize ISIN format on both sides before merge
    df = normalize_isin_column(df, isin_col="ISIN")
    df_mapping = normalize_isin_column(df_mapping, isin_col="ISIN")

    df = merge_isin_series(df, df_mapping)

    print("Columns after merge:")
    print(df.columns.tolist())
    print()

    print("Series value counts after merge:")
    print(df["Series"].value_counts(dropna=False).head(20))
    print()

    df = filter_nominal_only(df)

    print("Shape after nominal-only filter:", df.shape)
    print(df.head())

    available_dates = get_available_dates_for_country(df, country_code)
    print(f"Number of available dates for {country_code}: {len(available_dates)}")
    print("First 5 dates:", available_dates[:5])
    print("Last 5 dates:", available_dates[-5:])
    print()

    # On prend quelques dates espacées dans le temps
    selected_dates = [
        str(available_dates[0])[:10],
        str(available_dates[len(available_dates) // 2])[:10],
        str(available_dates[-1])[:10],
    ]

    results = []

    for date_str in selected_dates:
        result = run_dynamic_nss_for_date(
            df=df,
            country_code=country_code,
            date_str=date_str,
            show_plot=True
        )
        if result is not None:
            results.append(result)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for res in results:
        print(
            f"{res['country']} | "
            f"date={res['date']} | "
            f"n_obs={res['n_obs']} | "
            f"NSS RMSE={res['rmse']:.6f}"
        )


if __name__ == "__main__":
    main()