import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

from gp_model import fit_gp_residuals, predict_gp_residuals


def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


def coverage_rate(y_true, y_pred, y_std, z=1.96):
    """
    Empirical coverage of predictive intervals.
    Returns the share of observations lying inside y_pred ± z * y_std.
    """
    lower = y_pred - z * y_std
    upper = y_pred + z * y_std
    inside = (y_true >= lower) & (y_true <= upper)
    return np.mean(inside)


def make_blocked_folds_by_maturity(tau, n_splits=5):
    """
    Sort observations by maturity and split into contiguous blocks.
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


def build_ridge_pipeline():
    numeric_features = ["Residual_Maturity"]
    categorical_features = ["Cntry of Risk", "BICS Level 1", "BBG Composite"]

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_features),
            ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
        ]
    )

    model = Pipeline(
        steps=[
            ("preprocessor", preprocessor),
            ("ridge", Ridge(alpha=1.0))
        ]
    )

    return model


def plot_test_predictions_with_uncertainty(tau_test, y_test, y_pred, y_std, fold_id, selected_date):
    """
    Plot observed test yields, predictive mean, and 95% predictive intervals.
    """
    order = np.argsort(tau_test.flatten())
    tau_sorted = tau_test.flatten()[order]
    y_test_sorted = y_test[order]
    y_pred_sorted = y_pred[order]
    y_std_sorted = y_std[order]

    plt.figure(figsize=(9, 5))
    plt.scatter(tau_sorted, y_test_sorted, label="Observed test yields")
    plt.plot(tau_sorted, y_pred_sorted, label="Ridge + GP prediction")
    plt.fill_between(
        tau_sorted,
        y_pred_sorted - 1.96 * y_std_sorted,
        y_pred_sorted + 1.96 * y_std_sorted,
        alpha=0.2,
        label="95% predictive interval"
    )
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"Corporate test uncertainty - Fold {fold_id} - {selected_date}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def run_corporate_predictive_cv(file_path, selected_date="2026-01-02", n_splits=5, show_uncertainty_plot=True):
    df = pd.read_csv(file_path)

    # -------------------------
    # CLEANING
    # -------------------------
    df = df[(df["Residual_Maturity"] > 0) & (df["Residual_Maturity"] <= 50)]
    df = df[(df["Yield"] > 0) & (df["Yield"] <= 10)]

    # -------------------------
    # DATE FILTER
    # -------------------------
    df = df[df["Date"] == selected_date].copy()

    features = [
        "Residual_Maturity",
        "Cntry of Risk",
        "BICS Level 1",
        "BBG Composite"
    ]
    target = "Yield"

    df = df.dropna(subset=features + [target]).reset_index(drop=True)

    X = df[features].copy()
    y = df[target].to_numpy()
    tau = df["Residual_Maturity"].to_numpy()

    folds = make_blocked_folds_by_maturity(tau, n_splits=n_splits)

    fold_results = []

    print("=" * 90)
    print(f"Corporate predictive validation | Date: {selected_date} | Blocked CV with {n_splits} folds")
    print("=" * 90)
    print(f"Number of observations: {len(df)}")
    print()

    for fold_id, (train_idx, test_idx) in enumerate(folds, start=1):
        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx].copy()

        y_train = y[train_idx]
        y_test = y[test_idx]

        tau_train = X_train["Residual_Maturity"].to_numpy().reshape(-1, 1)
        tau_test = X_test["Residual_Maturity"].to_numpy().reshape(-1, 1)

        # -------------------------
        # BASELINE: RIDGE
        # -------------------------
        ridge_model = build_ridge_pipeline()
        ridge_model.fit(X_train, y_train)

        y_ridge_train = ridge_model.predict(X_train)
        y_ridge_test = ridge_model.predict(X_test)

        ridge_test_rmse = rmse(y_test, y_ridge_test)

        # -------------------------
        # GP ON RIDGE RESIDUALS
        # -------------------------
        residuals_train = y_train - y_ridge_train

        gp = fit_gp_residuals(tau_train, residuals_train)

        gp_mean_test, gp_std_test = predict_gp_residuals(gp, tau_test)
        y_ridge_gp_test = y_ridge_test + gp_mean_test

        ridge_gp_test_rmse = rmse(y_test, y_ridge_gp_test)
        ridge_gp_coverage = coverage_rate(y_test, y_ridge_gp_test, gp_std_test, z=1.96)
        mean_pred_std = np.mean(gp_std_test)

        fold_results.append({
            "fold": fold_id,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "tau_test_min": float(np.min(tau_test)),
            "tau_test_max": float(np.max(tau_test)),
            "ridge_test_rmse": ridge_test_rmse,
            "ridge_gp_test_rmse": ridge_gp_test_rmse,
            "ridge_gp_coverage": ridge_gp_coverage,
            "mean_pred_std": mean_pred_std,
            "gp_kernel": str(gp.kernel_),
        })

        print(f"Fold {fold_id}")
        print(f"  Test maturity block: [{np.min(tau_test):.2f}, {np.max(tau_test):.2f}]")
        print(f"  Ridge test RMSE:      {ridge_test_rmse:.6f}")
        print(f"  Ridge + GP test RMSE: {ridge_gp_test_rmse:.6f}")
        print(f"  95% interval coverage: {ridge_gp_coverage:.3f}")
        print(f"  Mean predictive std:   {mean_pred_std:.6f}")
        print(f"  GP kernel: {gp.kernel_}")
        print()

        # Plot uncertainty only on fold 1 by default
        if show_uncertainty_plot and fold_id == 1:
            plot_test_predictions_with_uncertainty(
                tau_test=tau_test,
                y_test=y_test,
                y_pred=y_ridge_gp_test,
                y_std=gp_std_test,
                fold_id=fold_id,
                selected_date=selected_date
            )

    results_df = pd.DataFrame(fold_results)

    print("-" * 90)
    print("AVERAGE RESULTS")
    print("-" * 90)
    print(f"Mean Ridge test RMSE:       {results_df['ridge_test_rmse'].mean():.6f}")
    print(f"Mean Ridge + GP test RMSE:  {results_df['ridge_gp_test_rmse'].mean():.6f}")
    print(f"Mean 95% interval coverage: {results_df['ridge_gp_coverage'].mean():.3f}")
    print(f"Mean predictive std:        {results_df['mean_pred_std'].mean():.6f}")
    print()

    return results_df


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")

    selected_date = "2026-01-02"
    n_splits = 5

    results_df = run_corporate_predictive_cv(
        file_path=file_path,
        selected_date=selected_date,
        n_splits=n_splits,
        show_uncertainty_plot=True
    )

    print("=" * 90)
    print("FINAL RESULTS TABLE")
    print("=" * 90)
    print(results_df)


if __name__ == "__main__":
    main()