import os
import numpy as np
import pandas as pd

from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

from gp_model import fit_gp_residuals, predict_gp_residuals


def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


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


def build_preprocessor():
    numeric_features = ["Residual_Maturity"]
    categorical_features = ["Cntry of Risk", "BICS Level 1", "BBG Composite"]

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_features),
            ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
        ]
    )

    return preprocessor


def build_ridge_pipeline():
    return Pipeline(
        steps=[
            ("preprocessor", build_preprocessor()),
            ("ridge", Ridge(alpha=1.0))
        ]
    )


def build_rf_pipeline():
    return Pipeline(
        steps=[
            ("preprocessor", build_preprocessor()),
            ("rf", RandomForestRegressor(
                n_estimators=300,
                max_depth=12,
                min_samples_leaf=5,
                random_state=0,
                n_jobs=-1
            ))
        ]
    )


def run_corporate_model_comparison(file_path, selected_date="2026-01-02", n_splits=5):
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
    print(f"Corporate model comparison | Date: {selected_date} | Blocked CV with {n_splits} folds")
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
        # RIDGE
        # -------------------------
        ridge_model = build_ridge_pipeline()
        ridge_model.fit(X_train, y_train)

        y_ridge_train = ridge_model.predict(X_train)
        y_ridge_test = ridge_model.predict(X_test)
        ridge_test_rmse = rmse(y_test, y_ridge_test)

        # -------------------------
        # RIDGE + GP
        # -------------------------
        residuals_train = y_train - y_ridge_train

        gp = fit_gp_residuals(tau_train, residuals_train)
        gp_mean_test, gp_std_test = predict_gp_residuals(gp, tau_test)

        y_ridge_gp_test = y_ridge_test + gp_mean_test
        ridge_gp_test_rmse = rmse(y_test, y_ridge_gp_test)

        # -------------------------
        # RANDOM FOREST
        # -------------------------
        rf_model = build_rf_pipeline()
        rf_model.fit(X_train, y_train)

        y_rf_test = rf_model.predict(X_test)
        rf_test_rmse = rmse(y_test, y_rf_test)

        fold_results.append({
            "fold": fold_id,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "tau_test_min": float(np.min(tau_test)),
            "tau_test_max": float(np.max(tau_test)),
            "ridge_test_rmse": ridge_test_rmse,
            "ridge_gp_test_rmse": ridge_gp_test_rmse,
            "rf_test_rmse": rf_test_rmse,
            "gp_kernel": str(gp.kernel_),
        })

        print(f"Fold {fold_id}")
        print(f"  Test maturity block: [{np.min(tau_test):.2f}, {np.max(tau_test):.2f}]")
        print(f"  Ridge test RMSE:      {ridge_test_rmse:.6f}")
        print(f"  Ridge + GP test RMSE: {ridge_gp_test_rmse:.6f}")
        print(f"  RF test RMSE:         {rf_test_rmse:.6f}")
        print()

    results_df = pd.DataFrame(fold_results)

    print("-" * 90)
    print("AVERAGE RESULTS")
    print("-" * 90)
    print(f"Mean Ridge test RMSE:      {results_df['ridge_test_rmse'].mean():.6f}")
    print(f"Mean Ridge + GP test RMSE: {results_df['ridge_gp_test_rmse'].mean():.6f}")
    print(f"Mean RF test RMSE:         {results_df['rf_test_rmse'].mean():.6f}")
    print()

    return results_df


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")

    selected_date = "2026-01-02"
    n_splits = 5

    results_df = run_corporate_model_comparison(
        file_path=file_path,
        selected_date=selected_date,
        n_splits=n_splits
    )

    print("=" * 90)
    print("FINAL RESULTS TABLE")
    print("=" * 90)
    print(results_df)


if __name__ == "__main__":
    main()