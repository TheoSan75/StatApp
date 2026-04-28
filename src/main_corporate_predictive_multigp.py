import os
import numpy as np
import pandas as pd

from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

from gp_model import fit_gp_residuals_multivariate, predict_gp_residuals


def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


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


def rating_to_score(series):
    """
    Convert rating strings into ordinal scores.
    Lower score = better credit quality.
    Unknown ratings become NaN.
    """
    rating_map = {
        "AAA": 1,
        "AA+": 2,
        "AA": 3,
        "AA-": 4,
        "A+": 5,
        "A": 6,
        "A-": 7,
        "BBB+": 8,
        "BBB": 9,
        "BBB-": 10,
        "BB+": 11,
        "BB": 12,
        "BB-": 13,
        "B+": 14,
        "B": 15,
        "B-": 16,
        "CCC+": 17,
        "CCC": 18,
        "CCC-": 19,
        "CC": 20,
        "C": 21,
        "D": 22,
    }

    s = series.astype(str).str.strip().str.upper()
    return s.map(rating_map)


def run_corporate_predictive_cv(file_path, selected_date="2026-01-02", n_splits=5):
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

    # Rating score for multivariate GP
    df["Rating_Score"] = rating_to_score(df["BBG Composite"])

    features = [
        "Residual_Maturity",
        "Cntry of Risk",
        "BICS Level 1",
        "BBG Composite"
    ]
    target = "Yield"

    df = df.dropna(subset=features + [target, "Rating_Score"]).reset_index(drop=True)

    X = df[features].copy()
    y = df[target].to_numpy()
    tau = df["Residual_Maturity"].to_numpy()
    rating_score = df["Rating_Score"].to_numpy()

    folds = make_blocked_folds_by_maturity(tau, n_splits=n_splits)

    fold_results = []

    print("=" * 90)
    print(f"Corporate predictive validation with multivariate GP | Date: {selected_date} | Blocked CV with {n_splits} folds")
    print("=" * 90)
    print(f"Number of observations: {len(df)}")
    print()

    for fold_id, (train_idx, test_idx) in enumerate(folds, start=1):
        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx].copy()

        y_train = y[train_idx]
        y_test = y[test_idx]

        tau_train = tau[train_idx]
        tau_test = tau[test_idx]

        rating_train = rating_score[train_idx]
        rating_test = rating_score[test_idx]

        # Ridge baseline
        ridge_model = build_ridge_pipeline()
        ridge_model.fit(X_train, y_train)

        y_ridge_train = ridge_model.predict(X_train)
        y_ridge_test = ridge_model.predict(X_test)

        ridge_test_rmse = rmse(y_test, y_ridge_test)

        residuals_train = y_train - y_ridge_train

        # Multivariate GP input = [maturity, rating_score]
        X_gp_train = np.column_stack([tau_train, rating_train])
        X_gp_test = np.column_stack([tau_test, rating_test])

        gp = fit_gp_residuals_multivariate(X_gp_train, residuals_train)

        gp_mean_test, gp_std_test = predict_gp_residuals(gp, X_gp_test)
        y_ridge_gp_test = y_ridge_test + gp_mean_test

        ridge_gp_test_rmse = rmse(y_test, y_ridge_gp_test)

        fold_results.append({
            "fold": fold_id,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "tau_test_min": float(np.min(tau_test)),
            "tau_test_max": float(np.max(tau_test)),
            "ridge_test_rmse": ridge_test_rmse,
            "ridge_gp_test_rmse": ridge_gp_test_rmse,
            "gp_kernel": str(gp.kernel_),
        })

        print(f"Fold {fold_id}")
        print(f"  Test maturity block: [{np.min(tau_test):.2f}, {np.max(tau_test):.2f}]")
        print(f"  Ridge test RMSE:             {ridge_test_rmse:.6f}")
        print(f"  Ridge + multivariate GP RMSE:{ridge_gp_test_rmse:.6f}")
        print(f"  GP kernel: {gp.kernel_}")
        print()

    results_df = pd.DataFrame(fold_results)

    print("-" * 90)
    print("AVERAGE RESULTS")
    print("-" * 90)
    print(f"Mean Ridge test RMSE:               {results_df['ridge_test_rmse'].mean():.6f}")
    print(f"Mean Ridge + multivariate GP RMSE:  {results_df['ridge_gp_test_rmse'].mean():.6f}")
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
        n_splits=n_splits
    )

    print("=" * 90)
    print("FINAL RESULTS TABLE")
    print("=" * 90)
    print(results_df)


if __name__ == "__main__":
    main()