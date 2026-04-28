import os
import pandas as pd
import numpy as np

from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

from gp_model import fit_gp_residuals, predict_gp_residuals


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")

    df = pd.read_csv(file_path)

    # -------------------------
    # CLEANING
    # -------------------------
    df = df[(df["Residual_Maturity"] > 0) & (df["Residual_Maturity"] <= 50)]
    df = df[(df["Yield"] > 0) & (df["Yield"] <= 10)]

    # -------------------------
    # SELECT DATE
    # -------------------------
    selected_date = "2026-01-02"
    df = df[df["Date"] == selected_date].copy()

    print("=" * 80)
    print(f"Corporate GP - Date: {selected_date}")
    print("=" * 80)
    print(f"Number of observations: {len(df)}")
    print()

    # -------------------------
    # FEATURES
    # -------------------------
    features = [
        "Residual_Maturity",
        "Cntry of Risk",
        "BICS Level 1",
        "BBG Composite"
    ]

    target = "Yield"

    df = df.dropna(subset=features + [target])

    X = df[features]
    y = df[target].values

    tau = df["Residual_Maturity"].values.reshape(-1, 1)

    # -------------------------
    # BASELINE (Ridge)
    # -------------------------
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

    model.fit(X, y)
    y_pred = model.predict(X)

    ridge_rmse = np.sqrt(np.mean((y - y_pred) ** 2))

    # -------------------------
    # RESIDUALS
    # -------------------------
    residuals = y - y_pred

    # -------------------------
    # GP ON RESIDUALS
    # -------------------------
    gp = fit_gp_residuals(tau, residuals)

    # in-sample prediction
    gp_mean, _ = predict_gp_residuals(gp, tau)
    y_total = y_pred + gp_mean

    gp_rmse = np.sqrt(np.mean((y - y_total) ** 2))

    # -------------------------
    # OUTPUT
    # -------------------------
    print(f"Ridge RMSE: {ridge_rmse:.6f}")
    print()
    print("Fitted GP kernel:")
    print(gp.kernel_)
    print()
    print(f"Ridge + GP RMSE: {gp_rmse:.6f}")


if __name__ == "__main__":
    main()