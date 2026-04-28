import os
import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "processed", "corporate_long.csv")

    df = pd.read_csv(file_path)

    # -------------------------
    # CLEANING
    # -------------------------
    df = df[df["Residual_Maturity"] <= 50]
    df = df[df["Yield"] <= 10]

    # -------------------------
    # SELECT DATE
    # -------------------------
    selected_date = "2026-01-02"
    df = df[df["Date"] == selected_date].copy()

    print("=" * 80)
    print(f"Corporate dataset - Date: {selected_date}")
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

    # -------------------------
    # PREPROCESSING
    # -------------------------
    numeric_features = ["Residual_Maturity"]
    categorical_features = ["Cntry of Risk", "BICS Level 1", "BBG Composite"]

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_features),
            ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
        ]
    )

    # -------------------------
    # MODEL
    # -------------------------
    model = Pipeline(
        steps=[
            ("preprocessor", preprocessor),
            ("ridge", Ridge(alpha=1.0))
        ]
    )

    model.fit(X, y)

    y_pred = model.predict(X)

    rmse = np.sqrt(np.mean((y - y_pred) ** 2))

    print(f"Ridge RMSE: {rmse:.6f}")
    print()

    return model, df


if __name__ == "__main__":
    main()