import os
import matplotlib.pyplot as plt

from data_loader import load_sovereign_data
from preprocessing import (
    filter_country_date,
    prepare_curve_data,
    filter_nominal_curve_kmeans
)


def plot_observed_curve(X, y, country, date):
    plt.figure(figsize=(9, 5))
    plt.scatter(X.flatten(), y)
    plt.xlabel("Residual Maturity (years)")
    plt.ylabel("Yield (%)")
    plt.title(f"Observed Sovereign Curve - {country} - {date}")
    plt.grid(True)
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

    df = load_sovereign_data(file_path)

    country = "FR"
    date = "2025-01-02"

    df_curve = filter_country_date(df, country=country, date=date)
    X, y, df_curve = prepare_curve_data(df_curve)

    print(f"Selected country: {country}")
    print(f"Selected date: {date}")
    print(f"Number of observations: {len(df_curve)}")
    print()
    print(df_curve.head())
    print()
    print(df_curve.describe())

    if len(df_curve) == 0:
        print("No observations found for this country/date.")
        return

    plot_observed_curve(X, y, country, date)


if __name__ == "__main__":
    main()