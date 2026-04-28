import pandas as pd
import numpy as np


import pandas as pd
import numpy as np


def filter_country_date(df: pd.DataFrame, country: str, date: str) -> pd.DataFrame:
    target_date = pd.to_datetime(date)

    out = df.copy()
    out = out[out["Country"] == country]
    out = out[out["Date"].dt.normalize() == target_date.normalize()]

    out = out.sort_values("Residual_Maturity").reset_index(drop=True)
    return out

def prepare_curve_data(df_curve: pd.DataFrame):
    """
    Prepare X and y for curve fitting.

    Parameters
    ----------
    df_curve : pd.DataFrame
        Dataframe containing one sovereign curve.

    Returns
    -------
    X : np.ndarray
        Residual maturity as 2D array of shape (n, 1).
    y : np.ndarray
        Yield values as 1D array.
    df_curve : pd.DataFrame
        Cleaned dataframe.
    """
    df_curve = df_curve.copy()

    # Keep only usable rows
    df_curve = df_curve.dropna(subset=["Residual_Maturity", "Yield"])

    # Remove non-positive maturities
    df_curve = df_curve[df_curve["Residual_Maturity"] > 0]

    X = df_curve["Residual_Maturity"].to_numpy().reshape(-1, 1)
    y = df_curve["Yield"].to_numpy()

    return X, y, df_curve

from sklearn.cluster import KMeans
import numpy as np


def filter_nominal_curve_kmeans(X, y, random_state=0):
    """
    Use KMeans clustering to separate nominal vs inflation-linked bonds.

    Parameters
    ----------
    X : np.ndarray
        Residual maturities, shape (n, 1)
    y : np.ndarray
        Yields, shape (n,)
    random_state : int
        For reproducibility

    Returns
    -------
    X_filtered : np.ndarray
    y_filtered : np.ndarray
    mask : np.ndarray (bool)
        Mask of selected observations (nominal curve)
    """

    data = np.hstack([X, y.reshape(-1, 1)])

    kmeans = KMeans(n_clusters=2, random_state=random_state, n_init=10)
    labels = kmeans.fit_predict(data)

    # Compute average yield per cluster
    cluster_means = [
        y[labels == i].mean() for i in range(2)
    ]

    # Nominal curve = higher yield cluster
    nominal_cluster = np.argmax(cluster_means)

    mask = labels == nominal_cluster

    return X[mask], y[mask], mask

def filter_static_nominal_bonds(df: pd.DataFrame, country_code: str) -> pd.DataFrame:
    """
    Filter nominal sovereign bonds in grid1 static sheets using the Series column.
    """

    out = df.copy()

    if "Series" not in out.columns:
        raise ValueError("Column 'Series' not found in dataframe.")

    if country_code == "FR":
        out = out[out["Series"] == "OAT"]

    elif country_code == "IT":
        # Exclude inflation-linked Italian bonds
        out = out[~out["Series"].isin(["CPI", "ICPI"])]

    elif country_code == "DE":
        # Keep all for now; later we can decide how to treat TWIN
        out = out.copy()

    return out.reset_index(drop=True)

def compute_residual_maturity_from_reference_date(
    df: pd.DataFrame,
    reference_date: str
) -> pd.DataFrame:
    """
    Compute residual maturity in years from a reference date.
    """
    out = df.copy()

    ref_date = pd.to_datetime(reference_date)

    out = out.dropna(subset=["Maturity", "Mid Yield to Convention"])
    out["Residual_Maturity"] = (out["Maturity"] - ref_date).dt.days / 365.25
    out = out[out["Residual_Maturity"] > 0]

    return out.sort_values("Residual_Maturity").reset_index(drop=True)

def prepare_static_curve_data(df_curve: pd.DataFrame):
    """
    Prepare X and y from a static sovereign sheet.
    """
    out = df_curve.copy()

    out = out.dropna(subset=["Residual_Maturity", "Mid Yield to Convention"])
    out = out.sort_values("Residual_Maturity").reset_index(drop=True)

    X = out["Residual_Maturity"].to_numpy().reshape(-1, 1)
    y = out["Mid Yield to Convention"].to_numpy()

    return X, y, out

def get_available_dates_for_country(df: pd.DataFrame, country: str):
    """
    Return sorted available dates for one country.
    """
    out = df[df["Country"] == country].copy()
    dates = sorted(out["Date"].dropna().dt.normalize().unique())
    return dates

def filter_country_date_min_obs(df: pd.DataFrame, country: str, date: str, min_obs: int = 10) -> pd.DataFrame:
    """
    Filter one country and one date, and ensure enough observations.
    """
    out = filter_country_date(df, country=country, date=date)
    out = out.dropna(subset=["Residual_Maturity", "Yield"])
    out = out[out["Residual_Maturity"] > 0]
    out = out.sort_values("Residual_Maturity").reset_index(drop=True)

    if len(out) < min_obs:
        return pd.DataFrame(columns=out.columns)

    return out

def merge_isin_series(df_cross, df_mapping):
    """
    Merge CrossTempo data with ISIN -> Series mapping.
    """
    df = df_cross.merge(df_mapping, on="ISIN", how="left")
    return df

def filter_nominal_only(df: pd.DataFrame):
    """
    Keep only nominal bonds (exclude inflation-linked).
    """
    return df[df["Series"] == "OAT"].copy()

def normalize_isin_column(df: pd.DataFrame, isin_col: str = "ISIN") -> pd.DataFrame:
    """
    Normalize ISIN codes by keeping only the first token before any space.
    Example: 'FR0010916924 Govt' -> 'FR0010916924'
    """
    out = df.copy()

    out[isin_col] = (
        out[isin_col]
        .astype(str)
        .str.strip()
        .str.split()
        .str[0]
    )

    return out