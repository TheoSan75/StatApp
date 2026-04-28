import pandas as pd


def load_sovereign_data(path: str) -> pd.DataFrame:
    """
    Load the sovereign dataset from a CSV file.

    Parameters
    ----------
    path : str
        Path to the CSV file.

    Returns
    -------
    pd.DataFrame
        Loaded dataframe with cleaned column names and parsed dates.
    """
    df = pd.read_csv(path)

    # Clean column names
    df.columns = [col.strip() for col in df.columns]

    # Convert date column if present
    if "Date" in df.columns:
        df["Date"] = pd.to_datetime(df["Date"], errors="coerce")

    return df

def load_grid1_sheet(path: str, sheet_name: str) -> pd.DataFrame:
    """
    Load one sovereign static sheet from grid1 workbook.
    """
    df = pd.read_excel(path, sheet_name=sheet_name)
    df.columns = [col.strip() for col in df.columns]

    if "Maturity" in df.columns:
        df["Maturity"] = pd.to_datetime(df["Maturity"], dayfirst=True, errors="coerce")

    return df

def load_isin_mapping(path: str) -> pd.DataFrame:
    """
    Load ISIN -> Series mapping from Feuil3.
    """
    df = pd.read_excel(path, sheet_name="Feuil3")

    df = df.rename(columns={
        "Unnamed: 1": "ISIN",
        "Unnamed: 5": "Series"
    })

    df = df[["ISIN", "Series"]]
    df = df.dropna(subset=["ISIN", "Series"])
    df = df.drop_duplicates(subset=["ISIN"]).reset_index(drop=True)

    return df