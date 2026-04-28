import os
import pandas as pd


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, "data", "raw", "souverains2(1).xlsx")

    xls = pd.ExcelFile(file_path)

    print("Sheets found:")
    print(xls.sheet_names)
    print()

    for sheet in xls.sheet_names:
        print("=" * 100)
        print(f"SHEET: {sheet}")
        print("=" * 100)

        df = pd.read_excel(file_path, sheet_name=sheet)

        print("Columns:")
        print(df.columns.tolist())
        print()

        print("First 5 rows:")
        print(df.head())
        print()

        # Print possible useful columns if they exist
        candidate_cols = [
            "ISIN",
            "Country",
            "Series",
            "Issuer Name",
            "Description",
            "Type",
            "Instrument",
            "Inflation",
            "Linked"
        ]

        existing = [col for col in candidate_cols if col in df.columns]
        if existing:
            print("Potentially useful columns found:")
            print(existing)
            print()

            for col in existing:
                print(f"Unique values preview for column: {col}")
                vals = df[col].dropna().astype(str).unique()[:20]
                print(vals)
                print()

        print("\n")


if __name__ == "__main__":
    main()