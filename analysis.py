import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# =================================================================
# 1. PARAMÈTRES
# =================================================================
EXCEL_FILE = 'grid1_uzd01ebt.xlsx'
VALUATION_DATE = pd.to_datetime('2024-05-22')
COUNTRIES = ['France', 'Allemagne', 'Italie']

def parse_tenor(tenor):
    if pd.isna(tenor): return None
    tenor = str(tenor).strip().upper()
    try:
        if 'M' in tenor: return float(tenor.replace('M', '')) / 12.0
        if 'Y' in tenor: return float(tenor.replace('Y', ''))
        return float(tenor)
    except: return None

def identify_bond_type(row):
    series = str(row['Series']).strip().upper() if pd.notna(row['Series']) else ""
    name = str(row['Issuer Name']).strip().upper()
    country = row['Country']

    if country == 'France':
        if 'OATI' in series or 'OATE' in series or 'INDEX' in name: return 'Inflation'
        return 'Nominal'
    if country == 'Italie':
        if series in ['CPI', 'ICPI'] or 'BTP ITALIA' in name: return 'Inflation'
        if series in ['FUT', 'VALR', 'PIU']: return 'Retail/Other'
        return 'Nominal'
    if country == 'Allemagne':
        if 'I/L' in series or 'INFLATION' in name: return 'Inflation'
        return 'Nominal'
    return 'Nominal'

def run_pipeline():
    print(f"--- DÉMARRAGE DU PIPELINE (Source: {EXCEL_FILE}) ---")
    
    # 1. Swap Curve
    try:
        df_swap = pd.read_excel(EXCEL_FILE, sheet_name='Swap')
    except:
        print("Erreur: Onglet Swap introuvable.")
        return None
        
    df_swap['TTM'] = df_swap['Tenor'].apply(parse_tenor)
    df_swap = df_swap.dropna(subset=['TTM', 'Yield']).sort_values('TTM').drop_duplicates(subset=['TTM'])
    swap_curve = interp1d(df_swap['TTM'], df_swap['Yield'], kind='linear', bounds_error=False, fill_value="extrapolate")
    print("-> Swap Curve : OK")

    # 2. Pays
    all_data = []
    for country in COUNTRIES:
        try:
            df = pd.read_excel(EXCEL_FILE, sheet_name=country)
            df['Country'] = country
        except: continue

        df = df.rename(columns={'Mid Yield to Convention': 'Yield', 'Amt Out': 'Amount', 'Cpn': 'Coupon'})
        df['Maturity'] = pd.to_datetime(df['Maturity'], dayfirst=True, errors='coerce')
        df['Yield'] = pd.to_numeric(df['Yield'], errors='coerce')
        df['Amount'] = pd.to_numeric(df['Amount'], errors='coerce')
        df['Coupon'] = pd.to_numeric(df['Coupon'], errors='coerce')

        df['TTM'] = (df['Maturity'] - VALUATION_DATE).dt.days / 365.25
        df['Bond_Type'] = df.apply(identify_bond_type, axis=1)
        df = df[df['TTM'] > 0].dropna(subset=['Yield'])
        
        # Calculs
        df['Swap_Rate'] = swap_curve(df['TTM'])
        df['Spread_bps'] = (df['Yield'] - df['Swap_Rate']) * 100

        # === CORRECTION ICI : Ajout de 'Swap_Rate' dans la liste ===
        cols = ['Country', 'Issuer Name', 'Series', 'Bond_Type', 'Maturity', 'TTM', 
                'Yield', 'Swap_Rate', 'Spread_bps', 'Coupon', 'Amount']
        
        # Sécurité pour ne sélectionner que les colonnes qui existent vraiment
        cols = [c for c in cols if c in df.columns]
        
        all_data.append(df[cols])
        print(f"-> {country} : Traité")

    # 3. Export
    df_master = pd.concat(all_data, ignore_index=True)
    df_nominal = df_master[df_master['Bond_Type'] == 'Nominal'].copy()
    
    df_master.to_csv('Master_Dataset_Full.csv', index=False)
    df_nominal.to_csv('ML_Ready_Nominal_Dataset.csv', index=False)
    
    print(f"\nCorrection Appliquée : Fichier généré avec la colonne 'Swap_Rate'.")
    return df_nominal

if __name__ == "__main__":
    run_pipeline()