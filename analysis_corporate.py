import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# =================================================================
# 1. PARAMÈTRES ET MAPPING
# =================================================================
CORPO_FILE = 'NewData_Corpo_clean.csv'
SWAP_FILE = 'grid1_uzd01ebt.xlsx' # Ton fichier contenant l'onglet Swap
TARGET_DATE = 'YLD_16032026'

RATING_MAP = {
    'AAA': 1, 'AA+': 2, 'AA': 3, 'AA-': 4, 'A+': 5, 'A': 6, 'A-': 7,
    'BBB+': 8, 'BBB': 9, 'BBB-': 10, 'BB+': 11, 'BB': 12, 'BB-': 13,
    'B+': 14, 'B': 15, 'B-': 16, 'CCC+': 17, 'CCC': 18, 'CCC-': 19,
    'CC': 20, 'C': 21, 'D': 22
}

def parse_tenor(tenor):
    """Fonction utilitaire pour convertir les durées '6M', '2Y' en années (float)"""
    if pd.isna(tenor): return None
    tenor = str(tenor).strip().upper()
    try:
        if 'M' in tenor: return float(tenor.replace('M', '')) / 12.0
        if 'Y' in tenor: return float(tenor.replace('Y', ''))
        return float(tenor)
    except: return None

def prepare_corporate_dataset():
    print("--- 1. CALIBRATION DU BENCHMARK SWAP ---")
    # Chargement de la courbe de Swap
    df_swap = pd.read_excel(SWAP_FILE, sheet_name='Swap')
    df_swap['TTM'] = df_swap['Tenor'].apply(parse_tenor)
    df_swap = df_swap.dropna(subset=['TTM', 'Yield']).sort_values('TTM').drop_duplicates(subset=['TTM'])
    
    # Création de la fonction d'interpolation continue pour la courbe Swap
    swap_curve = interp1d(df_swap['TTM'], df_swap['Yield'], kind='linear', bounds_error=False, fill_value="extrapolate")
    print("Benchmark Swap calibré et interpolé.")

    print("\n--- 2. TRAITEMENT DE LA BASE CORPORATE ---")
    df_corpo = pd.read_csv(CORPO_FILE)
    
    # Nettoyage de base
    df_corpo = df_corpo.dropna(subset=[TARGET_DATE, 'Residual_Maturity', 'BBG Composite'])
    df_corpo = df_corpo.rename(columns={'Residual_Maturity': 'TTM', TARGET_DATE: 'Yield'})
    
    # Mapping de la note de crédit en valeur numérique
    df_corpo['Rating_Num'] = df_corpo['BBG Composite'].map(RATING_MAP)
    df_corpo = df_corpo.dropna(subset=['Rating_Num']) # Vire les NR (Not Rated)
    
    # Filtrer les obligations à option
    if 'Mty Type' in df_corpo.columns:
        df_corpo = df_corpo[df_corpo['Mty Type'] != 'CALLABLE']

    print(f"Obligations Corporate exploitables : {len(df_corpo)}")

    print("\n--- 3. CALCUL DU SPREAD (CORPO - SWAP) ---")
    # On calcule le taux Swap exact pour chaque obligation corporate
    df_corpo['Swap_Rate'] = swap_curve(df_corpo['TTM'])
    
    # Calcul final du Z-Spread (approximé) en points de base
    df_corpo['Spread_bps'] = (df_corpo['Yield'] - df_corpo['Swap_Rate']) * 100
    
    # Export de la base propre
    df_corpo.to_csv('ML_Ready_Corporate_Dataset.csv', index=False)
    print("Base Corporate exportée avec succès : ML_Ready_Corporate_Dataset.csv")
    
    return df_corpo

if __name__ == "__main__":
    df_ready = prepare_corporate_dataset()