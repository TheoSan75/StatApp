import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C, WhiteKernel

class HybridGaussianProcessCurve:
    """
    Modèle Hybride : Régression Linéaire (Tendance) + Processus Gaussien Matérn 5/2 (Résidus).
    Version simple et stable pour obtenir des premiers résultats cohérents.
    """
    
    def __init__(self):
        # 1. Modèle déterministe (La "moyenne" demandée par l'encadrant)
        self.lr_model = LinearRegression()
        
        # 2. Modèle probabiliste (Le GP)
        # On utilise le noyau Matérn 5/2 (nu=2.5) pour éviter les zigzags du RBF
        kernel = C(100.0, (1.0, 1e4)) * Matern(length_scale=[10.0, 5.0], 
                                               length_scale_bounds=(3.0, 30.0), 
                                               nu=2.5) \
                 + WhiteKernel(noise_level=10.0, noise_level_bounds=(1e-1, 100.0))
        
        self.gp_model = GaussianProcessRegressor(kernel=kernel, 
                                                 n_restarts_optimizer=10, 
                                                 normalize_y=False,
                                                 random_state=42)
        self.is_fitted = False

    def prepare_features(self, df):
        """Prépare les données X (TTM et Rating)"""
        df_clean = df.dropna(subset=['TTM', 'Rating_Num']).copy()
        X = df_clean[['TTM', 'Rating_Num']].values
        return X, df_clean

    def fit(self, df_train):
        """Entraînement en deux temps : Tendance puis correction GP"""
        X_train, df_train_clean = self.prepare_features(df_train)
        y_train = df_train_clean['Spread_bps'].values
        
        # Étape 1 : On calcule la droite de tendance
        self.lr_model.fit(X_train, y_train)
        y_base_pred = self.lr_model.predict(X_train)
        
        # Étape 2 : On calcule l'écart (résidus)
        y_residuals = y_train - y_base_pred
        
        # Étape 3 : Le GP apprend à corriger la droite
        self.gp_model.fit(X_train, y_residuals)
        
        self.is_fitted = True
        print(f"✅ Modèle Hybride prêt.")

    def predict(self, df_input, return_std=False):
        """Prédiction : Tendance + Ajustement local"""
        if not self.is_fitted:
            raise Exception("Modèle non entraîné.")
            
        X_input, _ = self.prepare_features(df_input)
        y_base = self.lr_model.predict(X_input)
        
        if return_std:
            y_res_pred, sigma = self.gp_model.predict(X_input, return_std=True)
            return y_base + y_res_pred, sigma
        else:
            y_res_pred = self.gp_model.predict(X_input, return_std=False)
            return y_base + y_res_pred