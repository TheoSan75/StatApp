import numpy as np
from scipy.optimize import minimize
import pandas as pd
from scipy.interpolate import make_smoothing_spline
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
from scipy.optimize import least_squares

class NelsonSiegelCurve:
    """
    Implémentation du modèle Nelson-Siegel.
    Permet de calibrer les paramètres (beta0, beta1, beta2, lambda) 
    sur des données de marché et de prédire des taux.
    """
    
    def __init__(self):
        # Paramètres initialisés à None
        self.beta0 = None
        self.beta1 = None
        self.beta2 = None
        self.lam = None
        self.params = None
        self.rmse = None # Root Mean Square Error du dernier fit

    def formula(self, t, b0, b1, b2, lam):
        """
        La formule mathématique brute de Nelson-Siegel.
        t : Time to Maturity (en années)
        """
        # Sécurité pour éviter la division par zéro si t=0
        t = np.maximum(t, 1e-6)
        
        term1 = (1 - np.exp(-lam * t)) / (lam * t)
        term2 = term1 - np.exp(-lam * t)
        
        return b0 + (b1 * term1) + (b2 * term2)

    def fit(self, t_obs, y_obs):
        """
        Calibre le modèle sur les données observées (Optimisation).
        t_obs : Array des maturités (TTM)
        y_obs : Array des rendements observés (Yields)
        """
        # 1. Définition de la fonction d'erreur (Somme des Carrés)
        def objective(params):
            b0, b1, b2, lam = params
            if lam <= 0: return 1e10 # Pénalité si lambda est négatif ou nul
            
            y_pred = self.formula(t_obs, b0, b1, b2, lam)
            return np.sum((y_pred - y_obs) ** 2)

        # 2. Points de départ (Guess initial)
        # b0 (Long terme) : moyenne des taux longs
        # b1 (Court terme) : écart taux court - taux long
        # b2 (Courbure) : petit par défaut
        # lambda : 0.5 est un standard de départ (pic autour de 2-3 ans)
        initial_guess = [np.mean(y_obs), -2.0, 1.0, 0.5]

        # 3. Optimisation (Algorithme Nelder-Mead ou L-BFGS-B)
        result = minimize(objective, initial_guess, method='L-BFGS-B', 
                          bounds=((0, 15), (-15, 15), (-15, 15), (0.01, 5)))

        # 4. Stockage des résultats
        self.params = result.x
        self.beta0, self.beta1, self.beta2, self.lam = result.x
        
        # Calcul de la performance (RMSE)
        residuals = self.predict(t_obs) - y_obs
        self.rmse = np.sqrt(np.mean(residuals**2))
        
        return self.params

    def predict(self, t_vec):
        """
        Prédit les taux pour une liste de maturités donnée.
        """
        if self.params is None:
            raise Exception("Le modèle n'est pas calibré ! Lancez .fit() d'abord.")
            
        return self.formula(t_vec, self.beta0, self.beta1, self.beta2, self.lam)

    def get_params_dict(self):
        """Retourne les paramètres sous forme de dictionnaire lisible."""
        return {
            "Level (Beta0)": round(self.beta0, 4),
            "Slope (Beta1)": round(self.beta1, 4),
            "Curvature (Beta2)": round(self.beta2, 4),
            "Lambda": round(self.lam, 4),
            "RMSE": round(self.rmse, 5)
        }
    

class SvenssonCurve:
    def __init__(self):
        # Svensson a 6 paramètres : Beta0, Beta1, Beta2, Beta3, Tau1, Tau2
        self.params = None 

    def nss_function(self, t, b0, b1, b2, b3, tau1, tau2):
        # Sécurité pour éviter la division par zéro
        t = np.maximum(t, 1e-6)
        tau1 = np.maximum(tau1, 1e-6)
        tau2 = np.maximum(tau2, 1e-6)
        
        term1 = (1 - np.exp(-t/tau1)) / (t/tau1)
        term2 = term1 - np.exp(-t/tau1)
        term3 = (1 - np.exp(-t/tau2)) / (t/tau2) - np.exp(-t/tau2)
        
        return b0 + b1*term1 + b2*term2 + b3*term3

    def fit(self, t_obs, y_obs):
        # Fonction d'erreur à minimiser (différence entre modèle et réalité)
        def residuals(p):
            return self.nss_function(t_obs, p[0], p[1], p[2], p[3], p[4], p[5]) - y_obs
        
        # Initialisation "intelligente" (b0, b1, b2, b3, tau1, tau2)
        # On part des valeurs typiques de marché
        p0 = [3.0, -1.0, -1.0, 1.0, 2.0, 5.0] 
        
        # Contraintes : Les Tau doivent être positifs
        bounds = ([-np.inf, -np.inf, -np.inf, -np.inf, 0.1, 0.1], 
                  [np.inf, np.inf, np.inf, np.inf, 30.0, 30.0])
        
        res = least_squares(residuals, p0, bounds=bounds, loss='soft_l1')
        self.params = res.x

    def predict(self, t_vec):
        if self.params is None: return np.zeros_like(t_vec)
        # On "déballe" les 6 paramètres pour la fonction
        return self.nss_function(t_vec, *self.params)

    def get_params_dict(self):
        if self.params is None: return {}
        return {
            "Beta0 (Long Terme)": round(self.params[0], 4),
            "Beta1 (Court Terme)": round(self.params[1], 4),
            "Beta2 (Courbure 1)": round(self.params[2], 4),
            "Beta3 (Courbure 2)": round(self.params[3], 4),
            "Tau1 (Position Bosse 1)": round(self.params[4], 4),
            "Tau2 (Position Bosse 2)": round(self.params[5], 4)
        }
    



class RandomForestCurve:
    """
    Modèle de courbe des taux basé sur le Machine Learning (Random Forest).
    Il prédit le SPREAD (vs Swap) et non le Yield brut pour plus de stabilité.
    """
    
    def __init__(self, n_estimators=100, max_depth=None):
        self.model = RandomForestRegressor(n_estimators=n_estimators, 
                                           max_depth=max_depth, 
                                           random_state=42)
        self.encoder = LabelEncoder() # Pour transformer 'France' en 0, 'Italie' en 1
        self.is_fitted = False

    def prepare_features(self, df):
        """
        Transforme le tableau de données en matrice X (Features) pour l'IA.
        On utilise : TTM, Coupon, Amount, Country.
        """
        df_clean = df.copy()
        
        # Encodage du pays (Texte -> Nombre)
        # Note: Dans un vrai projet de prod, il faudrait sauvegarder l'encoder.
        # Ici on le refit à chaque fois pour simplifier l'exemple.
        if 'Country' in df_clean.columns:
            # Si l'encoder n'est pas encore calibré sur des pays
            try:
                df_clean['Country_Code'] = self.encoder.transform(df_clean['Country'])
            except:
                df_clean['Country_Code'] = self.encoder.fit_transform(df_clean['Country'])
        else:
            df_clean['Country_Code'] = 0 # Cas mono-pays
            
        # Sélection des features X
        X = df_clean[['TTM', 'Country_Code', 'Coupon', 'Amount']].fillna(0)
        return X

    def fit(self, df_train):
        """
        Entraîne le modèle.
        df_train : DataFrame contenant les colonnes TTM, Yield, Swap_Rate, Spread_bps...
        """
        X = self.prepare_features(df_train)
        y = df_train['Spread_bps'] # ON APPREND LE SPREAD, PAS LE YIELD !
        
        self.model.fit(X, y)
        self.is_fitted = True
        print(f"Random Forest entraîné sur {len(df_train)} obligations.")

    def predict(self, df_input):
        """
        Prédit le Yield Final.
        Formule : Prediction(Spread) + Taux_Swap_Reel
        """
        if not self.is_fitted:
            raise Exception("Modèle non entraîné.")
            
        X = self.prepare_features(df_input)
        
        # 1. L'IA prédit le spread
        predicted_spread_bps = self.model.predict(X)
        
        # 2. On reconstruit le Yield : Swap + Spread
        # Attention: on divise par 100 car le spread est en bps (points de base)
        final_yield_pred = df_input['Swap_Rate'] + (predicted_spread_bps / 100.0)
        
        return final_yield_pred
    
import xgboost as xgb
import pandas as pd
import numpy as np

class XGBoostCurve:
    def __init__(self, n_estimators=1000, max_depth=3, learning_rate=0.01):
        """
        XGBoost Regressor.
        - n_estimators : Nombre d'arbres (cycles de correction).
        - max_depth : Profondeur max (3 est souvent suffisant pour éviter l'overfitting).
        - learning_rate : Vitesse d'apprentissage (plus c'est petit, plus c'est précis mais lent).
        """
        self.model = xgb.XGBRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            objective='reg:squarederror',
            n_jobs=-1,
            random_state=42,
            # Paramètres anti-overfitting supplémentaires
            subsample=0.8,       # Utilise 80% des données par arbre
            colsample_bytree=0.8 # Utilise 80% des colonnes par arbre
        )

    def fit(self, df_train):
        # Gestion automatique du Country_Code si absent
        if 'Country_Code' not in df_train.columns:
            if 'Country' in df_train.columns:
                df_train = df_train.copy()
                df_train['Country_Code'] = pd.factorize(df_train['Country'])[0]
            else:
                raise ValueError("Il manque la colonne Country ou Country_Code")

        X = df_train[['TTM', 'Coupon', 'Amount', 'Country_Code']]
        y = df_train['Spread_bps']
        
        self.model.fit(X, y)

    def predict(self, df_input):
        # Préparation identique à fit
        if isinstance(df_input, pd.DataFrame):
            X = df_input.copy()
        else:
            return np.zeros(len(df_input))

        if 'Country_Code' not in X.columns and 'Country' in X.columns:
             X['Country_Code'] = pd.factorize(X['Country'])[0]
             
        X_pred = X[['TTM', 'Coupon', 'Amount', 'Country_Code']]
        
        return self.model.predict(X_pred)


class GaussianProcessCurve:
    """
    Modèle probabiliste non-paramétrique.
    Avantages : Lissage naturel + Intervalle de confiance (Incertitude).
    Inconvénients : Plus lent à calculer sur des gros datasets (O(N^3)).
    """
    
    def __init__(self):
        # Définition du "Noyau" (Kernel) - C'est le coeur du moteur GP
        # 1. RBF : Gère la forme lisse (la corrélation entre les points proches)
        # 2. WhiteKernel : Gère le "bruit" du marché (les petites erreurs de prix)
        kernel = C(1.0, (1e-3, 1e3)) * RBF(length_scale=10.0, length_scale_bounds=(1e-2, 1e3)) \
                 + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-10, 1e-1))
        
        # normalize_y=True est CRUCIAL car les spreads sont petits
        self.model = GaussianProcessRegressor(kernel=kernel, 
                                              n_restarts_optimizer=9, 
                                              normalize_y=True)
        self.is_fitted = False

    def prepare_features(self, df):
        """
        Pour le GP, on reste simple : on regarde surtout le TEMPS (TTM).
        On peut ajouter le Pays si on veut un modèle multi-courbes.
        Pour cet exemple, on se concentre sur une courbe 1D (TTM -> Spread).
        """
        # On redimensionne en matrice 2D (N, 1) car sklearn l'exige
        return df['TTM'].values.reshape(-1, 1)

    def fit(self, df_train):
        """Entraînement sur le Spread"""
        X = self.prepare_features(df_train)
        y = df_train['Spread_bps'].values # On apprend le Spread
        
        self.model.fit(X, y)
        self.is_fitted = True
        
        # On peut afficher les hyper-paramètres trouvés par le modèle
        print(f"GP Calibré. Noyau optimal : {self.model.kernel_}")

    def predict(self, df_input, return_std=False):
        """
        Prédit le Yield Final.
        Option : return_std=True renvoie aussi l'incertitude (écart-type).
        """
        if not self.is_fitted:
            raise Exception("Modèle non entraîné.")
            
        X = self.prepare_features(df_input)
        
        # Le GP renvoie le spread ET l'écart-type (sigma) si demandé
        if return_std:
            spread_pred, sigma = self.model.predict(X, return_std=True)
        else:
            spread_pred = self.model.predict(X, return_std=False)
            sigma = None
        
        # Reconstruction du Yield : Swap + Spread
        # spread_pred est en bps, on divise par 100
        yield_pred = df_input['Swap_Rate'].values + (spread_pred / 100.0)
        
        if return_std:
            # On renvoie aussi sigma converti en % (divisé par 100)
            return yield_pred, sigma / 100.0
        
        return yield_pred
    

class CubicSplineCurve:
    def __init__(self, smoothing_factor=1e-3, label="Spline"):
        self.lam = smoothing_factor 
        self.label = label # Pour les différencier dans le graph
        self.model = None

    def fit(self, t_obs, y_obs):
        idx = np.argsort(t_obs)
        t_sorted, y_sorted = t_obs[idx], y_obs[idx]
        df_tmp = pd.DataFrame({'t': t_sorted, 'y': y_sorted}).groupby('t').mean().reset_index()
        self.model = make_smoothing_spline(df_tmp['t'], df_tmp['y'], lam=self.lam)

    def predict(self, t_vec):
        if self.model is None: return np.zeros_like(t_vec)
        return np.clip(self.model(t_vec), -1, 10)
    
