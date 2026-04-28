import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, WhiteKernel, ConstantKernel

def fit_gp_residuals(X, residuals):
    """
    Fit a Gaussian Process on NSS residuals.

    Parameters
    ----------
    X : np.ndarray
        Residual maturities, shape (n, 1)
    residuals : np.ndarray
        NSS residuals, shape (n,)

    Returns
    -------
    gp : fitted GaussianProcessRegressor
    """
        # Default kernel (chosen after comparison): Matérn 3/2
    kernel = (
       ConstantKernel(1.0, (1e-3, 1e3))
        * Matern(length_scale=5.0, length_scale_bounds=(1e-2, 1e2), nu=1.5)
        + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-8, 1e0))
    )

    # Alternative: RBF kernel (uncomment to use)
    # kernel = (
    #     ConstantKernel(1.0, (1e-3, 1e3))
    #     * RBF(length_scale=5.0, length_scale_bounds=(1e-2, 1e2))
    #     + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-8, 1e0))
    # )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=10,
        random_state=0
    )

    gp.fit(X, residuals)
    return gp


def predict_gp_residuals(gp, X_grid):
    """
    Predict GP residuals on a grid.

    Parameters
    ----------
    gp : GaussianProcessRegressor
        Fitted GP model.
    X_grid : np.ndarray
        Prediction grid, shape (m, 1)

    Returns
    -------
    mean_pred : np.ndarray
    std_pred : np.ndarray
    """
    mean_pred, std_pred = gp.predict(X_grid, return_std=True)
    return mean_pred, std_pred


def fit_gp_residuals_with_kernel(X, residuals, kernel_type="matern52"):
    """
    Fit GP with different kernel choices.
    """
    if kernel_type == "rbf":
        kernel = ConstantKernel(1.0) * RBF(length_scale=1.0) + WhiteKernel()
    elif kernel_type == "matern32":
        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=1.5) + WhiteKernel()
    elif kernel_type == "matern52":
        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5) + WhiteKernel()
    else:
        raise ValueError(f"Unknown kernel_type: {kernel_type}")

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=5
    )

    gp.fit(X, residuals)

    return gp

def fit_gp_residuals_multivariate(X_multi, residuals):
    """
    Fit GP on multivariate inputs, e.g. [Residual_Maturity, Rating_Score].
    """
    kernel = (
        ConstantKernel(1.0, (1e-3, 1e3))
        * Matern(length_scale=[5.0, 2.0], length_scale_bounds=(1e-2, 1e2), nu=1.5)
        + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-8, 1e0))
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=10,
        random_state=0
    )

    gp.fit(X_multi, residuals)
    return gp

