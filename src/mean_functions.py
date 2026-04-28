import numpy as np
from scipy.optimize import minimize


def nelson_siegel_yield(tau, beta0, beta1, beta2, lambd):
    """
    Nelson-Siegel yield curve.

    Parameters
    ----------
    tau : array-like
        Residual maturities in years.
    beta0, beta1, beta2, lambd : float
        Nelson-Siegel parameters.

    Returns
    -------
    np.ndarray
        Fitted yields.
    """
    tau = np.asarray(tau, dtype=float)

    # Avoid division by zero
    tau_safe = np.where(tau == 0, 1e-6, tau)

    factor1 = (1 - np.exp(-lambd * tau_safe)) / (lambd * tau_safe)
    factor2 = factor1 - np.exp(-lambd * tau_safe)

    return beta0 + beta1 * factor1 + beta2 * factor2


def nelson_siegel_loss(params, tau, y):
    """
    Sum of squared errors loss for Nelson-Siegel.
    """
    beta0, beta1, beta2, lambd = params
    y_hat = nelson_siegel_yield(tau, beta0, beta1, beta2, lambd)
    return np.sum((y - y_hat) ** 2)


def fit_nelson_siegel(tau, y):
    """
    Fit Nelson-Siegel parameters by nonlinear optimization.

    Parameters
    ----------
    tau : array-like
        Residual maturities in years.
    y : array-like
        Observed yields.

    Returns
    -------
    dict
        Dictionary containing fitted parameters and fitted values.
    """
    tau = np.asarray(tau, dtype=float)
    y = np.asarray(y, dtype=float)

    initial_guess = [
        np.mean(y),   # beta0
        -1.0,         # beta1
        1.0,          # beta2
        0.2           # lambda
    ]

    bounds = [
        (None, None),   # beta0
        (None, None),   # beta1
        (None, None),   # beta2
        (1e-4, 5.0)     # lambda > 0
    ]

    result = minimize(
        nelson_siegel_loss,
        x0=initial_guess,
        args=(tau, y),
        method="L-BFGS-B",
        bounds=bounds
    )

    beta0, beta1, beta2, lambd = result.x
    y_fitted = nelson_siegel_yield(tau, beta0, beta1, beta2, lambd)

    return {
        "beta0": beta0,
        "beta1": beta1,
        "beta2": beta2,
        "lambda": lambd,
        "y_fitted": y_fitted,
        "success": result.success,
        "message": result.message,
        "loss": result.fun
    }

def fit_nelson_siegel_multi_start(tau, y):
    """
    Fit Nelson-Siegel with multiple starting points and keep the best solution.
    """
    tau = np.asarray(tau, dtype=float)
    y = np.asarray(y, dtype=float)

    initial_guesses = [
        [np.mean(y), -1.0, 1.0, 0.05],
        [np.mean(y), -1.0, 1.0, 0.10],
        [np.mean(y), -1.0, 1.0, 0.20],
        [y[-1], y[0] - y[-1], 1.0, 0.10],
        [y[-1], -2.0, 2.0, 0.15],
        [y[-1], -1.0, 0.5, 0.30],
        [np.median(y), -0.5, 0.5, 0.50],
        [np.mean(y), -2.0, 2.0, 1.00],
    ]

    bounds = [
        (-10.0, 10.0),   # beta0
        (-20.0, 20.0),   # beta1
        (-20.0, 20.0),   # beta2
        (1e-3, 3.0)      # lambda
    ]

    best_result = None
    best_loss = np.inf

    for guess in initial_guesses:
        result = minimize(
            nelson_siegel_loss,
            x0=guess,
            args=(tau, y),
            method="L-BFGS-B",
            bounds=bounds
        )

        if result.fun < best_loss:
            best_loss = result.fun
            best_result = result

    beta0, beta1, beta2, lambd = best_result.x
    y_fitted = nelson_siegel_yield(tau, beta0, beta1, beta2, lambd)

    return {
        "beta0": beta0,
        "beta1": beta1,
        "beta2": beta2,
        "lambda": lambd,
        "y_fitted": y_fitted,
        "success": best_result.success,
        "message": best_result.message,
        "loss": best_result.fun
    }

def fit_nelson_siegel_grid_lambda(tau, y, lambda_grid=None):
    """
    Stable Nelson-Siegel fit:
    - grid search over lambda
    - linear least squares for betas conditional on lambda
    """
    tau = np.asarray(tau, dtype=float)
    y = np.asarray(y, dtype=float)

    if lambda_grid is None:
        lambda_grid = np.linspace(0.02, 2.0, 200)

    best_sse = np.inf
    best_params = None
    best_y_fitted = None

    tau_safe = np.where(tau == 0, 1e-6, tau)

    for lambd in lambda_grid:
        factor1 = (1 - np.exp(-lambd * tau_safe)) / (lambd * tau_safe)
        factor2 = factor1 - np.exp(-lambd * tau_safe)

        X = np.column_stack([
            np.ones_like(tau),
            factor1,
            factor2
        ])

        # OLS conditional on lambda
        beta_hat, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        y_hat = X @ beta_hat
        sse = np.sum((y - y_hat) ** 2)

        if sse < best_sse:
            best_sse = sse
            best_params = (beta_hat[0], beta_hat[1], beta_hat[2], lambd)
            best_y_fitted = y_hat

    beta0, beta1, beta2, lambd = best_params

    return {
        "beta0": beta0,
        "beta1": beta1,
        "beta2": beta2,
        "lambda": lambd,
        "y_fitted": best_y_fitted,
        "success": True,
        "message": "Grid search over lambda + OLS for betas",
        "loss": best_sse
    }