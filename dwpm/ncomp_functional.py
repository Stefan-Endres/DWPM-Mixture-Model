#%% Critical imports
import copy
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
#from shgo import shgo  # !pip install shgo
from scipy.optimize import shgo
from scipy.optimize import linprog
import numpy as np
from scipy.optimize import linprog, shgo
from scipy.optimize import NonlinearConstraint

#%% Equilibrium solvers
import numpy as np
from scipy.optimize import linprog, shgo, NonlinearConstraint

##############################################################################
# 1) LBD function
##############################################################################
def lbd(x, g_func, Lambda, Z_0, g_func_args=()):
    """
    Evaluate the lower-bounding function for the dual problem.

    Parameters
    ----------
    x : 1D array_like
        Composition (length n). Must satisfy sum(x)=1, min_x <= x[i] <= max_x.
    g_func : callable
        Gibbs free energy function G(x), returning scalar.
    Lambda : 1D array_like
        Lagrange multipliers (length n).
    Z_0 : 1D array_like
        Feed composition (length n), sum(Z_0)=1.

    Returns
    -------
    float
        G(x) + Lambda dot (Z_0 - x).
    """
    if 0:
        print(f'Z_0 = {Z_0}')
        print(f'Lambda = {Lambda}')
        print(f'x = {x}')
        print(f'g_func(x) = {g_func(x)}')
        print(f'g_func(x) + np.dot(Lambda, (Z_0 - x)) = {g_func(x) + np.dot(Lambda, (Z_0 - x))}')

    n = len(Z_0)
  #  x_full = np.zeros(n + 1)
  #  x_full[-1] = 1 - np.sum(x)
  #  x_full[:-1] = x
    return g_func(x, *g_func_args) + np.dot(Lambda, (Z_0 - x))
    #  return g_func(x) + np.dot(Lambda, (Z_0 -  x_full))

##############################################################################
# 2) UBD problem construction
##############################################################################
def ubd(X_D, Z_0, g_func, G_P=0, lambda_bound=1e2, g_func_args=()):
    """
    Builds the linear program (LP) for the upper bounding problem:

      maximize  eta
      subject to:
        eta <= G(x_d) + Lambda^T (Z_0 - x_d)   for each x_d in X_D
        eta <= G(Z_0)
        -lambda_bound <= Lambda_i <= lambda_bound
        -∞ <= eta <= +∞

    We implement this in linprog by minimizing -eta with these constraints in
    a standard "A_ub @ x <= b_ub" format.

    Parameters
    ----------
    X_D : list of 1D array_like
        Already discovered composition points. length(X_D) = m, each x_d has dimension n.
    Z_0 : 1D array_like
        Feed composition, dimension n, sum(Z_0)=1, Z_0[i] >= 0.
    g_func : callable
        G(x) => float, the Gibbs free energy function.
    lambda_bound : float, optional
        Bound on each Lambda_i to avoid unbounded slopes. Default = 1e2.

    Returns
    -------
    c, A_ub, b_ub, bounds
        For scipy.optimize.linprog( ... , method='highs').

    Notes
    -----
    c = [0...0, -1], so we are minimizing -eta => maximizing eta.
    For each x_d:
      +1 * eta - Lambda^T (Z_0 - x_d) <= G(x_d).
    For the global feed:
      +1 * eta <= G(Z_0).
    Then we set bounds: Lambda_i in [-lambda_bound, lambda_bound], eta in [-∞, +∞].
    """

    import numpy as np
    n = len(Z_0)

    # c => minimize -eta => c[-1] = -1
    #c = np.zeros(n)
    c = np.zeros(n + 1)
    c[-1] = -1.0

    num_points = len(X_D)
    #A_ub = np.zeros((num_points + 1, n + 1))
    #A_ub = np.zeros((num_points + 1, n))
    A_ub = np.zeros((num_points + 1, n + 1))
    print(f'A_ub.shape = {A_ub.shape}')
    print(f'n = {n}')
    b_ub = np.zeros(num_points + 1)

    # (1) Constraint: eta <= G(Z_0)
    G_feed = G_P # = g_func(Z_0)
    print(f'G_feed = {G_feed}')
    #A_ub[num_points, -1] = +1.0   # +1 * eta
    A_ub[num_points, -1] = +1.0   # +1 * eta
    b_ub[num_points] = G_feed

    # (2) For each x_d:  eta <= G(x_d) + Lambda^T(Z_0 - x_d).
    # In "A_ub @ vars <= b_ub" form:
    #   +eta - sum_i[Lambda_i * (Z_0[i] - x_d[i])] <= G(x_d).
    # That means for row k:
    #   A_ub[k, i]   = - (Z_0[i] - x_d[i]) = (x_d[i] - Z_0[i]).
    #   A_ub[k, -1]  = +1.
    #   b_ub[k]      = G(x_d).

    for k, x_d in enumerate(X_D):
        #for i in range(n - 1):
        #print(f'list(range(n)) = {list(range(n))}')
        #print(f'x_d = {x_d}')
        #print(f' Z_0 = { Z_0}')
        for i in range(n):
            A_ub[k, i] = x_d[i] - Z_0[i]  # Equal to -(Z_0[i] - x_d[i])

        A_ub[k, -1] = +1.0  # A_ub[:, -1] = set all eta coefficients = 1
        if 0:
            print(f'b_ub[k] = {b_ub[k] } '
                  f'g_func(x_d) = {g_func(x_d, *g_func_args)}')
        b_ub[k] = g_func(x_d, *g_func_args)

    # (3) Bounds for decision vars [Lambda_1..Lambda_n, eta].
    # Limit Lambdas to [-lambda_bound, lambda_bound],  eta is unbounded.
    import math
    big_inf = 1.0e15
    bounds = [(-lambda_bound, lambda_bound)] * n + [(-big_inf, big_inf)]
    #bounds = [(1e-10, lambda_bound)] * n + [(-np.inf, np.inf)]
    #bounds = [(1e-10, lambda_bound)] * (n - 1) + [(-np.inf, np.inf)]
    if 0:
        print(f'c = {c}')
        print(f'A_ub = {A_ub}')
        print(f'b_ub = {b_ub}')
        print(f'bounds = {bounds}')
    return c, A_ub, b_ub, bounds


##############################################################################
# 3) Main Dual Solver
##############################################################################
def solve_dual_equilibrium(
    g_func,
    Z_0,
    g_func_args=(),
    shgo_n=10,
    tol=1e-9,
    max_iter=20,
    min_x=1e-10,
    max_x=1-1e-10,
    lambda_bound=1e2
):
    """
    Solves for a dual hyperplane using the Mitsos/Barton approach, ensuring:
      - x[i] in [min_x, max_x]
      - sum(x)=1
      - Lambda_i in [-lambda_bound, lambda_bound]

    This prevents vertical hyperplanes (x[i]=0 or x[i]=1) and
    unbounded slopes (Lambda -> +/- inf).

    Parameters
    ----------
    g_func : callable
        G(x) => float, x has dimension n, sum(x)=1, x[i]>=0 normally.
    Z_0 : 1D array_like
        Feed composition of length n, sum(Z_0)=1.
    shgo_n : int
        Number of samples for shgo in each LBD subproblem.
    tol : float
        Convergence tolerance for (UBD - LBD).
    max_iter : int
        Maximum dual iterations.
    min_x : float, optional
        Lower bound for each x[i], default=1e-10.
    max_x : float, optional
        Upper bound for each x[i], default=1-1e-10.
    lambda_bound : float, optional
        Bound for each Lagrange multiplier in the LP, default=1e2.

    Returns
    -------
    x_sol : 1D array
        Composition that yields the final LBD.
    Lambda_sol : 1D array
        The final Lagrange multipliers for the dual tangent plane.
    history : dict
        Tracking info: 'iterations', 'UBD', 'LBD', 'Lambda', 'X_star'.

    Notes
    -----
    - By bounding x in [min_x, max_x] we avoid x=0 or x=1.
    - By bounding Lambda in [-lambda_bound, lambda_bound], we avoid
      infinite slopes in the LP solution.
    - If feed is truly unstable, you'll eventually see a nonzero slope.
      Otherwise, if feed is stable or no lower G(x) is found,
      solver might converge with Lambda ~ 0.
    """
    #TODO: Add g_func_args
    Z_0 = np.asarray(Z_0, dtype=float)
    n = len(Z_0)
    assert abs(Z_0.sum() - 1.0) < 1e-12, "Z_0 must sum to 1"
    assert np.all(Z_0 >= 0), "Z_0 must be >= 0 in each component"

    Z_0 = Z_0[:-1]

    # Build initial set X_D with corner-ish points + feed
    X_D = []
    eps = 1e-8

    # "Almost pure" corners: x[i] near 1, distribute small remainder to others
    # NOTE: Should be reduced by to n-1, but we need to keep n for the LP.
    if 0:  # Old code below:
        for i in range(n):
            x_corner = np.full(n, min_x)
            print(f'x_corner = {x_corner}')
            x_corner[i] = 1.0 - eps * (n - 1)
            # etc., but ensure x_corner is length n, i.e. old (n-1).
            # Clip, normalize, etc.
            X_D.append(x_corner)
            # near-zero in i-th component => x[i]=min_x
            # distribute leftover among the other (n-1) in [min_x, max_x]
            if n == 1:
                # Fallback: there's no "other" component
                x_low = np.array([1.0])
            else:
                x_low = np.full(n, (1.0 - min_x * n) / (n - 1))
            x_low[i] = min_x
            # Clip & re-normalize
            x_low = np.clip(x_low, min_x, max_x)
            s_ = x_low.sum()
            if s_ < 1e-14:
                # fallback
                x_low[i] = 1.0
            else:
                x_low /= s_
            X_D.append(x_low)

    #print(f'X_D = {X_D}')
    # New code attempt:
    print(f'list(range(n-1))={list(range(n-1))}')
    if 1:
        lp_bounds = [None] * (n - 1)  # we'll store (L_i, U_i) for each i

        for i in range(n - 1):
            # Generate corner "bar" point:
            x_bar_i = copy.copy(Z_0)
            x_bar_i[i] = 1.0 - np.sum(np.delete(Z_0, i))
            x_bar_i = np.clip(x_bar_i, min_x, max_x)

            # Generate corner "hat" point:
            x_hat_i = copy.copy(Z_0)
            x_hat_i[i] = 0.0
            x_hat_i = np.clip(x_hat_i, min_x, max_x)

            # Evaluate G-values:
            G_feed = g_func(Z_0, *g_func_args)
            G_bar = g_func(x_bar_i, *g_func_args)
            G_hat = g_func(x_hat_i, *g_func_args)

            # Build candidate slopes => each corner yields one slope, e.g.:
            #    slope_bar = [ G_bar - G_feed ] / [ Z_0[i] - x_bar_i[i] ]
            #    slope_hat = ...
            # or whichever formula you prefer.  We'll do a quick example:

            slopes_min = []
            slopes_max = []

            denom_bar = Z_0[i] - x_bar_i[i]
            if abs(denom_bar) > 1e-12:
                slope_bar = (G_bar - G_feed) / denom_bar
                if denom_bar > 0:
                    slopes_max.append(slope_bar)
                else:
                    slopes_min.append(slope_bar)

            denom_hat = Z_0[i] - x_hat_i[i]
            if abs(denom_hat) > 1e-12:
                slope_hat = (G_hat - G_feed) / denom_hat
                if denom_hat > 0:
                    slopes_max.append(slope_hat)
                else:
                    slopes_min.append(slope_hat)

            # Combine candidate min / max:
            lambda_min = max(slopes_min) if slopes_min else -1e6
            lambda_max = min(slopes_max) if slopes_max else +1e6

            lp_bounds[i] = (lambda_min, lambda_max)
    else:
        for i in range(n - 1):
            print(f'i = {i}')
            x_bar_i = copy.copy(Z_0)
            x_bar_i[i] = 1 - np.sum(np.delete(Z_0, i))
            x_bar_i = np.clip(x_bar_i, min_x, max_x)
            print(f'x_bar_i = {x_bar_i}')
            x_bar_i = np.atleast_1d(x_bar_i)
            X_D.append(x_bar_i)

           # x_hat_i = np.zeros(n-1)
            x_hat_i = copy.copy(Z_0)#[i]
            x_hat_i[i] = 0.0
            x_hat_i = np.clip(x_hat_i, min_x, max_x)
            x_hat_i = np.atleast_1d(x_hat_i)
            X_D.append(x_hat_i)

    #X_D.append(np.array([0.0000001]))  # tes
    #X_D.append(np.array([0.9999999]))
    #print(f'X_D = {X_D}')
    # Add feed
    Z0_clipped = np.clip(Z_0, min_x, max_x)  # Removed because this should never be needed
    #Z_0 = np.clip(Z_0, min_x, max_x)
    print(f'Z0_clipped = {Z0_clipped}')
    print(f'Z_0 = {Z_0}')
    #print(f'Z0_clipped = {Z0_clipped}')
    sum_feed = Z_0.sum()
    if sum_feed < 1e-14:
        Z_0[0] = 1.0  # fallback
    else:
        Z_0 /= sum_feed

    #X_D.append(np.array([0.0]))

    X_D.append(Z0_clipped)

    print(f'X_D = {X_D}')
    # Evaluate initial LBD & UBD
    LBD = -1e15
    UBD = g_func(Z0_clipped, *g_func_args)

    history = {
        'iterations': [],
        'UBD': [],
        'LBD': [],
        'Lambda': [],
        'X_star': []
    }

    G_P = g_func(Z_0, *g_func_args)  # G_feed remains unchanged throughout the iterations
    for iteration in range(1, max_iter+1):

        # (A) Solve UBD via LP
        c, A_ub, b_ub, lp_bounds = ubd(X_D, Z0_clipped, g_func, lambda_bound=lambda_bound, G_P=G_P,
                                       g_func_args=g_func_args)
        lp_res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=lp_bounds, method='highs')
        if not lp_res.success:
            print(f"[Iter {iteration}] LP failed: {lp_res.message}")
            break
        print(f'lp_res.x = {lp_res.x}')
        #Lambda_sol = lp_res.x[:n]  # the multipliers
        Lambda_sol = lp_res.x[:-1]  # the multipliers
        print(f'Lambda_sol = {Lambda_sol}')
        eta_sol = lp_res.x[-1]
        cur_UBD = -lp_res.fun

        # (B) Solve LBD with shgo
        if 0:  # no longer needed
            def lbd_obj(xprime, gf, lam, feed):
                x = np.append(xprime, 1 - np.sum(xprime))
                return lbd(x, gf, lam, feed, g_func_args=g_func_args)

        # Bounds: x in [min_x, max_x]
        print(f'n = {n}')
        print(f'min_x = {min_x}')
        print(f'max_x = {max_x}')
        shgo_bounds = [(min_x, max_x)] * (n - 1)

        # Prepare options dict:
        #TODO: Add to solver after fixing
        options = {}
        symm = np.zeros(n)
        options['symmetry'] = symm

        # sum(x)=1 => eq constraint
        def sum_constraint(x):
            return np.sum(x) # <= 1.0

        nonlin_con = NonlinearConstraint(sum_constraint, -np.inf, 1.0)

        # Run main solver
        #print(f'lbd_obj = {;}')
        print(f'shgo_bounds = {shgo_bounds}')
        res = shgo(#lbd_obj,
                   func=lbd,
                   bounds=shgo_bounds,
                   args=(g_func, Lambda_sol, Z0_clipped, g_func_args),
                   constraints=nonlin_con,
                   n=shgo_n
        )
        print(f'res = {res}')
        #TODO: The speed of the method can be improved by adding all solutions of shgo
        #     (with res.funl values near zero) to the X_D set instead of just x_star.
        if not res.success:
            print(f"[Iter {iteration}] SHGO LBD fail: {res.message}")
            break

        #x_star = np.append(res.x, 1 - np.sum(res.x))
        x_star = res.x
        val_star = res.fun
        cur_LBD = val_star

        # Store iteration results
        history['iterations'].append(iteration)
        history['UBD'].append(cur_UBD)
        history['LBD'].append(cur_LBD)
        history['Lambda'].append(Lambda_sol)
        history['X_star'].append(x_star)

        # Add discovered point to X_D
        X_D.append(x_star)

        # Termination check
        gap = abs(cur_UBD - cur_LBD)
        print(f"[Iter {iteration}] UBD={cur_UBD:.6g}, LBD={cur_LBD:.6g}, gap={gap:.3g}")
        if gap < tol:
            print(f"Converged at iteration {iteration} with gap={gap:.3g}")
            break

        # Prepare next iteration
        LBD = cur_LBD
        UBD = cur_UBD

    return x_star, Lambda_sol, history




# %% Dual plane eta for optim visualization
def dual_plane_sol(x, G_sol, Lambda_sol, x_star):
    """
    Computes the value of the tangent plane (dual solution hyperplane) at a
    composition x, given a reference composition x_star, its Gibbs free
    energy G_sol, and multipliers Lambda_sol.

    Parameters
    ----------
    x : array_like of length n
        Composition at which the hyperplane is to be evaluated.
    G_sol : float
        Gibbs free energy at x_star, i.e. G(x_star).
    Lambda_sol : array_like of length n
        Final Lagrange multipliers for the dual tangent plane.
        Typically found via solve_dual_equilibrium.
    x_star : array_like of length n
        The "reference" composition used to define the plane. Often the
        composition from the final LBD subproblem in the dual approach.

    Returns
    -------
    float
        The plane value T(x) = G_sol + sum_i [Lambda_sol[i] * (x[i] - x_star[i])].

    Notes
    -----
    - This function is dimension-independent: n can be 2 for binary, 3 for ternary, etc.
    - `G_sol` must be consistent with `x_star` and `g_func(x_star)`.
    """
    x = np.asarray(x, dtype=float)
    x_star = np.asarray(x_star, dtype=float)
    Lambda_sol = np.asarray(Lambda_sol, dtype=float)

    return float(G_sol + np.dot(Lambda_sol, (x - x_star)))


def dual_lagrange(x, g_func, Lambda_sol, Z_0):
    """
    Evaluates the "Lagrangian" expression at composition x for the dual plane,
    i.e. G(x) + sum(Lambda_sol * (Z_0 - x)).

    Parameters
    ----------
    x : array_like of length n
        Composition at which to evaluate. Must satisfy sum(x) = 1, x[i] >= 0
        if in a physical simplex.
    g_func : callable
        A Gibbs free energy function G(x) => float, where x is length n.
    Lambda_sol : array_like of length n
        The final Lagrange multipliers for the dual tangent plane.
    Z_0 : array_like of length n
        The feed composition (must satisfy sum(Z_0)=1, Z_0[i]>=0).

    Returns
    -------
    float
        Value = G(x) + sum_i [Lambda_sol[i] * (Z_0[i] - x[i])].

    Notes
    -----
    - This is typically used to visualize or verify the dual hyperplane
      at various x.
    - In the Mitsos/Barton dual approach, the feed composition Z_0 is the
      "primal" reference point, though you are free to interpret it otherwise.
    """
    x = np.asarray(x, dtype=float)
    Z_0 = np.asarray(Z_0, dtype=float)
    Lambda_sol = np.asarray(Lambda_sol, dtype=float)

    return float(g_func(x) + np.dot(Lambda_sol, (Z_0 - x)))
