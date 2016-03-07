#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" execfile('tgo.py')
"""
#from __future__ import division, print_function, absolute_import
from UQToolbox.sobol_lib import i4_sobol_generate
 # TODO: Replace with latinhypercube sampling used in differentialevolution.py
import numpy
import scipy.spatial
import scipy.optimize



# %% Define tgo class
class TGO(object):
    """
    This class implements the tgo routine
    """

    def __init__(self, func, bounds, args=(), g_func=None, g_args=(), n=100,
                 skip=1, k_t=None, callback=None, minimizer_kwargs=None,
                 disp=False):

        self.func = func
        self.bounds = bounds
        self.args = args
        self.g_func = g_func
        self.g_args = g_args
        self.n = n
        self.skip = skip
        self.k_t = k_t
        self.callback = callback
        self.minimizer_kwargs = minimizer_kwargs
        self.disp = disp

        # Initialize return object
        self.res = scipy.optimize.OptimizeResult()
        self.res.nfev = n  # Include each sampling point as func evaluation
        self.res.nlfev = 0  # Local function evals for all minimisers
        self.res.nljev = 0  # Local jacobian evals for all minimisers
        #self.res.func

    # %% Define funcs # TODO: Create tgo class to wrap all funcs
   def sampling(self):
        """
        Generates uniform sampling points in a hypercube and scales the points
        to the bound limits.
        """
        # Generate sampling points.
        #  TODO Assert if func output matches dims. found from bounds
        self.m = len(self.bounds)  # Dimensions

        # Generate uniform sample points in R^m
        self.B = i4_sobol_generate(self.m, self.n, self.skip)
        self.C = numpy.column_stack([self.B[i] for i in range(self.m)])

        # Distribute over bounds
        # TODO: Find a better way to do this
        for i in range(len(self.bounds)):
            self.C[:,i] = (self.C[:,i] *
                           (self.bounds[i][1] - self.bounds[i][0])
                           + self.bounds[i][0] )

        return self.C


    def t_matrix(H, F):
        """
        Returns the topographical matrix with True boolean values indicating
        positive entries and False ref. values indicating negative values.
        """
        return (H.T > F.T).T


    def k_t_matrix(T, k):  # TODO: Replace delete with simpler array access
        """Returns the k-t topograph matrix"""
        return numpy.delete(T, numpy.s_[k:numpy.shape(T)[1]], axis=-1)


    def minimizers(K):
        """Returns the minimizer indexes of a k-t matrix"""
        Minimizers = numpy.all(K, axis=-1)
        # Find data point indexes of minimizers:
        return numpy.where(Minimizers)[0]


    def K_optimal(T): # TODO: Recheck correct implementation, compare with HS19
        """
        Returns the optimimal k-t topograph with the semi-empircal correlation
        proposed by Henderson et. al. (2015)
        """
        K_1 = k_t_matrix(T, 1)  # 1-t topograph
        k_1 = len(minimizers(K_1))
        k_i = k_1
        i = 2
        while k_1 == k_i:
            K_i = k_t_matrix(T, i)
            k_i = len(minimizers(K_i))
            i += 1

        ep = i * k_i / (k_1 - k_i)
        k_c = numpy.floor((-(ep - 1) + numpy.sqrt((ep - 1.0)**2 + 80.0 * ep))
                          / 2.0)

        k_opt = int(k_c + 1)
        if k_opt > numpy.shape(T)[1]:  # If size of k_opt exceeds t-graph size.
            k_opt = int(numpy.shape(T)[1])

        K_opt = k_t_matrix(T, k_opt)
        return K_opt



def tgo(func, bounds, args=(), g_func=None, g_args=(), n=100, skip=1, k_t=None,
        callback=None, minimizer_kwargs=None, disp=False):
    """
    Finds the global minima of a function using topograhphical global
    optimisation.

    Parameters TODO: REVIEW DOC
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.

    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.

    args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        objective function.

    g_func : callable, optional
        Function used to define a limited subset to defining the feasible set
        of solutions in R^n in the form g(x) <= 0 applied as g : R^n -> R^m
        ex 1. To impose the constraint  x[0] + x[1] + x[2] - 1 <= 0
              Use the function definition

              def g_func(A):  # A is an input array of sample points in R^m*n
                  return numpy.sum(A, axis=1) - 1.0#

        ex 2. To impose the constraints:
              -(x[0] - 5)**2 - (x[1] - 5)**2  - 100 <= 0
              (x[0] - 6)**2 - (x[1] - 5)**2  - 82.81 <= 0
              Use the function definition

              def g_func(A):
                  return ((-(C[:,0] - 5)**2 - (C[:,1] - 5)**2  - 100.0 <= 0.0)
                         & ((C[:,0] - 6)**2 - (C[:,1] - 5)**2  - 82.81 <= 0.0))


    g_args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        feasible set function.

    n : int, optional
        Number of sampling points used in the construction of the topography
        matrix.

    skip : int, optional
        Tthe number of initial points to skip in the Sobol generation sequence.

    k : int, optional
        Defines the number of columns constructed in the k-t matrix. The higher
        k is the lower the amount of minimizers will be used for local search
        routines. If None the empircal model of Henderson et. al. (2015) will
        be used. (Note: Lower k values decrease performance, but could
        potentially be more robust due to testing more local minimizers in the
        function hypersuface)

# TODO:    disp : bool, optional
        Display status messages

# TODO:    callback : callable, `callback(xk, convergence=val)`, optional:
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns `True`, then the minimization
        is halted (any polishing is still carried out).

# TODO:    minimizer_kwargs : dict, optional
        Extra keyword arguments to be passed to the minimizer
        ``scipy.optimize.minimize()`` Some important options could be:

            method : str
                The minimization method (e.g. ``"L-BFGS-B"``)
            args : tuple
                Extra arguments passed to the objective function (``func``) and
                its derivatives (Jacobian, Hessian).

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If `polish`
        was employed, then OptimizeResult also contains the `jac` attribute.

    Notes
    -----
    Global optimization using the Topographical Global Optimization (TGO)
    method first proposed by Törn (1990) [1] with
    [ TODO: EXPLAIN k-t matrix use]
    Hendorson et. al. (2015) [2].

    The TGO is a clustering method that uses graph theory to generate good
    starting points for local search methods from points distributed uniformly
    in the interior of the feasible set. These points are generated using the
    Sobol (1967) [3] sequence.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] Törn, A (1990) "Topographical global optimization", Reports on
           Computer Science and Mathematics Ser. A, No 199, 8p. Abo Akademi
           University, Sweden
    .. [2] Henderson, N, de Sá Rêgo, M, Sacco, WF, Rodrigues, RA Jr. (2015) "A
           new look at the topographical global optimization method and its
           application to the phase stability analysis of mixtures",
           Chemical Engineering Science, 127, 151-174
    .. [3] Sobol, IM (1967) "The distribution of points in a cube and the
           approximate evaluation of integrals. USSR Comput. Math. Math. Phys.
           7, 86-112.

    """
    # Initiate TGO class
    TGOc = TGO(func, bounds, args=(), g_func=None, g_args=(), n=100, skip=1,
               k_t=None, callback=None, minimizer_kwargs=None, disp=False)

    # Generate sampling points
    TGOc.sampling()

    # Find subspace of feasible points
    if g_func is not None: # TODO: Improve
        C =  C[g_func(C)]  # Subspace of feasible points.

    Y = scipy.spatial.distance.cdist(C, C, 'euclidean')
    Z = numpy.argsort(Y, axis=-1)
    A = numpy.delete(Z, 0, axis=-1)  # Topographical matrix without signs

    # %% Obj. function returns to be used as reference table.:
    F = numpy.zeros(numpy.shape(C)[0])
    for i in range(numpy.shape(C)[0]):
        F[i] = func(C[i,:], *args)
    # TODO: see scipy.spatial.KDTree for F lookup?

    # %% Create float value and bool topograph:
    H = F[A] # This replaces all index values in A with the function result
    #
    T = t_matrix(H, F)  # Topograph with Boolean entries
    # %% Find the optimial k+ topograph
    # Find epsilon_i parameter for current system
    if k_t is None:
        K_opt = K_optimal(T)
    else:
        K_opt = k_t
    # %% Local Search: Find the minimzer float values and

#TODO: IMPROVE

    Min_ind = minimizers(K_opt)
    x_vals = []
    Func_min = numpy.zeros_like(Min_ind)
    for i, ind in zip(range(len(Min_ind)), Min_ind):
        # Find minimum x vals
        x_min = scipy.optimize.minimize(func, C[ind,:],
                                        method='L-BFGS-B',
                                        bounds=bounds,
                                        args=args)['x']
        x_vals.append(x_min)
        # Find func float vals
        Func_min[i] = func(x_min, *args)

    # Find global of all minimizers
    i_glob = numpy.argsort(Func_min)[0]
    x_global_min = x_vals[i_glob]


    # Initialize return object
   # self.res = scipy.optimize.OptimizeResult()
  #  self.res.nfev = n # Include each sampling point as func evaluation

    return x_global_min

    
if __name__ == '__main__':
    pass
    
    
    
 