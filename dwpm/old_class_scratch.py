# %% Define Gibbs energy functions for VdW EoS (cubic with 2 volume roots)
def g_R_k_i(s, p, k='x', i=1):  # (Validated)
    """
    Residual Gibbs energy g_R_k_i((T,P) of a single component where k is the
    specified phase and i is the component in phase k.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of the specified component.

    p : class
        Contains the dictionary describing the component parameters.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'.

    i : integer, optional
        The component number in phase k to calculate. ex. 1

    Dependencies
    ------------
    math

    Returns
    -------
    g_R_k_i : scalar output.
    """
    from numpy import log
    if k == 'y':  # 'y' = Vapour phase standard
        V = s.c[i]['V_v']
    else:  # Assume all phases other than vapour are liquid, ex. 'x'
        V = s.c[i]['V_l']

    return (s.c[i]['P'] * V / (p.c[i]['R'] * s.c[i]['T']) - 1.0
            - log(s.c[i]['P'] / (p.c[i]['R'] * s.c[i]['T']))
            - log(V - s.c[i]['b'])
            - s.c[i]['a'] / (p.c[i]['R'] * s.c[i]['T'] * V))


def g_R_mix_i(s, p, k='x'):  # (Validated)
    """
    Residual Gibbs energy g_R_mix_i((T,P) of a mixture where k is the specified
    phase.


    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of the specified component.

    p : class
        Contains the dictionary describing the component parameters.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'.

    Dependencies
    ------------
    math

    Returns
    -------
    g_R_k_i : scalar output.
    """
    from numpy import log
    if k == 'y':  # 'y' = Vapour phase standard
        V = s.m['V_v']
    else:  # Assume all phases other than vapour are liquid, ex. 'x'
        V = s.m['V_l']

    if V < s.m['b']:
        import logging
        # V = (1.0 + 100000*(s.m['b'] - V)) * s.m['b']
        if V > 0:
            V = s.m['b'] + 1e-15 / (s.m['b'] - V)  # * (1.0 + 1e-30)
        if V < 0:
            V = s.m['b'] + 1e-15 / abs(V)

        # logging.warning("V < b_mix in g_R_mix_i, setting to V = {}".format(V))
        # logging.warn("V < b_mix in g_R_mix_i, setting to V = {}".format(V))

    return (s.m['P'] * V / (p.m['R'] * s.m['T']) - 1.0
            - log(s.m['P'] / (p.m['R'] * s.m['T']))
            - log(V - s.m['b'])
            - s.m['a'] / (p.m['R'] * s.m['T'] * V))


def g_IG_k(s, p, k='x'):  # (Validated)
    """
    Change in gibbs energy for mixing of ideal gasses.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'.

    Dependencies
    ------------
    math

    Returns
    -------
    g_IG_k : scalar output.
    """
    from numpy import log

    Sigma_g_IG_k = 0.0  # Sum of ideal gas terms
    for i in range(1, p.m['n'] + 1):
        if s.c[i][k] == 0.0:  # Prevent math errors from zero log call.
            pass  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0
        else:
            if s.c[i][k] < 0.0:
                # TODO: This should never step outside bounds, found out why
                # s.c[2][k] is ofter < 0
                print('s.c[{}][{}] = {}'.format(i, k, s.c[i][k]))
                # s.c[i][k] = abs(s.c[i][k])
            Sigma_g_IG_k += s.c[i][k] * log(s.c[i][k])
    return Sigma_g_IG_k


def g_mix(s, p, k=None, ref='x', update_system=False):  # (Validated)
    """
    Returns the gibbs energy at specified composition relative to a reference
    phase pure component gibbs energy.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    k : string, optional # TODO UNFINISHED
        Force phase to be calculated, ex. liquid phase 'x'.

    ref : string, optional
          Selected reference phase. Note that is common practice to choose a
          liquid phase 'x' as the reference phase.

    update_system : boolean, optional # TODO UNFINISHED
                    This updates the system state to the current P, T and
                    composition conditions. Only use if the system dictionary
                    has not been updated to the current independent variables.

    Dependencies
    ------------
    math

    Returns
    -------
    s : class output.
        Contains the following values (or more if more phases are chosen):
          s.m['g_mix']['t'] : scalar, Total Gibbs energy of mixing.
          s.m['g_mix']['x'] : scalar, Gibbs energy of mixing for liquid phase.
          s.m['g_mix']['y'] : scalar, Gibbs energy of mixing for vapour phase.
          s.m['g_mix']['ph min'] : string, contains the phase/volume root of
                                           with lowest Gibbs energy.
          s.s['Math Error'] : boolean, if True a math error occured during
                                       calculations. All other values set to 0.
    """
    import logging
    if update_system:  # TODO TEST
        Xvec = [[]]  # Construct update vector
        for i in range(1, p.m['n']):  # for n-1 independent components
            Xvec[0].append(s.c[i][k])

        s.update_state(s, p, P=s.m['P'], T=s.m['T'], phase=k,
                       X=Xvec)

    s.update_state(s, p)  # Update Volumes and activity coeff.

    # try:
    Sigma_g_ref = 0.0
    for i in range(1, p.m['n'] + 1):
        Sigma_g_ref -= s.c[i][ref] * g_R_k_i(s, p, k=ref, i=i)

    s.m['g_mix^R'] = {}
    for ph in p.m['Valid phases']:
        s.m['g_mix^R'][ph] = g_R_mix_i(s, p, k=ph) + Sigma_g_ref

    s.m['g_mix'] = {}
    g_min = []
    g_abs_min = numpy.inf  # long  # "inf" large int
    for ph in p.m['Valid phases']:
        s.m['g_mix'][ph] = s.m['g_mix^R'][ph] + g_IG_k(s, p, k=ph)
        g_min.append(s.m['g_mix'][ph])
        if s.m['g_mix'][ph] < g_abs_min:  # Find lowest phase string
            s.m['g_mix']['ph min'] = ph
            g_abs_min = s.m['g_mix'][ph]

    s.m['g_mix']['t'] = min(g_min)

    s.s['Math Error'] = False

    # except(ValueError, ZeroDivisionError):
    #     import numpy
    #     s.m['g_mix'] = {}
    #     s.s['Math Error'] = True
    #     logging.error('Math Domain error in g_mix(s,p)')
    #     for ph in p.m['Valid phases']:
    #             s.m['g_mix'][ph] = numpy.nan#0.0
    #
    #     s.m['g_mix']['t'] = numpy.nan#0.0

    return s


# %% Duality formulation
def ubd(X_D, Z_0, g_x_func, s, p, k=None):
    """
    Returns the arguments to be used in the optimisation of the upper bounding
    problem with scipy.optimize.linprog.

    used
    Parameters
    ----------
    X_D : vector (1xn array)
          Contains the current composition point in the overall dual
          optimisation. Constant for the upper bounding problem.

    Z_0 : vector (1xn array)
          Feed composition. Constant.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    k : list, optional (TODO)
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.

    Returns
    -------
    c : array_like
        Coefficients of the linear objective function to be minimized.

    A : A_eq : array_like, optional
        2-D array which, when matrix-multiplied by x, gives the values of the
        upper-bound inequality constraints at x.

    b : array_like, optional
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in.
    """
    import numpy
    # Coefficients of UBD linear objective function
    c = numpy.zeros([p.m['n']])  # linrpog maximize/minimize? D
    # Documentation is contradictory across version; check
    c[p.m['n'] - 1] = -1.0  # -1.0 change max --> min problem

    # Coefficients of Lambda inequality constraints

    A = numpy.zeros([len(X_D) + 1,  # rows = for all X_D + Global ineq
                     p.m['n']]  # cols = n comps + eta
                    )
    b = numpy.zeros(len(X_D) + 1)

    # Global problem bound (Fill last row of A and last element in b
    # G_p (Primal problem Z_0_i - x_i = 0 for all i)
    # TODO: Move outside function and outside loop in dual
    s = s.update_state(s, p, X=Z_0, Force_Update=True)
    G_P = g_x_func(s, p).m['g_mix']['t']
    A[len(X_D), p.m['n'] - 1] = 1  # set eta to 1
    b[len(X_D)] = G_P

    # Bounds for all X_d in X_D
    A[:, p.m['n'] - 1] = 1  # set all eta coefficients = 1
    for X, k in zip(X_D, range(len(X_D))):
        # Find G(X_d)
        s = s.update_state(s, p, X=X, Force_Update=True)
        # TODO: This only needs to be evaluated once for every x \in X^D
        G_d = g_x_func(s, p).m['g_mix']['t']
        b[k] = G_d
        for i in range(p.m['n'] - 1):
            A[k, i] = -(Z_0[i] - X_D[k][i])

    if False:
        print('c shape = {}'.format(numpy.shape(c)))
        print('A shape = {}'.format(numpy.shape(A)))
        print('b shape = {}'.format(numpy.shape(b)))
        print('c = {}'.format(c))
        print('A = {}'.format(A))
        print('b = {}'.format(b))

    return c, A, b


def lbd(X, g_x_func, Lambda_d, Z_0, s, p, k=['All']):
    """
    Returns the lower bounding problem of the dual extremum.

    Parameters
    ----------
    X : vector (1xn array)
        Contains the current composition point in the overall dual
        optimisation to be optimised to the minimum value of the lbd.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    Lambda_d : vector (1xn array)
               Contains the diality multipliers Lambda \in R^m.
               Constant for the lower bounding problem.

    Z_0 : vector (1xn array)
          Feed composition. Constant.

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.

    Dependencies
    ------------
    numpy.array
    math.e

    Returns
    -------
    lbd : scalar
          Value of the lower bounding problem at X.
    """
    # Update system to new composition.
    s.update_state(s, p, X=X, phase=k, Force_Update=True)

    return g_x_func(s, p).m['g_mix']['t'] + sum(Lambda_d * (Z_0 - X))


def dual_equal(s, p, g_x_func, Z_0, k=None, P=None, T=None, tol=1e-9, n=100):
    """
    Dev notes and TODO list
    -----------------------
    TODO: -The material bounds is too high since (Mitsos')  \hat{X} is the true
            upper limit for a given feedpoint
          -Look into using X_bounds scheme of Pereira instead for low mole Z_0
          -Add valid phases option.
          -Add constraints to x_i =/= 0 or 1 for all i to prevent vertical
            hyperplanes.
          -Construct bounds which is a list of tuples in [0,1] \in R^n

    NOTES: -Strictly the composition in all phases in should be specified in
            X_d, refer to older versions of this script when different comp.
            spaces need to be used.
    -----------------------
    Find the phase equilibrium solution using the daul optimization algorithm.
    Ref. Mitsos and Barton (2007)

    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.

    P : scalar, optional
        Pressure (Pa), if unspecified the current state pressure will be used.

    T : scalar, optional
        Temperature (K), if unspecified  the current state temperature will be
        used.

    Z_0 : vector
          Contains the feed composition point (must be and unstable point to
          find multiphase equilibria).

    tol : scalar, optional
          Tolerance, if epsilon >= UBD - LBD that will terminate the routine.

    n : scalar, optional
        Number of sampling points used in the tgo routine in solving LBD of the
        dual problem.
        Note: It is recommended to use at least ``100 + p.m['n'] * 100``

    Dependencies
    ------------
    numpy

    Returns
    -------
    X_sol : vector
            Contains the first optimised equilibrium point of the dual problem

    Lambda_sol : vector
                 Contains the optimised lagrange multipliers (partial chemical
                 potential) sol. of the dual solution hyperplane

    d_res : optimisation object return
            Contains the final solution to the dual problem with the
            following values:
                d_res.fun : lbd plane solution at equil point
                d_res.xl : Other local composition solutions from final tgo
                d_res.funl : lbd plane at local composition solutions

    """
    import numpy
    from scipy.optimize import linprog
    from tgo import tgo

    def x_lim(X):  # limiting function used in TGO defining material constraints
        import numpy
        # return -numpy.sum(X, axis=-1) + 1.0
        return -numpy.sum(X, axis=0) + 1.0

    if k is None:
        k = p.m['Valid phases']

    # Initialize
    Z_0 = numpy.array(Z_0)
    LBD = - 1e300  # -inf
    s.update_state(s, p, X=Z_0, phase=k, Force_Update=True)

    # G_p (Primal problem Z_0_i - x_i = 0 for all i):
    UBD = g_x_func(s, p).m['g_mix']['t']

    # X bounds used in UBD optimization
    X_bounds = [[],  # Upper bound (bar x)
                []  # Lower bound (hat x)
                ]

    for i in range(p.m['n'] - 1):
        # Append an independent coordinate point for each bound
        X_bounds[0].append(numpy.zeros(shape=(p.m['n'] - 1)))
        X_bounds[1].append(numpy.zeros(shape=(p.m['n'] - 1)))

        # Set upper bound coordinate point i
        Sigma_ind = 0.0  # Sum of independent components excluding i
        # (lever rule)
        for k_ind in range(p.m['n'] - 1):  # Note: k var name is used as phase
            # Set k != i  (k==i changed at end of loop)
            X_bounds[0][i][k_ind] = Z_0[k_ind]
            if k_ind != i:
                Sigma_ind += Z_0[k_ind]  # Test, use numpy.sum if working

                # Set Lower bound coordinate point i
                X_bounds[1][i][k_ind] = Z_0[k_ind]
                # (Remaining bound coordinate points kept at zero)

        X_bounds[0][i][i] = 1.0 - Sigma_ind  # change from Z_0[k]

    # Construct physical bounds x \in [0, 1] for all independent components
    Bounds = []
    L_bounds = []  # Lambda inf bounds used in linprog.
    for i in range(p.m['n'] - 1):
        Bounds.append((1e-10, 0.99999999))
        L_bounds.append((-numpy.inf, numpy.inf))

    L_bounds.append((-numpy.inf, numpy.inf))  # Append extra bound set for eta

    # Update state to random X to initialise state class.
    X_sol = numpy.array(X_bounds[1][0])  # set X_sol to lower bounds
    s.update_state(s, p, X=X_sol, phase=k, Force_Update=True)

    X_D = []  # set empty list
    # Add every X_bounds to X_D list to use in linprog
    for i in range(p.m['n'] - 1):
        X_D.append(X_bounds[0][i])
        X_D.append(X_bounds[1][i])

    if True:  # Lambda estimates using differentials at Z_0
        # NOTE on CO2-ethane test this cut the iterations down to 6 from 9
        Lambda_d = numpy.zeros_like(Z_0)
        s.update_state(s, p, X=Z_0, phase=k, Force_Update=True)
        for z in range(1, p.m['n']):
            Lambda_d[z - 1] = FD(g_mix, s, p, d=1, z=z, gmix=True)
            # print('Lambda_d from init FD est. = {}'.format(Lambda_d))

        # Solve LBD for first cutting plane
        d_res = tgo(lbd, Bounds, args=(g_x_func, Lambda_d, Z_0, s, p, k),
                    g_cons=x_lim,
                    n=n,
                    # k_t=2,
                    # n = 100 + 100*(p.m['n'] - 1),
                    )  # skip=2)

        X_sol = d_res.x
        X_D.append(X_sol)
        # print('X_sol from init FD est. = {}'.format(X_sol))
        if len(d_res.xl) > 0:
            for i in range(len(d_res.xl)):
                # print('d_res.xl{}'
                #      ' from init FD est. = {}'.format(i, d_res.xl[i]))
                X_D.append(d_res.xl[i])

    # print('X_D at init = {}'.format(X_D))

    # %% Normal calculation of daul problem if Z_0 is unstable.
    iteration = 0
    # X_D.append(numpy.array([ 0.19390632]))
    while abs(UBD - LBD) >= tol:
        iteration += 1
        # Solve UBD
        # Find new bounds for linprog
        c, A, b = ubd(X_D, Z_0, g_x_func, s, p)
        # Find mulitpliers with max problem.
        lp_sol = linprog(c, A_ub=A, b_ub=b, bounds=L_bounds)
        Lambda_sol = numpy.delete(lp_sol.x, numpy.shape(lp_sol.x)[0] - 1)
        # If float convert back to 1x1 array
        Lambda_sol = numpy.array(Lambda_sol)

        UBD = -lp_sol.fun  # Final func value is neg. of minimised max. problem

        # dual stepping plots
        if 0:
            print('Iteration number: {}'.format(iteration))
            # print('Lambda_sol: {}'.format(Lambda_sol))
            print('X_sol: {}'.format(X_sol))
            print('X_D: {}'.format(X_D))
            x_r = 1000
            # Lagrange function surface
            plane_args = (Lambda_sol, Z_0, g_x_func, s, p, ['All'])
            plot.plot_ep(dual_lagrange, x_r, s, p, args=plane_args)

            # Dual plane
            if p.m['n'] == 2:
                s.update_state(s, p, X=X_sol[0], Force_Update=True)
                G_sol = g_x_func(s, p).m['g_mix']['t']
                print('G_sol : {}'.format(G_sol))
                plot.plot_g_mix(s, p, g_x_func, Tie=[[Z_0, X_sol]], x_r=1000,
                                plane_func=dual_plane_sol,
                                # plan_args=(G_sol, -Lambda_sol, Z_0, X_sol))
                                plan_args=(G_sol, Lambda_sol, Z_0, X_sol))

        # Solve LBD
        d_res = tgo(lbd, Bounds, args=(g_x_func, Lambda_sol, Z_0, s, p, k),
                    g_cons=x_lim,
                    n=n,
                    # n = 100 + 100*(p.m['n'] - 1),
                    )  # skip=2)

        X_sol = d_res.x
        X_D.append(X_sol)

        # if True:  # NOTE: Reduced iterations from 6 to 3 !
        if False:
            if len(d_res.xl) > 0:
                for i in range(len(d_res.xl)):
                    # print('X_D = {}'.format(X_D))
                    # print('d_res.xl) = {}'.format(d_res.xl))
                    X_D.append(d_res.xl[i])

        # Calculate LBD
        LBD = lbd(X_sol, g_x_func, Lambda_sol, Z_0, s, p, k)
        # End if tol

        if True:  # dual stepping plots
            print('Iteration number: {}'.format(iteration))
            # print('Lambda_sol: {}'.format(Lambda_sol))
            # print('X_sol: {}'.format(X_sol))
            print('X_D: {}'.format(X_D))
            x_r = 1000
            # Lagrange function surface
            # plane_args = (Lambda_sol, Z_0, g_x_func, s, p, ['All'])
            # plot.plot_ep(dual_plane, x_r, s, p, args=plane_args)

            # Dual plane
            if p.m['n'] == 2:
                s.update_state(s, p, X=X_sol[0], Force_Update=True)
                G_sol = g_x_func(s, p).m['g_mix']['t']
                print('G_sol: {}'.format(G_sol))
                plot.plot_g_mix(s, p, g_x_func, Tie=[[Z_0, X_sol]], x_r=1000,
                                plane_func=dual_plane_sol,
                                # plan_args=(G_sol, -Lambda_sol, Z_0, X_sol))
                                plan_args=(G_sol, Lambda_sol, Z_0, X_sol))

    if False:  # Print results optional
        print('Final UBD = {}'.format(UBD))
        print('Final LBD = {}'.format(LBD))
        print('Final UBD - LBD = {}'.format(UBD - LBD))
        print('Final Z_eq = {}'.format(X_sol))
        print('Final Lambda_d = {}'.format(Lambda_d))

    if False:  # Feed point plane estimate dev
        x_r = 1000
        # Suppose data_solutions at
        # [[array([ 0.1939063]), array([ 0.30898849])]]
        X_I = numpy.array([0.1939063])
        s.update_state(s, p, X=X_I, Force_Update=True)
        G_sol_I = g_x_func(s, p).m['g_mix']['t']
        X_II = numpy.array([0.30898849])
        s.update_state(s, p, X=X_II, Force_Update=True)
        G_sol_II = g_x_func(s, p).m['g_mix']['t']
        print(' G_sol_I = {}'.format(G_sol_I))
        print(' G_sol_II = {}'.format(G_sol_II))
        # Plane estimates
        # (NOTE: These lambda estimates need to be done for each component
        # in higher dimensions)
        Lambda_sol_est = (G_sol_II - G_sol_I) / (X_II - X_I)
        Z_0 = (X_I + X_II) / 2.0  # (Estimate of feed point)
        plot.plot_g_mix(s, p, g_x_func, Tie=[[Z_0, X_sol]], x_r=1000,
                        plane_func=dual_plane_sol,
                        plan_args=(G_sol_I, Lambda_sol_est, X_II, X_I)
                        # plan_args=(G_sol_I, Lambda_sol_est, Z_0, X_I)
                        )
    # Returns
    return X_sol, Lambda_sol, d_res

