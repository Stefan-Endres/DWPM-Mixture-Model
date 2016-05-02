#!/usr/bin/env python

class TopShiftParam:
    def __init__(self, p, rs=True, kij=False, rskij=False):
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import dual_equal
        if p.m['Model'] == 'DWPM':
            self.param_func = self.vdw_dwpm_params
            self.rs = rs
            self.kij = kij
            self.rskij = rskij

    def optimise(self, s, p, g_x_func, Z_0,
                 method_d='L-BFGS-B',
                 method_eq='L-BFGS-B',
                 bounds=None):
        """
        optimise the objective function using the dual lagrange shifting method
        to provide a starting value for the total equilibrium optimisation over
        all data points
        """
        import scipy
        import logging
        from tgo import tgo

        # Make logging.info verbose in CLI
        logging.basicConfig(level=logging.DEBUG)

        # Find local minima of dual func
        tsp_args = (s, p, g_x_func, False)

        if method_d == 'tgo':
            res_d = tgo(self.tsp_objective_function, bounds=bounds,
                        args=tsp_args, n=3000)
        else:
            res_d = scipy.optimize.minimize(self.tsp_objective_function,
                                            Z_0,
                                            args=tsp_args,
                                            method=method_d
                                            )

        # Verbose report before computationally expensive equilibrium
        # optimization over all data points
        print('='*100)
        logging.info(res_d)
        if self.rs:
            print("-r {} -s {}".format(res_d.x[0], res_d.x[1]))
        if self.kij:
            print("kij {} {}".format(res_d.x[0], res_d.x[1]))

        if self.rskij:
            print("-r {} -s {} -kij {} {}".format(res_d.x[0],
                                                  res_d.x[1],
                                                  res_d.x[2],
                                                  res_d.x[3]))
        # Optimize over all equilibrium points using dual plane minima as a
        # starting point
        Z_0 = res_d.x
        print("Z_0 = ")
        tsp_args = (s, p, g_x_func, True)

        if method_eq == 'tgo':
            res_eq = tgo(self.tsp_objective_function, bounds=bounds,
                         args=tsp_args)
        else:
            res_eq = scipy.optimize.minimize(self.tsp_objective_function,
                                             Z_0,
                                             args=tsp_args,
                                             method=method_eq
                                             )
        print('=' * 100)
        logging.info(res_eq)

        return res_eq

    def vdw_dwpm_params(self, params, p):
        """
        This function accepts a params vector and feeds updates to the
        parameter class of the VdW-dwpm EOS (TODO: Move to van_der_waals module

        DWPM-EOS: [p.m['r'], p.m['s'],
                   p.m['k'][1][i], ... for all i \in N <= p.m['n']
                   p.m['k'][2][i], ... for all i \in N <= p.m['n']
                   ...
                   p.m['k'][n][i], ... for all i \in N <= p.m['n']
                  ]
        """
        import logging
        if self.rs or self.rskij:
            p.m['r'], p.m['s'] = params[0], params[1]
            #print('r = {}'.format(params[0]))
            #print('s = {}'.format(params[1]))

            if abs(p.m['r']) <= 1e-10:  # Avoid singularities
                logging.warning("r parameter close to singularity, "
                                "setting to r"
                                " = 1e-10")
                p.m['r'] = 1e-10
            if abs(p.m['s']) <= 1e-10:
                logging.warning("s parameter close to singularity, "
                                "setting to s"
                                " = 1e-10")
                p.m['s'] = 1e-10

        if self.kij or self.rskij:
            if self.rs or self.rskij:
                pint = 2
            else:
                pint = 0
            for i in range(1, p.m['n'] + 1):
                for j in range(1, p.m['n'] + 1):
                    if i == j:
                        pass
                    else:
                        p.m['k'][i][j] = params[pint]
                        #print('k_{}{} = {}'.format(i, j, p.m['k'][i][j]))
                        #if abs(1 - p.m['k'][i][j]) <= 1e-10:  # Avoid singularities
                        #    logging.warning(
                        #        "k_{0}{1} parameter close to singularity, "
                         #       "setting to k_{0}{1}"
                        #        " = 0.999".format(i, j))
                        #   p.m['k'][i][j] = 0.999
                        pint += 1
                        #print('p[k] =')
                        #print(p.m['k'])
        return p

    def d_points(self, N, X_I, X_II):  # Validated for binary

        """
        Generates a set of N \in R^n vector points uniformly distributed on the
        directional vector between data points X_I and X_II in the unstable
        region

        These sampling points are used in the objective function to find the
        "degree of convexity" between the solution plane and the Gibbs free
        energy surface.

        TODO: Use Sobol seq. stretched over X_I and X_II bounds
        """
        import numpy
        X_D = []
        #TODO: Use list of equilibrium point vectors and iterate
        #X_sum = X_I + X_II
        X_diff = X_I - X_II
        for n in range(N):
            #X_D.append((n + 1) * X_sum / (N + 1))
            X_D.append(((n + 1) / (N + 1.0)) * X_diff + X_II)

        return X_D


    def o_points(self, N, X_I, X_II):  # Validated for binary

        """
        Generates a set \in R^n vector points uniformly distributed on the
        directional vector between data points X_I and X_II in the stable
        region.

        These sampling points are used in the objective function to find the
        "degree of instability in the stable regions" between the solution
        plane and the Gibbs free energy surface.

        TODO: Use Sobol seq. stretched over X bounds
        """
        import numpy
        X_o = []
        X_diff = X_I - X_II
        X_I_u_diff = 1.0 - X_I
        X_II_u_diff = 1.0 - X_II
        X_I_l_diff = X_I
        X_II_l_diff = X_II

        # Generate stable points from X_I to bound
        for n in range(N):
            point = []
            for i in range(len(X_diff)):
                if X_diff[i] < 0:
                    point_i = -((n + 1) / (N + 1.0)) * X_I_l_diff[i] + X_I[i]
                    point.append(point_i)
                if X_diff[i] > 0:
                    point_i = ((n + 1) / (N + 1.0)) * X_I_u_diff[i] + X_I[i]
                    point.append(point_i)

            point = numpy.array(point)

            X_o.append(point)

        # Generate stable points from X_II to bound
        for n in range(N):
            point = []
            for i in range(len(X_diff)):
                if X_diff[i] < 0:
                    point_i = ((n + 1) / (N + 1.0)) * X_II_u_diff[i] + X_II[i]
                    point.append(point_i)
                if X_diff[i] > 0:
                    point_i = -((n + 1) / (N + 1.0)) * X_II_l_diff[i] + X_II[i]
                    point.append(point_i)

            point = numpy.array(point)

            X_o.append(point)

        return X_o

    def d_plane(self, g_x_func, s, p, X_I, X_II):  # Validated
        """
        Generates a solution plane function that can be used to find a
        scalar output value of a
        _ at any point X.
        """
        import numpy
        # Find gibbs surface values at each solution point
        s.update_state(s, p, X=X_I, Force_Update=True)
        G_I = g_x_func(s, p).m['g_mix']  # Solution of all phases
        G_sol_I = G_I['t']  # minimum Gibbs energy of all phase surfaces
        s.update_state(s, p, X=X_II, Force_Update=True)
        G_II = g_x_func(s, p).m['g_mix']
        G_sol_II = G_II['t']
        Lambda_sol_est = []
        X_I = numpy.array(X_I)
        X_II = numpy.array(X_II)

        #print('G_sol_I = {}'.format(G_sol_I))
        #print('G_sol_II = {}'.format(G_sol_II))
        # Find duality multipliers at set
        if len(X_I) > 1:
            for x_i, x_ii in X_I, X_II:
                Lambda_sol_est.append((G_sol_II - G_sol_I) / (x_ii - x_i))
        else:
            Lambda_sol_est = [(G_sol_II - G_sol_I) / (X_II - X_I)]

        # Define dual lagrange function on solutions:
        def d_plane_sol(X):
            return G_sol_I + sum(Lambda_sol_est * (X - X_I))

        G_sol = [G_sol_I, G_sol_II]
        G_all = [G_I, G_II]
        return d_plane_sol, Lambda_sol_est, G_sol, G_all

    def d_Z_0(self, X_I, X_II):
        """
        Estimate a feeding point between two equilibrium points
        """
        return (X_I + X_II)/2.0

    def dual_gap(self, g_x_func, plane, X_D, s, p):
        """
        Finds the duality gap between the solution plane and the Gibbs surface
        for each sequenced point in X_D, returns vector of len(X_D)
        """
        # Find the difference between the dual solution plane and the
        # gibbs surface values at each solution point
        f_dual_gap = []
        for X in X_D:
            s.update_state(s, p, X=X, Force_Update=True)
            G = g_x_func(s, p).m['g_mix']
            G_surface = G['t']
            #TODO: Check normalization with stable regions
            #norm = max(abs(G_surface), abs(plane(X)))
            #print('G_surface = {}'.format(G_surface))
            #print('plane(X) = {}'.format(plane(X)))
            f_dual_gap.append(plane(X) - G_surface)

        return f_dual_gap


    def dual_gap_error_sum(self, f_dual_gap):
        """
        Add an error if the dual plane is not below the Gibbs surface at
        all sampled points
        """
        import numpy
        epsilon_d = 0.0
        for fdg in f_dual_gap:
            # If the duality gap does not exist/convexity
            try:
                if numpy.float(fdg) > 0.0:
                    epsilon_d += fdg
            except ValueError:
                if numpy.float(fdg[0]) > 0.0:
                    epsilon_d += numpy.float(fdg[0])
        # (If duality gap exists/concavity at the point then we add no penalty)
        return epsilon_d

    def surface_gap_error_sum(self, f_dual_gap_s):
        """
        Add an error if the dual plane is above the Gibbs surface at any
        sampled point
        """
        import numpy
        epsilon_s = 0.0
        for fdg in f_dual_gap_s:
            # If the duality gap does not exist/convexity
            try:
                if numpy.float(fdg) > 0.0:
                    epsilon_s += abs(fdg)
            except ValueError:
                if numpy.float(fdg[0]) > 0.0:
                    epsilon_s += numpy.float(abs(fdg[0]))
        # (If dual plane is BELOW the Gibbs surface no error is added
        return epsilon_s

    def phase_error(self, G_all, ph):
        """
        Returns an error if the lowest calculated model phase is not equal to
        the data phase for multiphase systems.

        Parameters
        ----------
        G_sol

        Returns
        -------

        """
        epsilon_ph = 0.0

        if G_all['ph min'] == ph:
            epsilon_ph = 0.0
        else:
            epsilon_ph = abs(G_all[ph] - G_all[G_all['ph min']])
            #print 'G_all[ph] = {}'.format(G_all[ph])
            #print 'G_all[G_all[\'ph min\']] = {}'.format(G_all[G_all['ph
            # min']])

        return epsilon_ph

    def norm_eta_sum(self, X_D, Lambda_sol_est, X_I, X_II, G_sol):
        """
        Normalized difference between each plane function eta at Z_0 \in X_D
        and max(eta)
        """
        epsilon_e = 0.0

        for X in X_D:
            eta_I = G_sol[0] + sum(Lambda_sol_est * (X - X_I))
            eta_II = G_sol[1] + sum(Lambda_sol_est * (X - X_II))
            epsilon_e += abs((eta_I - eta_II)/max(eta_I, eta_II))

        return epsilon_e


    def data_error(self, X_data, ph_data, X_D, g_x_func, s, p):
        """
        Calculate the error for a single data point (by summing each element
        Sigma n-1 at a single phase equilibrium point).

        (should be in loop:
        for ph, ind in zip(p.m['Data phases'],
                   range(len(p.m['Data phases']))):

        """
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import dual_equal
        import numpy
        #Z_0 = sorted(X_D)[len(X_D) // 2]
        # TODO: Estimate in middle of all equilibrium points in X_data
        Z_0 = self.d_Z_0(X_data[0], X_data[1])
        #dual_equal(s, p, g_x_func, Z_0)
        X_eq, g_eq, phase_eq = pec(s, p, g_x_func, Z_0)

        if len(X_eq) < len(X_data):
            epsilon_x = len(X_data)  # Sigma n - 1 elements of data points

        elif len(X_eq) == len(X_data):
            #Epsilon_x = []
            epsilon_x = 0.0
            data_ind = numpy.argsort(ph_data)
            model_ind = numpy.argsort(phase_eq)
            #X_data[data_ind]
            for i in range(len(X_data)):
                epsilon_x += abs(X_data[data_ind[i]] - X_eq[model_ind[i]])
                #for X_sol, ph in zip(X_eq, phase_eq):
            #for X_sol in X_eq:
                #for ph in p.m['Data phases']:
                #    epsilon_x +=
                #for ph_eq in phase_eq:
                #     for ph_dat in ph_data:
                #         if ph_eq == ph_dat:
                #             for i in range(p.m['n']):
                #                 Epsilon_x.append(abs(X_sol[i] - X_data[i]))
        #
        # elif len(X_eq) > len(X_data):
        #     Epsilon_x = []
        #     for X in X_eq:
        #         for X_sol, ph in zip(X, phase_eq):
        #             if ph == ph_data:
        #                 for i in range(p.m['n']):
        #                     Epsilon_x.append(abs(X_sol[i] - X_data[i]))
        #     epsilon_x = min(Epsilon_x)

        elif len(X_eq) > len(X_data):
            epsilon_x = len(X_data)

        return epsilon_x


    def tsp_objective_function(self, params, s, p, g_x_func,
                               dp_pec=True,
                               dual_s=True,
                               N=7):
        """
        Objective function to minimize, accept some parameter set as first
        argument and calculated error over range of data points.
        """
        import numpy
        import logging
        # Update parameters
        self.param_func(params, p)

        Epsilon = 0.0
        # Error weights (TODO Assess need):
        #  a = 1.0  # Duality gap does not exist errors

        l_d = float(len(p.m['P']))
        l_ph = float(len(p.m['Data phases']))
       # a = 1/180.0 # Stable surfaces
        a = 1#/(N*2.0) # Stable surfaces
        b = 1.0/l_d # 1.0  # Lagrangian plane errors
        c = 1.0/(l_d*l_ph) # Equilibrium point errors

        # Stores for plots
        self.Epsilon_d = 0.0
        self.Epsilon_s = 0.0
        self.Epsilon_ph = 0.0
        self.Epsilon_e = 0.0
        self.Epsilon_x = 0.0

        # Loop through all data points:
        for i in range(len(p.m['T'])):
            p.m['T'][i]

            try: #TODO: Deal with failures in Vroot here
                s.update_state(s, p, P=p.m['P'][i], T=p.m['T'][i],
                               Force_Update=True)
            except numpy.linalg.linalg.LinAlgError:  # , IndexError):
                pass

            # Loop through all phases in equilibrium to find points in
            # equilibrium
            X_eq_data = []  # Data container that contains X_I, X_II, ...
            X_ph_data = []  # Data container that contains phases of
                            #  X_I, X_II, ... (used to skip LLE points in
                            # epsilon_ph routines)

            for ph in p.m['Data phases']:
                X_eq_d = []
                for n in range(1, p.m['n']):
                    X_eq_d.append(p.m[ph][n][i])

                X_eq_data.append(numpy.array(X_eq_d))
                X_ph_data.append(ph)

            # Find error at each point
            # TODO: Update this two phase equilibrium to arbitrarily high
            # NOTE: It will probably be better to keep all functions with only
            # two data points as input since the sequenced sampling points is
            # generated on vectors between the points anyway; so instead
            # iterate through all data points.
            X_I = X_eq_data[0]
            X_II = X_eq_data[1]
            ph_I = X_ph_data[0]
            ph_II = X_ph_data[1]

            if (X_II - X_I).all() == 0.0:  # Skip pure points
                continue

            # Generate (Note, this only needs to be done once and saved in a
            # set for the current data points)

            X_D = self.d_points(N, X_I, X_II)
            plane, Lambda_sol_est, G_sol, G_all = self.d_plane(g_x_func, s,
                                                               p, X_I, X_II)

            f_dual_gap_d = self.dual_gap(g_x_func, plane, X_D, s, p)

            # Find dual gap error if it does not exist
            epsilon_d = self.dual_gap_error_sum(f_dual_gap_d)

            # Find surface errors in the stable region
            if dual_s:
                # Generate (Note, this only needs to be done once and saved in
                # a set for the current data points)
                X_O = self.o_points(N, X_I, X_II)
                f_dual_gap_s = self.dual_gap(g_x_func, plane, X_O, s, p)

                # Find dual gap error if it does not exist
                epsilon_s = self.surface_gap_error_sum(f_dual_gap_s)
            else:
                epsilon_s = 0.0

            # Find phase error if ph_data is not <= ph_model
            if len(p.m['Valid phases']) > 1:
                for ph, i in zip([ph_I, ph_II], range(2)):
                    # (might extend to X_I, X_II, ...)
                    if ph not in p.m['Valid phases']:
                        continue  #TODO: Test if working by skipping LLE
                    else:
                        epsilon_ph = self.phase_error(G_all[i], ph)
            else:
                epsilon_ph = 0.0

            if False:
                print('='*100)
                print('X_I = {}'.format(X_I))
                print('X_II = {}'.format(X_II))
                print('X_D = {}'.format(X_D))
                print('X_O = {}'.format(X_O))
                print('='*100)

            # Find dual plane and phase equilibrium error
            if dp_pec:#epsilon_d == 0.0:
                try:
                    #TODO: Add timeouts
                    epsilon_e = self.norm_eta_sum(X_D, Lambda_sol_est,
                                                  X_I, X_II,
                                                  G_sol)
                    epsilon_x = self.data_error([X_I, X_II],
                                                p.m['Data phases'],
                                                X_D, g_x_func, s, p)

                except numpy.linalg.linalg.LinAlgError:#, IndexError):
                    logging.warning("LinAlgError in phase equil calculation"
                                    " setting epsilons to maximum")
                    epsilon_e = 1.0  # (max normalized plane error)
                    epsilon_x = 1.0 * len(p.m['Data phases'])  # (max eq error)

                    # Remove nans from dict
                    s.update_state(s, p, X=X_I, Force_Update=True)

            # elif epsilon_d < 1e-3:
            #     try:
            #         epsilon_e = self.norm_eta_sum(X_D, Lambda_sol_est,
            #                                       X_I, X_II,
            #                                       G_sol)
            #         epsilon_x = self.data_error([X_I, X_II],
            #                                     p.m['Data phases'],
            #                                     X_D, g_x_func, s, p)
            #     except(numpy.linalg.linalg.LinAlgError, IndexError):
            #         logging.warning(
            #             "LinAlgError in phase equil calculation"
            #             " setting epsilons to maximum")
            #         epsilon_e = 1.0  # (max normalized plane error)
            #         epsilon_x = 1.0 * len(
            #             p.m['Data phases'])  # (max eq error)
            #
            #         # Remove nans from dict
            #         s.update_state(s, p, X=X_I, Force_Update=True)

            else:  # if duality gap does not exist set max error for data
                   # point
                epsilon_e = 1.0  #(max normalized plane error)
                epsilon_x = 1.0 * len(p.m['Data phases']) # (max eq error)

            # Convert all nested values to floats:
            epsilon_d = numpy.float(epsilon_d)
            epsilon_s = numpy.float(epsilon_s)
            epsilon_ph = numpy.float(epsilon_ph)
            epsilon_e = numpy.float(epsilon_e)
            epsilon_x = numpy.float(epsilon_x)

            # SUM all errors
            #print 'Epsilon = {}'.format(Epsilon)
            # Epsilon += a * epsilon_e + b * epsilon_e + c * epsilon_x
            #Epsilon += epsilon_d + b * epsilon_e + c * epsilon_x
            #Epsilon += epsilon_d + epsilon_s + b * epsilon_e + c * epsilon_x

            # Sum stores for plots

            self.Epsilon_d += epsilon_d
            self.Epsilon_s += epsilon_s
            self.Epsilon_ph +=  epsilon_ph
            self.Epsilon_e += epsilon_e
            self.Epsilon_x += epsilon_x


        # Normalize
        #max_surf = max(self.Epsilon_d, self.Epsilon_s)
        #a = 1.0 / max_surf
        self.Epsilon_d = self.Epsilon_d
        self.Epsilon_s = a * self.Epsilon_s
        self.Epsilon_ph = self.Epsilon_ph
        self.Epsilon_e = b * self.Epsilon_e #/1.5
        self.Epsilon_x = c * self.Epsilon_x #/1.5

        # SUM all errors
        Epsilon += self.Epsilon_d + self.Epsilon_s  + self.Epsilon_ph\
                   + self.Epsilon_e + self.Epsilon_x

        if True:
            print('r = {}'.format(p.m['r']))
            print('s = {}'.format(p.m['s']))
            print('k_12 = {}'.format(p.m['k'][1][2]))
            print('k_21 = {}'.format(p.m['k'][2][1]))
            print('self.Epsilon_d = {}'.format(self.Epsilon_d))
            print('self.Epsilon_s = {}'.format(self.Epsilon_s))
            print('self.Epsilon_ph = {}'.format(self.Epsilon_ph))
            print('self.Epsilon_e = {}'.format(self.Epsilon_e))
            print('self.Epsilon_x = {}'.format(self.Epsilon_x))
            print('Epsilon = {}'.format(Epsilon))

        return Epsilon

    def obj_func_range(self, func, bounds=[(0.8, 1.2), (0.8, 1.2)], x_r=6,
                       args=(), comps=None, vary='rs'):
        import numpy
        from tinydb import TinyDB
        import logging

        x_range = numpy.linspace(bounds[0][0], bounds[0][1], x_r)
        y_range = numpy.linspace(bounds[1][0], bounds[1][1], x_r)
        xg, yg = numpy.meshgrid(x_range, y_range)
        func_r = numpy.zeros((x_r, x_r))
        func_ed = numpy.zeros((x_r, x_r))
        func_es = numpy.zeros((x_r, x_r))
        func_eph = numpy.zeros((x_r, x_r))
        func_ee = numpy.zeros((x_r, x_r))
        func_ex = numpy.zeros((x_r, x_r))

        for i in range(xg.shape[0]):
            for j in range(yg.shape[0]):
                X = [xg[i, j], yg[i, j]]  # [x_1, x_2]

                f_out = func(X, *args)  # Scalar outputs
                func_r[i, j] = numpy.float64(f_out)
                func_ed[i, j] = numpy.float64(self.Epsilon_d)
                func_es[i, j] = numpy.float64(self.Epsilon_s)
                func_eph[i, j] = numpy.float64(self.Epsilon_ph)
                func_ee[i, j] = numpy.float64(self.Epsilon_e)
                func_ex[i, j] = numpy.float64(self.Epsilon_x)
                if False:
                    print('self.Epsilon_d = {}'.format(self.Epsilon_d))
                    print('self.Epsilon_s = {}'.format(self.Epsilon_s))
                    print('self.Epsilon_ph = {}'.format(self.Epsilon_ph))
                    print('self.Epsilon_e = {}'.format(self.Epsilon_e))
                    print('self.Epsilon_x = {}'.format(self.Epsilon_x))

        # Save results

        plot_kwargs = {'func_r' : func_r,
                       'xg' : xg,
                       'yg' : yg,
                       'func_ed' : func_ed,
                       'func_es' : func_es,
                       'func_eph' : func_eph,
                       'func_ee' : func_ee,
                       'func_ex' : func_ex,
                       }
        if False:
            # TODO: Find an easy way to store multidimensional numpy arrays
            db = TinyDB('.db/p_obj_db.json')
            db_store = {'dtype' : 'p_obj_r',
                        'comps' : comps,
                        'bounds' : bounds,
                        'x_r' : x_r,
                        'vary' : vary,
                        'plot_kwargs' : plot_kwargs}
            db.insert(db_store)
        return plot_kwargs

    def plot_ep(self, plot_kwargs, axis_labels=['r', 's']):
        """
        Plot the speficied single var input error function
        over a range size x_r
        """
        from matplotlib import pyplot as plot
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plot
        from matplotlib import cm
        import numpy
        # Plots
        fig = plot.figure()

        func_r = plot_kwargs['func_r' ]
        xg = plot_kwargs['xg']
        yg = plot_kwargs['yg']
        func_ed = plot_kwargs['func_ed']
        func_es = plot_kwargs['func_es']
        func_eph = plot_kwargs['func_eph']
        func_ee = plot_kwargs['func_ee']
        func_ex = plot_kwargs['func_ex']

        ax = fig.gca(projection='3d')
        X, Y = xg, yg

        # find annotation points
        xa = numpy.max(xg)
        ya = numpy.max(yg)
        za = func_r[-1, -1]

        # Gibbs phase surfaces
        if True:
            print 'numpy.min(Z) = {}'.format(numpy.nanmin(func_r))
            cset = ax.contourf(X, Y, func_r, zdir='z',
                               offset=numpy.nanmin(func_r)-0.05,
                               cmap=cm.coolwarm)
            ax.plot_surface(X, Y, func_r, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon$')
            ax.text(xa, ya, func_r[-1, -1], '$\epsilon$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_r[0, 0], '$\epsilon$')
            ax.plot_surface(X, Y, func_ed, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_D$')
            ax.text(xa, ya, func_ed[-1, -1], '$\epsilon_D$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_ed[0, 0], '$\epsilon_D$')
            ax.plot_surface(X, Y, func_es, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_S$')
            ax.text(xa, ya, func_es[-1, -1], '$\epsilon_S$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_es[0, 0], '$\epsilon_S$')
            ax.plot_surface(X, Y, func_eph, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_{ph}$')
            ax.text(xa, ya, func_eph[-1, -1], '$\epsilon_{ph}$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_eph[0, 0], '$\epsilon_{ph}$')
            ax.plot_surface(X, Y, func_ee, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_e$')
            ax.text(xa, ya, func_ee[-1, -1], '$\epsilon_e$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_ee[0, 0], '$\epsilon_e$')
            ax.plot_surface(X, Y, func_ex, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_x$')
            ax.text(xa, ya, func_ex[-1, -1], '$\epsilon_x$')
            ax.text(numpy.min(xg), numpy.min(xg),
                    func_ex[0, 0], '$\epsilon_x$')

        if False:
            surf = ax.plot_surface(X, Y, func_r, rstride=1, cstride=1,
                                   cmap=cm.coolwarm, linewidth=0,
                                   antialiased=True, alpha=0.5,
                                   shade = True)

            fig.colorbar(surf, shrink=0.5, aspect=5)


        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        if self.rs:
            ax.set_xlabel('$r$')
            ax.set_ylabel('$s$')

        if self.kij:
            ax.set_xlabel('$k_{12}$')
            ax.set_ylabel('$k_{21}$')

        ax.set_zlabel('$\epsilon$', rotation=90)
        #ax.legend()
        plot.show()

        return



if __name__ == '__main__':
    for i in range(10):
        if i == 5:
            continue








