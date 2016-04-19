#!/usr/bin/env python

class TopShiftParam:
    def __init___(self):
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import dual_equal
        pass

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
        p.m['r'], p.m['s'] = params[0], params[1]

        if abs(p.m['r']) <= 1e-10:  # Avoid singularities
            logging.warning("r parameter close to singularity, setting to r"
                            " = 1e-3")
            p.m['r'] = 1e-3
        if abs(p.m['s']) <= 1e-10:
            logging.warning("s parameter close to singularity, setting to s"
                            " = 1e-3")
            p.m['s'] = 1e-3

        # for i in range(1, p.m['n'] + 1):
        #     for j in range(1, p.m['n'] + 1):
        #         if i == j:
        #             pass
        #         else:
        #             p.m['k'][i][j] = Params[pint]
        #             pint += 1
        return p

    def d_points(self, N, X_I, X_II):  # Validated for binary

        """
        Generates a set of N \in R^n vector points uniformly distributed on the
        directional vector between data points X_I and X_II.

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

    def d_plane(self, g_x_func, s, p, X_I, X_II):  # Validated
        """
        Generates a solution plane function that can be used to find a
        scalar output value of a
        _ at any point X.
        """
        import numpy
        # Find gibbs surface values at each solution point
        s.update_state(s, p, X=X_I, Force_Update=True)
        G_sol_I = g_x_func(s, p).m['g_mix']['t']
        s.update_state(s, p, X=X_II, Force_Update=True)
        G_sol_II = g_x_func(s, p).m['g_mix']['t']
        Lambda_sol_est = []
        X_I = numpy.array(X_I)
        X_II = numpy.array(X_II)

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
        return d_plane_sol, Lambda_sol_est, G_sol


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
            f_dual_gap.append(plane(X) - g_x_func(s, p).m['g_mix']['t'])

        return f_dual_gap


    def dual_gap_error_sum(self, f_dual_gap):
        import numpy
        epsilon_d = 0.0
        for fdg in f_dual_gap:
            # If the duality gap does not exist/convexity
            if numpy.float(fdg) > 0.0:
                epsilon_d += fdg

        # (If duality gap exists/concavity at the point then we add no penalty)
        return epsilon_d

    def norm_eta_sum(self, X_D, Lambda_sol_est, X_I, X_II, G_sol):
        """
        Normalized difference between each plane function eta at Z_0 \in X_D
        and max(eta)
        """
        epsilon_e = 0.0

        for X in X_D:
            eta_I = G_sol[0] + sum(Lambda_sol_est * (X - X_I))
            eta_II = G_sol[1] + sum(Lambda_sol_est * (X - X_II))
            epsilon_e += (eta_I - eta_II)/max(eta_I, eta_II)

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
        Z_0 = sorted(X_D)[len(X_D) // 2]
        #
        #dual_equal(s, p, g_x_func, Z_0)

        X_eq, g_eq, phase_eq = pec(s, p, g_x_func, Z_0)


        print X_data
        print 'X_eq = {}'.format(X_eq)
        print 'phase_eq  = {}'.format(phase_eq )

        #for



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


    def tsp_objective_function(self):
        pass














