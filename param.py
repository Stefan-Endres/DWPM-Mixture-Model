#!/usr/bin/env python

class TopShiftParam:
    def __init___(self):
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

        return d_plane_sol


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

    def norm_eta(self):
        """
        Normalized difference between each plane function eta at Z_0 \in X_D
        """
        pass




















