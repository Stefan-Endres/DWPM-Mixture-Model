#!/usr/bin/env python

class TopShiftParam:
    def __init__(self, p):
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import dual_equal
        if p.m['Model'] == 'DWPM':
            self.param_func = self.vdw_dwpm_params

    def vdw_dwpm_params(self, params, p, rs=True, kij=False):
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
        if rs:
            p.m['r'], p.m['s'] = params[0], params[1]

            if abs(p.m['r']) <= 1e-10:  # Avoid singularities
                logging.warning("r parameter close to singularity, "
                                "setting to r"
                                " = 1e-3")
                p.m['r'] = 1e-3
            if abs(p.m['s']) <= 1e-10:
                logging.warning("s parameter close to singularity, "
                                "setting to s"
                                " = 1e-3")
                p.m['s'] = 1e-3

        if kij:
            if rs:
                pint = 2
            else:
                pint = 0
            for i in range(1, p.m['n'] + 1):
                for j in range(1, p.m['n'] + 1):
                    if i == j:
                        pass
                    else:
                        p.m['k'][i][j] = params[pint]
                        print('k_{}{} = {}'.format(i, j, p.m['k'][i][j]))
                        #if abs(1 - p.m['k'][i][j]) <= 1e-10:  # Avoid singularities
                        #    logging.warning(
                        #        "k_{0}{1} parameter close to singularity, "
                         #       "setting to k_{0}{1}"
                        #        " = 0.999".format(i, j))
                        #   p.m['k'][i][j] = 0.999
                        pint += 1
                        print('p[k] =')
                        print(p.m['k'])
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
            try:
                if numpy.float(fdg) > 0.0:
                    epsilon_d += fdg
            except ValueError:
                if numpy.float(fdg[0]) > 0.0:
                    epsilon_d += numpy.float(fdg[0])
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
        Z_0 = sorted(X_D)[len(X_D) // 2]
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


    def tsp_objective_function(self, params, s, p, g_x_func, dp_pec=True):
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
        b = 1e-2  # 1.0  # Lagrangian plane errors
        c = 5e-3  # 2.0  # Equilibrium point errors

        l_d = float(len(p.m['P']))
        l_ph = float(len(p.m['Data phases']))
        b = 1.0/l_d # 1.0  # Lagrangian plane errors
        c = 1.0/(l_d*l_ph)  # 2.0  # Equilibrium point errors

        # Stores for plots
        self.Epsilon_d = 0.0
        self.Epsilon_e = 0.0
        self.Epsilon_x = 0.0

        # Loop through all data points:
        for i in range(len(p.m['T'])):
            p.m['T'][i]
            s.update_state(s, p, P=p.m['P'][i], T=p.m['T'][i],
                           Force_Update=True)

            # Loop through all phases in equilibrium to find points in
            # equilibrium
            X_eq_data = []  # Data container that contains X_I, X_II, ...
            for ph in p.m['Data phases']:
                X_eq_d = []
                for n in range(1, p.m['n']):
                    X_eq_d.append(p.m[ph][n][i])

                X_eq_data.append(numpy.array(X_eq_d))

            # Find error at each point
            # TODO: Update this two phase equilibrium to arbitrarily high
            X_I = X_eq_data[0]
            X_II = X_eq_data[1]

            if (X_II - X_I).all() == 0.0:  # Skip pure points
                continue

            # Generate (Note, this only needs to be done once and saved in a
            # set for the current data points)
            X_D = self.d_points(5, X_I, X_II)
            plane, Lambda_sol_est, G_sol = self.d_plane(g_x_func, s,
                                                        p, X_I, X_II)
            f_dual_gap = self.dual_gap(g_x_func, plane, X_D, s, p)

            # Find dual gap error if it does not exist
            epsilon_d = self.dual_gap_error_sum(f_dual_gap)

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
            epsilon_e = numpy.float(epsilon_e)
            epsilon_x = numpy.float(epsilon_x)

            # SUM all errors
            #print 'Epsilon = {}'.format(Epsilon)
            # Epsilon += a * epsilon_e + b * epsilon_e + c * epsilon_x
            Epsilon += epsilon_d + b * epsilon_e + c * epsilon_x

            # Store for plots
            self.Epsilon_d += epsilon_d
            self.Epsilon_e += b * epsilon_e
            self.Epsilon_x += c * epsilon_x

        if True:
            print('r = {}'.format(p.m['r']))
            print('s = {}'.format(p.m['s']))
            print('k_12 = {}'.format(p.m['k'][1][2]))
            print('k_21 = {}'.format(p.m['k'][2][1]))
            print('self.Epsilon_d = {}'.format(self.Epsilon_d))
            print('self.Epsilon_e = {}'.format(self.Epsilon_e))
            print('self.Epsilon_x = {}'.format(self.Epsilon_x))
            print('Epsilon = {}'.format(Epsilon))

        return Epsilon

    def plot_ep(self, func,
                bounds=[(0.8, 1.2), (0.8, 1.2)], x_r=6, args=()):
        """
        Plot the speficied single var input error function
        over a range size x_r
        """
        from matplotlib import pyplot as plot
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plot
        from matplotlib import cm
        import numpy
        x_range = numpy.linspace(bounds[0][0],bounds[0][1], x_r)
        y_range = numpy.linspace(bounds[1][0],bounds[1][1], x_r)
        xg, yg = numpy.meshgrid(x_range, y_range)
        func_r = numpy.zeros((x_r, x_r))
        func_ed = numpy.zeros((x_r, x_r))
        func_ee = numpy.zeros((x_r, x_r))
        func_ex = numpy.zeros((x_r, x_r))

        for i in range(xg.shape[0]):
            for j in range(yg.shape[0]):
                X = [xg[i, j], yg[i, j]]  # [x_1, x_2]

                f_out = func(X, *args)  # Scalar outputs
                func_r[i, j] = numpy.float64(f_out)
                func_ed[i, j] = numpy.float64(self.Epsilon_d)
                func_ee[i, j] = numpy.float64(self.Epsilon_e)
                func_ex[i, j] = numpy.float64(self.Epsilon_x)


        # Plots
        fig = plot.figure()
        ax = fig.gca(projection='3d')
        X, Y = xg, yg

        # Gibbs phase surfaces
        Z = func_r

        if True:
            print 'numpy.min(Z) = {}'.format(numpy.nanmin(Z))
            cset = ax.contourf(X, Y, Z, zdir='z',
                               offset=numpy.nanmin(Z)-0.05,
                               cmap=cm.coolwarm)
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon$')
            ax.plot_surface(X, Y, func_ed, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_D$')
            ax.plot_surface(X, Y, func_ee, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_e$')
            ax.plot_surface(X, Y, func_ex, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm, label='$\epsilon_x$')

        if False:
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                   cmap=cm.coolwarm, linewidth=0,
                                   antialiased=True, alpha=0.5,
                                   shade = True)

            fig.colorbar(surf, shrink=0.5, aspect=5)

        ax.set_xlabel('$r$')
        #ax.set_xlabel('$k_12$')
        #ax.set_xlim(0, 2)
        ax.set_ylabel('$s$')
        #ax.set_ylabel('$k_21$')
        #ax.set_ylim(0, 2)
        ax.set_zlabel('$\epsilon$', rotation=90)
        plot.show()

        return



if __name__ == '__main__':
    for i in range(10):
        if i == 5:
            continue








