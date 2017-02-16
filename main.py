#!/usr/bin/env python
"""
python2 main.py -c acetone benzene water -p x y
python2 -i main.py -c acetone benzene water -p x y -T 273.15 -P 101e3 -r 1.0 -s 1.0

python2 main.py -c acetone benzene water -p x y -T 273.15 -P 101e3 -r 1.0 -s 1.0

-c carbon_dioxide ethane -p x y -P 24e5 -T 263.1 -r 1.0 -s 1.0 -plti -kij 0.124 0.124

"""

import data_handling
import ncomp
import pure

import argparse

if __name__ == '__main__':
    # Return basic data
    data = data_handling.ImportData()

    parser = argparse.ArgumentParser('Test argument parsing')
    # Positional arguments
    parser.add_argument('-c', '--compounds', nargs='+', required=True,
                        help='compounds to simulate for example'
                             'acetone benzene water')
    parser.add_argument('-p', '--phases', nargs='+', required=True,
                        help='List of valid phases in equilibrium'
                             'ex. for VLE use x y')
    # Optional arguments
    parser.add_argument('-e', '--eos', default='DWPM',
                        choices=['DWPM'],
                        help='Equation of State / Mixture rule')
    parser.add_argument('-r', type=float,
                        help='Force value of r')
    parser.add_argument('-s', type=float,
                        help='Force value of s')
    parser.add_argument('-m', '--model', nargs=1,
                        default="Adachi-Lu",
                        choices=['Adachi-Lu', 'Soave'],
                        help='Actvity coefficient model')

    #TODO Add k_ij array input and processing
    parser.add_argument('-kij', '--k_params', nargs='+', type=float,
                        help='Force value of interaction parameters')

    parser.add_argument('-T', '--temperature', type=float,
                        help='Temperature for point simulation')
    parser.add_argument('-P', '--pressure', type=float,
                        help='Pressure for point simulation')
    parser.add_argument('-z', type=float, nargs="+",
                        help='Composition for point simulation')

    #
#    parser.add_argument('-g_f', '--g_x_func' type=float, nargs="+",
#                        help='')

    parser.add_argument('-vle', '--vle_only',
                        action="store_true",
                        help='If specified then phase seperation of same '
                             'volume root instability will be ignored.')

    parser.add_argument('-lle', '--lle_only',
                        action="store_true",
                        help='Calculate only phase seperation of same volume '
                             'root')

    # Plots
    parser.add_argument('-pltg', '--plot_gibbs',
                        action="store_true",
                        help='plot gibbs energy (binary and ternary systems '
                             'only)')

    parser.add_argument('-pltit', '--plot_isotherms', nargs='+', type=float,
                        help='plot isotherm phase envelope (binary and '
                             'ternary systems only)')

    parser.add_argument('-pltib', '--plot_isobars', nargs='+', type=float,
                        help='plot isobar phase envelope (binary and '
                             'ternary systems only)')

    parser.add_argument('-pltp', '--plot_pure',
                        action="store_true",
                        help='Plot the pure vapour pressure model')

    parser.add_argument('-epsp', '--plot_epsilon',
                        action="store_true",
                        help='Plot the errors of the objective function')

    # Optimise
    parser.add_argument('-opt', '--optimise',
                        action="store_true",
                        help='Optimise the DWPM parameters')

    #  Save
    parser.add_argument('-save', nargs=1, type=bool,
                        default=False,
                        help='Save the results of the multi-component'
                             'optimisation')

    parser.add_argument('-save_pure', nargs=1, type=bool,
                        default=False,
                        help='Save the results of the pure component'
                             'optimisation')

    parser.add_argument('-force_pure_update', nargs=1, type=bool,
                        default=False,
                        help=' force a new optimisation for the m'
                             'parameter for the selected Model, to be '
                             'used if new vapour data is added')

    args = parser.parse_args()
    data.run_options(args)


    if len(data.comps) == 1:  # pure component simulation.
        # Load pure data
        data.load_pure_data() # Using data.comps

        # Find all specified outputs
        s, p = pure.pure_sim(data, i=0)

        # plot output
        if data.plot_pure:
            from plot import PsatPlots
            PP = PsatPlots(data.comps[0], data.model[0])

            if len(PP.DBr) == 0:
                P_sat_store, T_sat_store = PP.psat_range(s, p)

            PP.plot_Psat(data.comps[0], p)

    if len(data.comps) > 1:  # multi component simulation.
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import phase_seperation_detection as psd
        from ncomp import equilibrium_range as er
        # TODO: Select EOS model from input:
        from ncomp import g_mix as g_x_func
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()
        s, p = ncomp.n_comp_init(data)


        # Parameter optimisation
        if data.optimise:
            from param import TopShiftParam
            from ncomp import g_mix
            import numpy

            # Local+global routine
            s.update_state(s, p, P=24e5, T=263.1, X=[0.0],
                           Force_Update=True)

            #print('r = {}'.format(p.m['r']))
            #print('s = {}'.format(p.m['s']))
            TSP = TopShiftParam(p,
                                rs=True,
                                kij=True,
                                #rskij=True
                                )

            tsp_args = (s, p, g_mix)
            tsp_args = (s, p, g_mix, False, True, 6)

            Z_0 = [p.m['k'][1][2], p.m['k'][2][1]]
            Z_0 = [p.m['r'], p.m['s'], p.m['k'][1][2], p.m['k'][2][1]]

            ### Dev
            if False:
                TSP.tsp_objective_function(Z_0, s, p, g_x_func,
                                           dp_pec=False)
            ###

            #if False:
            TSP.optimise(s, p, g_x_func, Z_0,
                         method_d='tgo',
                         #method_d='L-BFGS-B',
                         #method_eq='L-BFGS-B',
                         method_eq='SLSQP',
                         bounds=[
                                 (-20.0, 20.0),
                                 (-20.0, 20.0),
                                 (-0.9, 0.99),
                                 (-0.9, 0.99),
                             #(-4.0, 5.0)
                                 ])

            if False:
                #TODO: Move this to a unittest
                s.update_state(s, p, P=24e5, T=263.1, X=[0.0], Force_Update=True)

                TSP = TopShiftParam(p)

                X_I = numpy.array([0.1939063])  # 'x'
                X_II = numpy.array([0.308988493])  # 'y'
                # X_I = numpy.array([0.1939063, 0.5])  # 'x'
                # X_II = numpy.array([0.308988493, 0.1])  # 'y'

                params = [1.0, 1.0]  # r and s
                TSP.vdw_dwpm_params(params, p)
                X_D = TSP.d_points(5, X_I, X_II)
                X_o = TSP.o_points(5, X_I, X_II)
                print 'X_o = {}'.format(X_o)

                plane, Lambda_sol_est, G_sol = TSP.d_plane(g_mix, s, p, X_I, X_II)
                f_dual_gap = TSP.dual_gap(g_mix, plane, X_D, s, p)
                epsilon_d = TSP.dual_gap_error_sum(f_dual_gap)
                print 'epsilon_d = {}'.format(epsilon_d)
                epsilon_e = TSP.norm_eta_sum(X_D, Lambda_sol_est, X_I, X_II, G_sol)
                epsilon_x = TSP.data_error([X_I, X_II], ['x', 'y'],
                                           X_D, g_mix, s, p)

                print 'epsilon_x = {}'.format(epsilon_x)
                Z_0 = TSP.d_Z_0(X_I, X_II)
                print 'Z_0 = {}'.format(Z_0)

        if data.plot_epsilon:
            from param import TopShiftParam
            from ncomp import g_mix
            import numpy
            s.update_state(s, p, P=24e5, T=263.1, X=[0.0],
                           Force_Update=True)

            # Plot
            tsp_args = (s, p, g_mix, False)
            TSP = TopShiftParam(p, rs=True, kij=False)
            #TSP = TopShiftParam(p, rs=False, kij=True)
            tsp_args = (s, p, g_mix, False, True, 5)
            bounds = [(-5.0, 4.0), (-5.0, 4.0)]
            #bounds = [(-2.0, 2.0), (-2.0, 2.0)]
            #bounds = [(-0.15, 0.99), (-0.15, 0.99)]
            x_r = 40


            plot_kwargs = TSP.obj_func_range(TSP.tsp_objective_function,
                                             bounds, x_r, tsp_args,
                                             comps=data.comps)

            TSP.plot_ep(plot_kwargs)


        # Simulate specifications
        if data.P is not None and data.T is not None and data.Z_0 is None:
            ph_eq, mph_eq, mph_ph = \
                psd(g_x_func, s, p, data.P, data.T, n=100,
                    #n_dual=1,
                    LLE_only=data.lle_only,
                    VLE_only=data.vle_only,
                    Plot_Results=True) # Tested/working

            print('ph_eq = {}'.format(ph_eq))
            print('mph_eq = {}'.format(mph_eq))
            print('mph_ph = {}'.format(mph_ph))

        if data.P is not None and data.T is not None and data.Z_0 is not None:
            pec(s, p, g_x_func, data.Z_0, k=None, P=data.P, T=data.T,
               tol=1e-9, Print_Results=True, Plot_Results=data.plot_gibbs)
                # Not tested


        if data.plot_gibbs: #TODO: No need for this if we plot in pec?
            pass # Need to add tie lines
            #options = plot.plot_options
            #plot.plot_g_mix(s, p, options, figno=None)
            #plot.plot_g_mix(s, p, g_x_func)#, figno=None)

        if data.plot_isotherms is not None:
            from ncomp import g_mix as g_x_func
            from plot import IsoDetection
            iso = IsoDetection(components=data.comps)
            import time
            start = time.time()
            # Use res = 30 for database standard
            iso.plot_iso(s, p, g_x_func, res=40, n=3000, T=data.plot_isotherms,
                         VLE_only=True, n_dual=300, Plot_Results=False)
            print("="*90)
            print('Done in {}'.format(time.time() - start))
            print("="*90)

            from matplotlib import pyplot as plot

            #plot.show()

        from matplotlib import pyplot as plot
        plot.show()