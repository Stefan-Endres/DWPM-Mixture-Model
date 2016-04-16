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
import plot
import logging
import os
import numpy
import scipy
import scipy.interpolate
from models import van_der_waals
VdW = van_der_waals.VdW()
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
                             'used if new vapour datais added')

    args = parser.parse_args()
    data.run_options(args)


    if len(data.comps) == 1:  # pure component simulation.
        # Load pure data
        data.load_pure_data() # Using data.comps

        # Find all specified outputs
        s, p = pure.pure_sim(data, i=0)

    if len(data.comps) > 1:  # multi component simulation.
        from ncomp import phase_equilibrium_calculation as pec
        from ncomp import phase_seperation_detection as psd
        from ncomp import equilibrium_range as er
        from ncomp import g_mix as g_x_func
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()
        s, p = ncomp.n_comp_init(data)

        # Parameter optimisation
        if data.optimise:
            p.m['r'] = 1.0
            p.m['s'] = 1.0
            p.m['k'][1][2] = 0.124
            p.m['k'][2][1] = 0.124
            # pass #TODO

            from tgo import tgo
            from ncomp import parameter_goal_func as pgf
            Bounds = [(-5.0, 5.0),  # r
                      (-5.0, 5.0)]#,  # s
                     # (0.0, 0.999),   # k12
                     # (0.0, 0.999),   # k21
                     # ]

            #p.m['r'], p.m['s'] = -3.75,  1.25
           # p.m['k'][1][2] = 0.874125
            #p.m['k'][2][1] = 0.874125

            #Params =  [-3.75, 1.25, 0.874125, 0.874125]

            #print pgf( Params, g_x_func, s, p, 200, # n
            #                     False,
             #                    True),def

            optimres = tgo(pgf, Bounds,
                           args=(g_x_func, s, p,
                                 200, # n
                                 False,
                                 True),  # VLE only
                                 #g_func=x_lim,
                                 n = 100,
                                 skip=2)

            print optimres.x


        # Simulate specifications
        if data.P is not None and data.T is not None and data.Z_0 is None:
            psd(g_x_func, s, p, data.P, data.T, n=100, LLE_only=data.lle_only,
                                   VLE_only=data.vle_only,
                                   Plot_Results=False) # Tested/working


        if data.P is not None and data.T is not None and data.Z_0 is not None:
            pec(s, p, g_x_func, data.Z_0, k=None, P=data.P, T=data.T,
               tol=1e-9, Print_Results=True, Plot_Results=data.plot_gibbs)
                # Not tested


        if data.plot_gibbs:
            options = plot.plot_options
            plot.plot_g_mix(s, p, options, figno=None)

        if data.plot_isotherms is not None:
            from ncomp import g_mix as g_x_func
            from plot import Iso
            iso = Iso()
            import time
            start = time.time()
            iso.plot_iso(s, p, g_x_func, res=4, n=1000, T=data.plot_isotherms,
                         VLE_only=True, n_dual=300, Plot_Results=True)
            print("="*90)
            print('Done in {}'.format(time.time() - start))
            print("="*90)

            from matplotlib import pyplot as plot

            plot.show()
