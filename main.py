#!/usr/bin/env python
"""
python2 main.py -c acetone benzene water -p x y
python2 -i main.py -c acetone benzene water -p x y -T 273.15 -P 101e3 -r 1.0 -s 1.0

python2 main.py -c acetone benzene water -p x y -T 273.15 -P 101e3 -r 1.0 -s 1.0

"""

import data_handling
import nComp
import pure
import plot
import logging
import os
import numpy
import scipy
import scipy.interpolate
import Van_der_Waals
VdW = Van_der_Waals.VdW()
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
    parser.add_argument('-kij', '--k_params', nargs='+',
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
                        help='Calculate oly phase seperation of same volume '
                             'root')

    # Plots
    parser.add_argument('-pltg', '--plot_gibbs',
                        action="store_true",
                        help='plot gibbs energy (binary and ternary systems '
                             'only)')

    parser.add_argument('-plti', '--plot_iso',
                        action="store_true",
                        help='plot phase seperations (binary and ternary '
                             'systems only) the most appropriate plot '
                             '(isotherms vs isobars) is determined from the'
                             ' data')

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
        from nComp import phase_equilibrium_calculation as pec
        from nComp import phase_seperation_detection as psd
        from nComp import equilibrium_range as er
        from nComp import g_mix as g_x_func
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()
        s, p = nComp.n_comp_init(data)

        # Parameter optimisation
        if data.optimise:
            p.m['r'] = 1.0
            p.m['s'] = 1.0
            p.m['k'][1][2] = 0.124
            p.m['k'][2][1] = 0.124

            ph_eq, mph_eq, mph_ph = psd(g_x_func, s, p,
                                        P=24e5, T=263.1,
                                        n=100,
                                       VLE_only=True,
                                       Plot_Results=False)

            P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph = er(g_x_func, s, p,
                                                               Data_Range=True)
            print 'P_range ='
            print P_range
            print 'T_range ='
            print T_range
            print 'r_ph_eq ='
            print r_ph_eq
            print 'r_mph_eq ='
            print r_mph_eq
            print 'r_mph_ph ='
            print r_mph_ph
            pass #TODO

        # Simulate specifications
        if data.P is not None and data.T is not None and data.Z_0 is None:
            psd(g_x_func, s, p, data.P, data.T, n=100, LLE_only=data.lle_only,
                                   VLE_only=data.vle_only) # Tested/working

        if data.P is not None and data.T is not None and data.Z_0 is not None:
            pec(s, p, g_x_func, Z_0, k=None, P=data.P, T=data.T, # Not tested
               tol=1e-9, Print_Results=True, Plot_Results=data.plot_gibbs)

        # TODO Move to plot.py and exec script?
        if data.plot_gibbs:
            options = plot.plot_options
            plot.plot_g_mix(s, p, options, figno=None)

        if data.plot_iso:
            pass #TODO