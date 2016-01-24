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
import argparse

VdW = Van_der_Waals.VdW()




#%% Define multi-component simulation function
def n_comp_sim(data):
    # Define parameter class
    p = data_handling.MixParameters()
    p.mixture_parameters(data.VLE, data)
    p.m['n'] = len(data.comps)  # Define system size
    for i in range(p.m['n']):  # Set params for all compounds
        p.parameters(data.c[i])  # Defines p.c[i]
        #p.parameter_build(data.c[i])
    p.m['R'] = p.c[1]['R']  # Use a component Universal gas constant

    # %% Initialize state variables
    s = nComp.state()
    s.mixed()  # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    for i in range(1, p.m['n']+1):
        s.pure(p, i)  # Call using ex. s.c[1]['key']

    p.m['R'] = p.c[1]['R']  # Use a component Universal gas constant

    # %% Initialize state variables
    s = nComp.state()
    s.mixed()  # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    for i in range(1, p.m['n']+1):
        s.pure(p, i)  # Call using ex. s.c[1]['key']

    return s, p




if __name__ == '__main__':
    # Return basic data
    data = data_handling.ImportData()

    parser = argparse.ArgumentParser('Test argument parsing')
    # Positional arguments
    parser.add_argument('-c', '--compounds', nargs='+', required=True,
                        help='compounds to simulate for example'
                             'acetone water phenol')
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
    parser.add_argument('-T', '--temperature', nargs=1, type=float,
                        help='Temperature for point simulation')
    parser.add_argument('-P', '--pressure', nargs=1, type=float,
                        help='Pressure for point simulation')
    parser.add_argument('-z', type=float, nargs="+",
                        help='Composition for point simulation')

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
    # Save
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
    # Plotting

    args = parser.parse_args()
    print args.plot_gibbs
    data.run_options(args)
    print data.plot_gibbs

    if len(data.comps) == 1:  # pure component simulation.
        # Load pure data
        data.load_pure_data() # Using data.comps

        # Find all specified outputs
        s, p = pure.pure_sim(data, i=0)

    if len(data.comps) > 1:  # multi component simulation.
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()

        # Load VLE and mixture parameter data
        data.load()

        s, p = n_comp_sim(data)


    #if data.plotting:
    #    from plot import *  # allow easier func calls from python shell