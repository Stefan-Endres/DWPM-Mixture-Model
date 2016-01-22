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

#%% Define pure simulation function
def pure_sim(data, i=0):
    """data = data class, i compound to simulate """
    # TODO: replace all c[0] with c[i]

    # Find model parameters if not defined
    p = data_handling.parameter_build(data.c[i])
    #p = parameters(data.c[0],I)
    ## (data.c[0], data.comps)
    s = {}
     # Find a_c, b_c  parameters if not available and test dimensional values
    if (data.c[i]['a_c (Pa m6 mol-2)'][0] == ''
        or data.c[i]['b_c (m3 mol-1)'][0] == ''):

        try: # Check if all params are available
            data.c[i]['b_c (m3 mol-1)'] = p['R']*p['T_c']/(8.0*p['P_c'])
            data.c[i]['a_c (Pa m6 mol-2)'] = 27*(p['R']**2)*(p['T_c']**2) \
                                             /(64*p['P_c'])
            if not numpy.round(p['b_c'], 8) == \
                   numpy.round(data.c[i]['b_c (m3 mol-1)'], 8):
                logging.warn('Calculated parameter for \'b_c\' does'
                              + 'not match stored data value')
                logging.warn('Changing data value for \'b_c\' from'
                             + ' \'b_c\' = {}'.format(p['b_c'])
                             + ' to \'b_c\' = {}'.format(data.c[0]['b_c'])
                             )

                p['b_c'] = data.c[0]['b_c']

            if not numpy.round(p['a_c'], 8) == \
                   numpy.round(data.c[i]['a_c (Pa m6 mol-2)'], 8):

                logging.warn(' Calculated parameter for \'a_c\' does'
                              + 'not match stored data value.')
                logging.warn(' Changing data value for \'a_c\' from'
                             + ' \'a_c\' = {}'.format(p['a_c'])
                             + ' to \'a_c\' = {}'.format(data.c[i]['a_c'])
                             )

                p['a_c'] = data.c[i]['a_c']

        except(NameError,KeyError):
            raise IOError('Missing \'P_c\' and/or \'T_c\' paramter')

     # Find 'm' parameter in a dependency model parameter if not available
    if data.c[i]['m ({})'.format(data.model)][0] == '' \
    or data.force_pure_update:  # Detect model params
                                # or force optimization if true
        try: # Check if vapour pressure data is available
            data.c[i]['P (Pa)']
            data.c[i]['T (K)']
        except(NameError, KeyError):
            raise IOError('No parameters or vapour pressure data detected.')

        #find_a_m() # Find params if data is available
        p['m'] = pure.optim_a_m(p)
        exec 'data.c[0][\'m ({})\'] = p[\'m\']'.format(data.model)
        data.c[i]['m ({})'.format(data.model)] = p['m']
        logging.info('New parameter found for'
                     ' {} model, m = {}'.format(data.model, p['m']))

    #%% Find phase equilibrium at specified Temperature point (T, V_v and V_l)
    if data.T: # Note that if data.T is > 0 then the boolean is 'True'
        s      = {}
        s['b'] = p['b_c'] # b = b_c
        s['a'] = p['a_c'] # First estimate
        s['T'] = data.T
        try:
            s['P'] = scipy.interpolate.interp1d(p['T'],p['P'])(s['T'])
        except(ValueError):
            raise IOError('Specified temperature {} K is larger than the criti'
                          'cal temperature {} K.'.format(s['T'],p['T_c']))

        s = VdW.Psat_V_roots(s, p)
        out_str =('VLE at {T} K: P_sat = {P} kPa, V_v = {Vv} m3 mol-1,'
               ' V_l = {Vl} m3 mol-1'.format(T  = data.T,
                                  P  = s['P_sat']/1000.0,
                                  Vv = s['V_v'],
                                  Vl = s['V_l']))
        print out_str
        logging.info(out_str)

    #%% Find phase equilibrium at specified Pressure point (P, V_v and V_l)
    try: #TODO:
        if data.P: # Note that if I['P'] is > 0 then the boolean is 'True'
            pass#VdW.Tsat_V_roots(s,p) # NOTE TODO!
    except(KeyError):
        pass

    #%% Save if True
    """ NOTE!!!: Save the data container and define first data entry [0]
    DO NOT save the dictionary container (or test first if save_to_dict.py is
    robust))
    """
    if data.save_pure:
        from csvDict import save_dict_as_csv
        # Order of headings to save in .csv
        Order = ['T (K)', 'P (Pa)', 'T_c (K)', 'P_c (Pa)', 'V_c (m3 mol-1)',
                 'Z_c', 'R (m3 Pa K-1 mol-1)' ,'w' ,'a_c (Pa m6 mol-2)',
                 'b_c (m3 mol-1)', 'm (Adachi-Lu)', 'm (Soave)', 'virialT',
                 'virialB','name']
        # Save path string
        sstr = os.path.join(data.datadir,
                            'Pure_Component',
                            '{}.csv'.format(data.c[i]['name'][0]))

        #sstr = 'Data/Pure_Component/{}.csv'.format(data.c[i]['name'][0])
        print 'Saving new results to {}'.format(sstr)
        save_dict_as_csv(data.c[i],sstr,Order)

    return s, p



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

    if len(data.comps) == 1:  # pure component simulation.
        # Load pure data
        data.load_pure_data() # Using data.comps
        # Find all specified outputs
        s, p = pure_sim(data, i=0)

    if len(data.comps) > 1:  # multi component simulation.
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()

        s, p = n_comp_sim(data)


    #if data.plotting:
    #    from plot import *  # allow easier func calls from python shell