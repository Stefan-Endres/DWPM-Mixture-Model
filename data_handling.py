#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Data read/write.

ex.
>>> dataset1.load_VLE(['acetone','water'])
>>> dataset2.load_VLE(['benzene','cyclohexane'])
"""
import ConfigParser
import os

#%% Load data classes
class ImportData:
    """Load data of coupounds specified with Data.load_pure_data()"""
    def __init__(self):
        self.c = []    # creates a new empty list for components

        # Read configuration file
        # TODO: Add argparse optionals to read from.
        config = ConfigParser.ConfigParser()
        configfile = 'config.cfg'

        if os.path.exists(configfile):
            config.read('config.cfg')
        else:
            message = ("Cannot find config file {0}. "
                       "Try copying sample_config.cfg to {0}.").format(
                                                                    configfile)
            raise EnvironmentError(message)

        # Load data vars into class and parse lists
        # TODO: There's probably a simpler way simpler way to do this l 33-79.
        compstr = config.get('inputs', 'compounds')
        self.comps = compstr.split('-')
        if isinstance(self.comps, basestring): #  1x1 list TODO: Not needed?
            self.comps = [self.comps]

        self.EOS = config.get('inputs', 'EOS')
        phasestr = config.get('inputs', 'valid_phases')
        self.phases = phasestr.split('-')
        self.model = config.get('inputs', 'model')
        self.r = config.get('inputs', 'r')
        if self.r == 'None': # TODO: Get None vals from empty config entry?
            self.r = None

        self.s = config.get('inputs', 's')
        if self.s == 'None':
            self.s = None

        self.T = config.get('inputs', 'T')
        if self.T == 'None':
            self.T = None

        self.P = config.get('inputs', 'P')
        if self.P == 'None':
            self.P = None

        self.Z_0 = config.get('inputs', 'Z_0')
        if self.Z_0 == 'None':
            self.Z_0 = None
        else:
            self.Z_0 = self.Z_0.split('-')

        self.save_results = config.getboolean('output settings',
                                              'save_results')
        self.force_update = config.getboolean('output settings',
                                               'force_update')
        self.plot_iso = config.getboolean('output settings',
                                          'plot_iso')
        self.plot_gibbs = config.getboolean('output settings',
                                            'plot_gibbs')
        self.plot_pure = config.getboolean('output settings',
                                           'plot_pure')
        self.save_pure = config.getboolean('output settings',
                                           'save_pure')
        self.force_pure_update = config.getboolean('output settings',
                                                   'force_pure_update')
        self.datadir =  config.get('paths',
                                   'datadir')

    def load_pure_data(self):
        from csvDict import load_csv_as_dict
        """
        Returns the pure component VLE data for specified components.
        """
        for i, comp in enumerate(self.comps):
            ldstr = os.path.join(self.datadir
                                 ,'Pure_Component','{}.csv'.format(comp))
            print ldstr
            #ldstr           = '/Pure_Component/{}.csv'.format(comp)
            try:
                Data = load_csv_as_dict(ldstr)
                Data['name'] = [comp,]
                Data['model'] = self.model
                self.c.append(Data)

            except IOError:  # Raise error if not found
                raise IOError('Data for specified component '  \
                               +'"{}" not found'.format(comp))

    def load(self):
        from csvDict import load_csv_as_dict
        # Find file name path for specified components
        filename = '_'.join(self.comps)
        ldstr = os.path.join(self.datadir, 'nComp_E', filename + '.csv')
        print 'TEST'
        try: # Load data from file path
            Data  = load_csv_as_dict(ldstr)
            self.VLE = Data
        except IOError: # Raise error if not found
            raise IOError('Phase data for '
                          'system "{}" not found'.format(filename))


#    def load_E(self):
#        """
#        Returns multicomponent equilibrium data from specified components.
#        """
#        self.load('Data/nComp_E/', self.comps)

    def test_internal(self): # TEST; DELETE
         self.test_int() # TEST; DELETE
         #print

    def test_int(self): # TEST; DELETE
         print 'Test Succesful!'


def parameter_build(Data):
    """
    Move data container to parameter output dictionary and find the
    critical Van der Waals contants if not defined.

    Parameters
    ----------
    Data : Dictionary containing data loaded from the stored .csv file.

    comp : string of definin

    """

    if len(Data): #TODO: Merge single and multi comp param defs here
        pass
    p = {'T'   : Data['T (K)'],
         'P'   : Data['P (Pa)'],
         'T_c' : Data['T_c (K)'][0],
         'P_c' : Data['P_c (Pa)'][0],
         'V_c' : Data['V_c (m3 mol-1)'][0],
         'Z_c' : Data['Z_c'][0],
         'R'   : Data['R (m3 Pa K-1 mol-1)'][0],
         'w'   : Data['w'][0],
         'a_c' : Data['a_c (Pa m6 mol-2)'][0],
         'b_c' : Data['b_c (m3 mol-1)'][0],
         'vT'  : Data['virialT'],
         'vB'  : Data['virialB'],
         'Model': Data['model'],
         'name': Data['name']
         }

    if Data['model'] == 'Adachi-Lu': # Find model params if not defined
        p['m'] = Data['m (Adachi-Lu)'][0]
    elif Data['model'] == 'Soave':
        p['m'] = Data['m (Soave)'][0]

    if p['a_c'] == '' or p['b_c'] == '':
        p['b_c'] = p['R']*p['T_c']/(8*p['P_c'])
        p['a_c'] = 27*(p['R']**2)*(p['T_c']**2)/(64.0*p['P_c'])
    else:
        pass

    for key, value in p.iteritems(): # Filter out '' values
        if not value.__class__ == float:
            p[key] = filter(lambda a: a != '', value)

    return p

class MixParameters:
    """
    Store mixture and pure parameters in the same class.

    Parameters
    ----------
    Data : Dictionary containing data loaded from the stored .csv file.

    """
    def __init__(self):
        self.c = []  # creates a new empty list for components
        self.c.append('nan')  # Define an empty set in index 0

    def mixture_parameters(self, data_VLE, data):
        """Mixture model parameters"""
        M = {'T'      : data_VLE['T (K)'],  # Temperature Pressure data
             'P'      : data_VLE['P (Pa)'],
             'n'      : len(data.comps),
             'phases' : len(data.phases),
             'Model'  : data.EOS,
             'Valid phases' : data.phases
             }

        # Define phases
        for i in range(len(data.phases)):
            M[data.phases[i]] = ['nan']

        # Define equilibria for each component in phase
        for j in range(1, M['n'] + 1):
            for i in range(len(data.phases)):
                M[data.phases[i]].append(data_VLE[data.phases[i] +
                                               '{}'.format(j)])

                # NOTE: This routine will change component strings the for the
                # equilibrium of each phase into a list simple list for each
                # phase ex. data_VLE['x1'], data_VLE['x2'] becomes M['x'][1], M['x'][2]

        # Define model paramters
        # Empty lists for model interaction paramters
        M['k'] = []
        [M['k'].append(['nan']) for j in range(M['n'] + 1)]

        if data.EOS == 'DWPM':
            # Find the interaction paramters between and put them into
            # component lists (ex. data_VLE['k12']  --> M['k'][1][2])

            for j in range(1, M['n'] + 1):
                for i in range(1, M['n'] + 1):
                    # Define empty list
                    M['k'][j].append('nan')
                    if i != j:  # Define interaction paramter
                        M['k'][j][i] = data_VLE['k{J}{I}'.format(J=j, I=i)][0]

            if data.r is None:
                M['r'] = data_VLE['r']
            else:
                M['r'] = data.r

            if data.s is None:
                M['s'] = data_VLE['s']
            else:
                M['s'] = data.s


        for key, value in M.iteritems():  # Filter out '' values
            if not value.__class__ == float:
            #  if key != 'x' and key != 'y' and key != 'k' and key != 'phases':
                try:
                    M[key] = filter(lambda a: a != '', value)
                except(TypeError):
                    pass

        self.m = M

    def parameters(self, data):
        # (Recieving data[i] dict, not full data class)
        p = parameter_build(data) # should be data_handling.parameter_build?
        #p['name'] = Data['name']
        self.c.append(p)


if __name__ == '__main__':
    data = ImportData()
    if len(data.comps) > 1:
        data.load()
    elif len(data.comps) == 1:
        data.load_pure_data()

    data.test_internal()
    #data = ImportData()
    #data.load_pure_data(Compounds) # Call ex. component 1: data.c[0]
    #if len(Compounds) > 1:
    #    data.load_VLE(Compounds)       # Call ex. component 1: data.c[0]
