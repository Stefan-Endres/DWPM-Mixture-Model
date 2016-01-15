#!/usr/bin/env python
# -*- coding: utf-8 -*-

def parameters(Data,I):
    """
    Move data container to parameter output dictionary and find the critical
    Van der Waals contants if not defined.

    Parameters
    ----------
    Data : Dictionary containing data loaded from the stored .csv file.

    I : Input dictionary, must contain key 'Model' defining the a parameter
        dependancy model.
    """
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
         'Model': I['Model']
         }

    if I['Model'] == 'Adachi-Lu': # Find model params if not defined
        p['m'] = Data['m (Adachi-Lu)'][0]
    elif I['Model'] == 'Soave':
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
