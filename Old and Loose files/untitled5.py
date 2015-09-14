#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""



    #%% Trim Azeotrope and single phase points
    
    p.m['x1 Pure vapour'], p.m['y1 Pure vapour'], p.m['P Pure vapour'], \
    p.m['T Pure vapour'], p.m['x1 Pure liquid'], p.m['y1 Pure liquid'], \
    p.m['P Pure liquid'], p.m['T Pure liquid'], p.m['x1 Azeo'], p.m['y1 Azeo']\
    , p.m['P Azeo'], p.m['T Azeo'] = [],[],[],[],[],[],[],[],[],[],[],[]
    

    def trim_conditions(p, i): 
        """Defined here to simplify proceeding loop readablity"""
        if p.m['x1'][i] > 0.9999 and p.m['y1'][i] > 0.9999:
            p.m['x1 Pure vapour'].append(p.m['x1'].pop(i))
            p.m['y1 Pure vapour'].append(p.m['y1'].pop(i))
            p.m['P Pure vapour'].append(p.m['P'].pop(i))
            p.m['T Pure vapour'].append(p.m['T'].pop(i))
        if p.m['x1'][i] < 1e-5 and p.m['y1'][i] < 1e-5:
            p.m['x1 Pure liquid'].append(p.m['x1'].pop(i))
            p.m['y1 Pure liquid'].append(p.m['y1'].pop(i))
            p.m['P Pure liquid'].append(p.m['P'].pop(i))
            p.m['T Pure liquid'].append(p.m['T'].pop(i))
        else:
            p.m['x1 Azeo'].append(p.m['x1'].pop(i))
            p.m['y1 Azeo'].append(p.m['y1'].pop(i))
            p.m['P Azeo'].append(p.m['P'].pop(i))
            p.m['T Azeo'].append(p.m['T'].pop(i))
        return p

    for i in range(len(p.m['P'])):
        # Stop loop when index > len(p.m['P']) (where p.m['P'] is trimed)
        print i
        if i > len(p.m['P'])-2:
            print 'test i = {}'.format(i)
            p = trim_conditions(p, i-1)
            break
        # WHEN CONDITION IS MET, Should redo (same i index now shifted back)
        while abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
            p = trim_conditions(p, i)





    #%% Trim Azeotrope and single phase points
    
    p.m['x1 Pure vapour'], p.m['y1 Pure vapour'], p.m['P Pure vapour'], \
    p.m['T Pure vapour'], p.m['x1 Pure liquid'], p.m['y1 Pure liquid'], \
    p.m['P Pure liquid'], p.m['T Pure liquid'], p.m['x1 Azeo'], p.m['y1 Azeo']\
    , p.m['P Azeo'], p.m['T Azeo'] = ([],)*12

    def trim_conditions(p, i): 
        """Defined here to simplify proceeding loop readablity"""
        if abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
            if p.m['x1'][i] > 0.9999 and p.m['y1'][i] > 0.9999:
                p.m['x1 Pure vapour'].append(p.m['x1'].pop(i))
                p.m['y1 Pure vapour'].append(p.m['y1'].pop(i))
                p.m['P Pure vapour'].append(p.m['P'].pop(i))
                p.m['T Pure vapour'].append(p.m['T'].pop(i))
                return p
            if p.m['x1'][i] < 1e-5 and p.m['y1'][i] < 1e-5:
                p.m['x1 Pure liquid'].append(p.m['x1'].pop(i))
                p.m['y1 Pure liquid'].append(p.m['y1'].pop(i))
                p.m['P Pure liquid'].append(p.m['P'].pop(i))
                p.m['T Pure liquid'].append(p.m['T'].pop(i))
                return p
            else:
                p.m['x1 Azeo'].append(p.m['x1'].pop(i))
                p.m['y1 Azeo'].append(p.m['y1'].pop(i))
                p.m['P Azeo'].append(p.m['P'].pop(i))
                p.m['T Azeo'].append(p.m['T'].pop(i))
                return p

    for i in range(len(p.m['P'])):
        # Stop loop when index > len(p.m['P']) (where p.m['P'] is trimed)
        if i > len(p.m['P'])-2:
            if abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
                p = trim_conditions(p, i-1)
                
            break
        # WHEN CONDITION IS MET, Should redo (same i index now shifted back)
        while abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
            p = trim_conditions(p, i)