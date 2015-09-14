#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""


#%% Error func (over all data points)
def error(k12, s, p):
    """
    Returns the error function at specified paramters.
    
    TO DO: ADD PUNISHMENT FOR MATH ERRORS (Every for loop?)
    
    (Previously insdie optim_mixture_parameter(s, p) )
    """
    # Init error over data set
    #s = args[0]
    #p = args[1]
    #print p.m['P'][0]
    if p.m['Model'] == 'VdW standard':
        p.m['k12'] = k12
        
    err = 0.0

    #% Errors related to pure phases and Azeotropes
    # TO DO; Are Pure phase checks needed since mix params do not depend
    # on pure VdW behavior?

#        for i in range(len(p.m['x1 Pure vapour'])): # Vapour Pres. penalty
#        # TO DO; Use del_g_v_min_l0(s,p) g_mix min ?? 
#        # INSPECT THEORETICAL VIABILITY
#            # Update values
#            s.s['P'], s.s['T'] = p.m['P Pure vapour'][i], \
#                                 p.m['T Pure vapour'][i]
#            s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
#            s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
#            s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
#            s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
#            
#            g_v_l = del_g_v_min_l0(s.c[1],p.c[1])
#                       
#            del_g_mix(s, p, x_1=p.m['x1 Pure vapour'][i], update_pure=True)
#            if 
    #% Azeotropes penatly    
    for i in range(len(p.m['P Azeo'])):
        # Update values
        s.s['P'], s.s['T'] = p.m['P Azeo'][i], \
                             p.m['T Azeo'][i]
        s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
        s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
        # Find difference in Gibbs energies between liq and vapour
        dg_v_l = del_g_mix(s, p, x_1=p.m['y1 Azeo'][i], \
        update_pure=True).m['del g_mix'] \
        - del_g_mix(s, p, x_1=p.m['x1 Azeo'][i], \
        update_pure=True).m['del g_mix'] 
        if dg_v_l > 1e-6:
            err += 10**dg_v_l
  
#            x1_diff = abs(p.m['x1'][i] - p.m['y1'][i]) 
#            if x1_diff < 1e-6 and True:
#                err += 10.0**x1_diff
#            #% Wrong pure phase penalty
#            if abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-6: 
#                err += 10.0**abs(p.m['x1'][i] - p.m['y1'][i])
    #% Find phase split errors at current parameters  
    for i in range(len(p.m['x1'])):
        # Uptate T, P INCLUDE THIS IN RELEVANT FUNCTION
        s.s['P'], s.s['T'] = p.m['P'][i], p.m['T'][i]
        s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
        s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3

        s, p = phase_split(s, p, est = [min(p.m['x1'][i], p.m['y1'][i]), \
                                        max(p.m['x1'][i], p.m['y1'][i])])

#            if i == 103:
#                print 'x = {}, y = {}'.format(s.m['x1'],s.m['y1'])
                                        
        err+= (p.m['x1'][i] - s.m['x1'])**2 + (p.m['y1'][i] - s.m['y1'])**2
        
    #print 'done' # DEBUGGING REMOVE
    return err

#%% Trim Azeotrope and single phase points    
def trim_azeo_and_pure(s, p):
    """Previously insdie optim_mixture_parameter(s, p) """
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
        if i > len(p.m['P'])-2:
            if abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
                p = trim_conditions(p, i-1)
                
            break
        # WHEN CONDITION IS MET, Should redo (same i index now shifted back)
        while abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
            p = trim_conditions(p, i)  

    return s, p
 
#%%    
def optim_mixture_parameter(s, p):
    """
    Stuff happens here, parameters come out.
    """
    from scipy.optimize import fmin_l_bfgs_b, anneal
    from numpy import array
    
    # b11 and b22 taken as the critical parameters
    s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
    
    #%% trim azeotrope and single comp. data points
    s, p = trim_azeo_and_pure(s, p)
    #%%    
    s.m['err test'] = error(0.0135,s,p)
    #error(0.011,s,p)

    #k12_0 = array([0.0135]) # Initial guess
    #Bounds = [(0,1.0)]
    #mininum = fmin_l_bfgs_b(error, k12_0, args=(s,p), bounds=Bounds,
    #                        approx_grad=True)

    #%% TEST PLOT DELETE                      
    from numpy import linspace

    k12r = linspace(0.01,0.015, num=1e2)
    k12r = linspace(0.0001,0.01, num=1e2)
    s.m['errstore'] = []
    for k12 in k12r:
        s.m['errstore'].append(error(k12,s,p))
        print len(s.m['errstore'])

    #%% TEST PLOT DELETE
    from matplotlib import rc
    from matplotlib import pyplot as plot
    from numpy import array
    plot.figure()
    plot.plot(k12r, s.m['errstore'], 'k', label='Error function')
    plot.xlabel(r"$k_{12}$", fontsize=14)
    plot.ylabel(r"$\epsilon$    ", fontsize=14, rotation='horizontal')
#%%    
    
    #s.m['minimum'] = minimum
    #anneal
    return s, p

def error_func(s,p):
    """Error function meshgrid plot for DWMP"""
    import numpy as np
    from scipy import optimize
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
   
    # b11 and b22 taken as the critical parameters
    s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
    # Trim raw data points
    s, p = trim_azeo_and_pure(s, p)
    p.m['k12'] = 0.0 # TO DO
    p.m['k21'] = 0.0 # TO DO
    # Get meshgrid for error func range
    rrange = np.linspace(-1, 1.5)
    srange = np.linspace(-1, 1.5)
    rg, sg = numpy.meshgrid(rrange, srange)
    zerr = numpy.zeros((50,50))

    # Find error values
    for i in range(rg.shape[0]):
        for j in range(sg.shape[0]):
            # Do stuff to find "z"
            p.m['r'] = rg[i, j]
            p.m['s'] = sg[i, j]
            print 'Progress = {} %'.format(((i+1)/50.0)*100.0)
            try:
                zer = error(0.0,s,p)
            except:
                zer = 100.0
            zerr[i,j] = zer
            
    # Plot Meshgrid        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(rg, sg, zerr, rstride=1, cstride=1,
                           cmap=plt.cm.jet, linewidth=0, antialiased=False)
    
    ax.set_xlabel('r')
    ax.set_ylabel('s')
    ax.set_zlabel('$\epsilon$',rotation='vertical')
    #ax.set_title('Six-hump Camelback function')