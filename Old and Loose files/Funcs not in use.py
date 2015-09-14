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

 
#%%    
def optim_mixture_parameter(s, p):
    """
    Stuff happens here, parameters come out.
    """
    from scipy.optimize import fmin_l_bfgs_b, anneal
    from numpy import array
    
    # b11 and b22 taken as the critical parameters
    s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
    
    #%%    
    #s.m['err test'] = error(0.0135,s,p)
    #error(0.011,s,p)

    #k12_0 = array([0.0135]) # Initial guess
    #Bounds = [(0,1.0)]
    #mininum = fmin_l_bfgs_b(error, k12_0, args=(s,p), bounds=Bounds,
    #                        approx_grad=True)

    #%% TEST PLOT DELETE                      
    from numpy import linspace

    k12r = linspace(0.01,0.015, num=1e2)
    k12r = linspace(0.0001,0.01, num=1e2)
    #k12r = linspace(0.0,0.2, num=1e3)
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
    """ OLD Error function meshgrid plot for DWMP"""
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