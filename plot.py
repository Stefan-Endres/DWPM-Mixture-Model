#!/usr/bin/env python
# (TODO)
from models import van_der_waals
VdW = van_der_waals.VdW()

#%% Plot pure Functions
def plot_Psat(s, p, options, figno=None):
    """

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.

    p : class
        Contains the dictionary describing the mixture parameters.

    options : dict
              Options to pass to matplotlib.pyplot.plot

    """
    import matplotlib.pyplot as plot
    from numpy import linspace, interp

    #TODO: Change these s dict values to normal local variables and test

    s['T_sat store'] = linspace(p['T'][0], p['T'][len(p['T'])-1])
    s['P_sat store'] = []
    s['P_est'] = interp(s['T_sat store'], p['T'], p['P'])
    i = 0
    for T, P in zip(s['T_sat store'][:len(s['T_sat store'])-1],
                    s['P_est'][:len(s['T_sat store'])-1]): # Trim crit.
        s['T'] = T # Solve P_sat at this Temperature
        s['P'] = P # P est
        s = VdW.Psat_V_roots(s,p,tol=1e-1)
        s['P_sat store'].append(s['P_sat'])

        s['P_sat store'].append(p['P_c']) # Append Critical point



    plot.figure(figno)
    plot.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    plot.rcParams.update(options)
    p['P'] = [Pa for Pa in p['P']] # Pa -> kPa
    s['P_sat store'] = [Pa for Pa in s['P_sat store']] # Pa -> kPa
    plot.plot(p['T'], p['P'], 'xr', label='Data points')
    plot.plot(s['T_sat store'],s['P_sat store'], '--r',
              label='Van der Waals EoS %s m = %s'% (p['Model'], p['m']))
    plot.xlabel("Temperature / K")
    plot.ylabel("Pressure$^{sat}$ / Pa")
    plot.title("Van der Waals EoS correlation for $%s$" \
                % (p.c[0]['name'][0]))
    plot.legend(loc=options['legend.loc'])
    return


# %% nComp Plots
class GibbsSurface: #TODO add g_mix plotes here and refactor local variables
    def __init__(self):
        pass

def plot_g_mix(s, p, g_x_func, Tie=None, x_r=1000, FigNo = None):
    """
    Plots the surface of the Gibbs energy of mixing function. Binary and
    Trenary plots only.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['m']

    Tie : list of vectors, optional
          Equilibrium tie lines (of n - 1 independent components) to be added
          to the plots.

          For a binary system specify 2 tie lines as
          [[x_1_a, x_1_b], [x_1_c, x_1_d]]

          For a trenary system specify the parameters in the plane construction
          as
          [[G_p, x_1, lambda_1, x_2, lambda_2]]
          where G_p is the solution to the global problem at [x_1, x_2] and
          lamda_i is the solution duality multipliers.



    x_r : int, optional
          Number of composition points to plot.

    FigNo : int or None, optional
            Figure number to plot in matplotlib. Specify None when plotting
            several figures in a loop.

    Dependencies
    ------------
    numpy
    matplotlib

    """

    import numpy as np

    #% Binary
    if p.m['n'] == 2:
        from matplotlib import pyplot as plot
        #% Initialize
        s.m['g_mix range'] = {} # Contains solutions for all phases
        s.m['g_mix range']['t'] = []
        for ph in p.m['Valid phases']:
            s.m['g_mix range'][ph] = []

        s.c[1]['x_range'] = np.linspace(0, 1, x_r)
        for i in range(len(s.c[1]['x_range']-1)):  # Calculate range of G_mix
            X = [s.c[1]['x_range'][i]]
            s = s.update_state(s, p, X = X, Force_Update=True)
            s = g_x_func(s, p)
            s.m['g_mix range']['t'].append(np.float64(s.m['g_mix']['t']))
            for ph in p.m['Valid phases']:
                s.m['g_mix range'][ph].append(np.float64(s.m['g_mix'][ph]))

        if FigNo is None:
            plot.figure()
        else:
            plot.figure(FigNo)

        for ph in p.m['Valid phases']:
            plot.plot(s.c[1]['x_range'], s.m['g_mix range'][ph], label=ph)

        plot.plot(s.c[1]['x_range'], s.m['g_mix range']['t'])
        if Tie is not None: # Add tie lines
            for point in Tie:
                X = [point[0]]
                s = s.update_state(s, p, X = X, Force_Update=True)
                s = g_x_func(s, p)
                G1 = np.float64(s.m['g_mix']['t'])
                X = [point[1]]
                s = s.update_state(s, p, X = X, Force_Update=True)
                s = g_x_func(s, p)
                G2 = np.float64(s.m['g_mix']['t'])
                plot.plot(point,[G1, G2])
                Slope = (G1 - G2)/ (point[0] - point[1])
                plot.annotate('slope = {}'.format(Slope),
                              xy=(point[1], G2 + G2*0.5) )

        plot.xlabel(r"$z_1$", fontsize=14)
        plot.ylabel(r"$\Delta$g", fontsize=14)
        plot.title("{}-{} Gibbs energy of mixing at T = {}, P = {}".format(
                    p.c[1]['name'][0],
                    p.c[2]['name'][0],
                    s.m['T'],
                    s.m['P']))
        plot.legend()
        plot.show()
        return

    #% Trenary
    if p.m['n'] == 3:
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plot
        from matplotlib import cm

        # Define Tie planes
        def tie(X, T):
            G_p, x_1, lambda_1, x_2, lambda_2 = T
            return G_p + lambda_1 * (X[0] - x_1) + lambda_2 * (X[1] - x_2)

        x_range = np.linspace(1e-15, 1.0, x_r)
        y_range = np.linspace(1e-15, 1.0, x_r)
        xg, yg = np.meshgrid(x_range, y_range)
        s.m['g_mix range'] = {}
        s.m['g_mix range']['t'] = np.zeros((x_r, x_r))
        for ph in p.m['Valid phases']:
                    s.m['g_mix range'][ph] = np.zeros((x_r, x_r))

        Tie_planes = []
        for z in range(len(Tie)):
            Tie_planes.append(np.zeros((x_r, x_r)))

        for i in range(xg.shape[0]):
            for j in range(yg.shape[0]):
                 X = [xg[i, j], yg[i, j]]  # [x_1, x_2]
                 if sum(X) > 1.0:#1.0:
                     for z in range(len(Tie)):
                         Tie_planes[z][i, j] = None

                     s.m['g_mix range']['t'][i, j] = None
                     for ph in p.m['Valid phases']:
                         s.m['g_mix range'][ph][i, j] = None

                 else:
                     # Tie lines
                     for z in range(len(Tie)):
                         Tie_planes[z][i, j] = tie(X, Tie[z])

                     # g_func
                     s = s.update_state(s, p, X = X, Force_Update=True)
                     s = g_x_func(s, p)
                     s.m['g_mix range']['t'][i, j] = np.float64(
                                                             s.m['g_mix']['t'])
                     for ph in p.m['Valid phases']:
                         s.m['g_mix range'][ph][i, j] = np.float64(
                                                             s.m['g_mix'][ph])

        if FigNo is None:
            fig = plot.figure()# plot.figure()
        else:
            fig = plot.figure(FigNo)


        # Plots
        ax = fig.gca(projection='3d')
        X, Y = xg, yg

        # Gibbs phase surfaces
        rst = 1
        cst = 1
        for ph in p.m['Valid phases']:
            Z = s.m['g_mix range'][ph]
            ax.plot_surface(X, Y, Z, rstride=rst, cstride=cst, alpha=0.1)
            cset = ax.contour(X, Y, Z, zdir='x', offset=-0.1, cmap=cm.coolwarm)
            cset = ax.contour(X, Y, Z, zdir='y', offset=-0.1, cmap=cm.coolwarm)

        # Gibbs minimum surface
        Z= s.m['g_mix range']['t']
        ax.plot_surface(X, Y, Z, rstride=rst, cstride=cst, alpha=0.1,color='r')
        if True:
            cset = ax.contour(X, Y, Z, zdir='x', offset=-0.1, cmap=cm.coolwarm)
            cset = ax.contour(X, Y, Z, zdir='y', offset=-0.1, cmap=cm.coolwarm)

        # Planes
        rst = 4
        cst = 4
        for z in range(len(Tie)):
            ax.plot_surface(X, Y, Tie_planes[z], rstride=rst, cstride=cst,
                        alpha=0.1)

        ax.set_xlabel('$x_1$')
        ax.set_xlim(0, 1)
        ax.set_ylabel('$x_2$')
        ax.set_ylim(0, 1)
        ax.set_zlabel('$\Delta g$', rotation = 90)
        plot.show()
        return s

    else:
        raise IOError('Too many independent components to plot hypersurface'
                      +'in R^3.')

class Iso:
    def __init__(self, T=None, P=None):
        self.T = T
        self.P = P
        pass

    def plot_iso(self, s, p, g_x_func, T=None, P=None, res=30, n=100,
                 tol=1e-9, gtol=1e-2, n_dual=100, phase_tol=1e-3):
        import numpy
        from ncomp import equilibrium_range as er
        P_data_f = numpy.array(p.m['P'])
        T_data_f = numpy.array(p.m['T'])

        if T is not None:
            for t in T:
                iso_ind = numpy.where(T_data_f == t)
                P_data = P_data_f[iso_ind]
                #T_data = T_data_f[iso_ind]
                PT_Range = [(min(P_data), max(P_data)),
                            (t, t)]

                P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph = er(g_x_func,
                                         s, p, PT_Range=PT_Range, n=n, res=res,
                                         tol=tol, gtol=gtol, n_dual=n_dual,
                                         phase_tol=phase_tol)

                # (Process results)
                data_x = {'x': [0,0]}
                # Plot resulting isotherm
                #self.plot_iso_t_bin(t, P_data,



    def plot_iso_t_bin(self, T, data_p, data_x, p, model_p=None, model_x=None
                       , k=['All'], FigNo=None, plot_options=None):
        """
        Plot binary isotherms for the specified data and model ranges

        Parameters
        ----------

        T : float
            Temperature of isotherm to plot.

        data_p : vector
                 Pressure data values at each point

        data_x : dict containing vectors
                 Contains the composition data points at each data_p for
                 every valid phase

                 ex. data_x = {'x': [p.m['x'][1][25:36],  # x_1
                                     p.m['x'][2][25:36]], # x_2

                               'y': [p.m['y'][1][25:36],  # y_1
                                     p.m['y'][2][25:36]]  # y_2
                               }

        p : class
            Contains the dictionary describing the parameters.

        data_p : vector, optional
                 Pressure model values at each point

        data_x : dict containing vectors, optional
                 Contains the simulated composition points at each data_p for
                 every valid phase, specified similarly to data_x.

        k : list, optional
            List contain valid phases for the current equilibrium calculation.
            ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.


        FigNo : int or None, optional
                Figure number to plot in matplotlib. Specify None when plotting
                several figures in a loop.

        options : dict
                  Options to pass to matplotlib.pyplot.plot

        """
        if k == ['All']:
            k = p.m['Valid phases']

        from matplotlib import pyplot as plot
        if FigNo is None:
            plot.figure()
        else:
            plot.figure(FigNo)

        # Plot data points:
        for ph in k:
            plot.plot(data_x[ph][1], data_p, 'x', label='{} data'.format(ph))

        # Plot model points
        if model_p is not None:
            for ph in p.m['Valid phases']:
                plot.plot(model_x[ph][1], data_p, '-',
                          label='{} model'.format(ph))

        plot.xlabel(r"$z_1$", fontsize=14)
        plot.ylabel(r"P (Pa)", fontsize=14)
        plot.title("{}-{} isotherm at T = {}".format(p.c[1]['name'][0],
                                                     p.c[2]['name'][0],
                                                     T))
        plot.legend()
        plot.show()
        return

    def plot_iso_p_bin(self):
        """Plot binary isobars for the specified data and model ranges"""
        pass

def plot_ep(func, X_r, s, p, args=()):
    """
    Plot the speficied single var input error function over a range X_r
    """
    from matplotlib import pyplot as plot
    #% Binary
    if p.m['n'] == 2:
        ep = []
        for X in X_r:
            s.update_state(s, p,  X = X, phase = ['All'], Force_Update=True)
            try:
                ep.append(func(X, *args)[0])
            except(IndexError):
                ep.append(func(X, *args)) # Scalar outputs
            #ep.append(func(X, *args))

        plot.figure()
        plot.plot(X_r, ep)
        plot.xlabel(r"$x_1$", fontsize=14)
        plot.ylabel(r"$\epsilon$", fontsize=16)
        plot.title("{}-{} at T = {}, P = {}".format(
                    p.c[1]['name'][0],
                    p.c[2]['name'][0],
                    s.m['T'],
                    s.m['P']))
        plot.show()
    return


if __name__ == '__main__':
     plot_options = {'text.usetex' : True, # Options for all plots
                     'font.size' : 12,
                     'font.family' : 'lmodern',
                     'text.latex.unicode': True
                     }