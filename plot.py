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
        g_mix_r = {} # Contains range of solutions for all phases
        g_mix_r['t'] = []
        for ph in p.m['Valid phases']:
            g_mix_r[ph] = []

        s.c[1]['x_range'] = np.linspace(0, 1, x_r)
        for i in range(len(s.c[1]['x_range']-1)):  # Calculate range of G_mix
            X = [s.c[1]['x_range'][i]]
            s = s.update_state(s, p, X = X, Force_Update=True)
            s = g_x_func(s, p)
            g_mix_r['t'].append(np.float64(s.m['g_mix']['t']))
            for ph in p.m['Valid phases']:
                g_mix_r[ph].append(np.float64(s.m['g_mix'][ph]))

        if FigNo is None:
            plot.figure()
        else:
            plot.figure(FigNo)

        for ph in p.m['Valid phases']:
            plot.plot(s.c[1]['x_range'], g_mix_r[ph], label=ph)

        plot.plot(s.c[1]['x_range'], g_mix_r['t'])
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
        #plot.show()
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
        g_mix_r = {}
        g_mix_r['t'] = np.zeros((x_r, x_r))
        for ph in p.m['Valid phases']:
            g_mix_r[ph] = np.zeros((x_r, x_r))

        Tie_planes = []
        for z in range(len(Tie)):
            Tie_planes.append(np.zeros((x_r, x_r)))

        for i in range(xg.shape[0]):
            for j in range(yg.shape[0]):
                 X = [xg[i, j], yg[i, j]]  # [x_1, x_2]
                 if sum(X) > 1.0:#1.0:
                     for z in range(len(Tie)):
                         Tie_planes[z][i, j] = None

                         g_mix_r['t'][i, j] = None
                     for ph in p.m['Valid phases']:
                         g_mix_r[ph][i, j] = None

                 else:
                     # Tie lines
                     for z in range(len(Tie)):
                         Tie_planes[z][i, j] = tie(X, Tie[z])

                     # g_func
                     s = s.update_state(s, p, X = X, Force_Update=True)
                     s = g_x_func(s, p)
                     g_mix_r['t'][i, j] = np.float64(s.m['g_mix']['t'])
                     for ph in p.m['Valid phases']:
                         g_mix_r[ph][i, j] = np.float64(s.m['g_mix'][ph])

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
            Z = g_mix_r[ph]
            ax.plot_surface(X, Y, Z, rstride=rst, cstride=cst, alpha=0.1)
            #cset = ax.contour(X, Y, Z, zdir='x', offset=-0.1, cmap=cm.coolwarm)
            #cset = ax.contour(X, Y, Z, zdir='y', offset=-0.1, cmap=cm.coolwarm)

        # Gibbs minimum surface
        Z = g_mix_r['t']
        ax.plot_surface(X, Y, Z, rstride=rst, cstride=cst, alpha=0.1,color='r')
        if False:
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
        #plot.show()
        return s

    else:
        import logging
        logging.warn('Too many independent components to plot hypersurface '
                      +'in R^3.')

class Iso:
    def __init__(self, T=None, P=None):
        self.T = T
        self.P = P
        pass

    def plot_iso(self, s, p, g_x_func, T=None, P=None, res=30, n=1000,
                 tol=1e-9, gtol=1e-2, n_dual=300, phase_tol=1e-3,
                 LLE_only=False, VLE_only=False, Plot_Results=False,
                 data_only=False):
        """
        Main function to call when for plotting isotherms/isobars for either
        binary or ternary systems.

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
                   Returns a class containing scalar value .m['g_mix']['t']

        T : list
            Isotherms to simulate

        P: list
            Isobars to simulate

        res : integer
              Specifies the number of data points to be simulated within the
              specified range.

        tol : scalar, optional
              Tolerance used in ``dual_equal``, if epsilon >= UBD - LBD that will
              terminate the routine.

        gtol : scalar, optional
              Minimum tolerance between hyperplane solution
              Note: The Dual solution is not perfect so a low tolerance is
              required, but a too low tol could potentially include points that do
              not truly lie on the equilibrium plane within the considered
              instability region.

        n_dual : scalar, optional
                Number of sampling points used in the tgo routine in solving LBD
                of the dual problem.
                Note: It is recommended to use at least ``100 + p.m['n'] * 100``

        phase_tol : scalar, optional
                    The minimum seperation between equilibrium planes required to
                    be considered a phase. Defaults to 0.001

        LLE_only : boolean, optional
                   If True then only phase seperation of same volume root
                   instability will be calculated.

        VLE_only : boolean, optional
                   If True then phase seperation of same volume root instability
                   will be ignored.

        Plot_Results : boolean, optional
                       If True the g_mix curve with tie lines will be plotted for
                       binary and ternary systems.

        data_only : boolean, optional
                    If True only data will be plotted with no model simulations
        """

        if T is not None:
            for t in T:
                # Find model results and data points
                (P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph, data_x_mph,
                 data_x_ph, data_t, data_p) =  self.iso_range(s, p,
                                   g_x_func, T=t, P=None, res=res,
                                   n=n, tol=tol, gtol=gtol, n_dual=n_dual,
                                   phase_tol=phase_tol, LLE_only=LLE_only,
                                   VLE_only=VLE_only, Plot_Results=True,
                                   data_only=data_only)

                # Process VLE points
                plot_kwargs = {}
                if not LLE_only:
                    if not data_only:
                        model_x_mph, model_p_mph, model_t_mph = \
                            self.process_VLE_range(p, P_range, T_range,
                                                   r_mph_eq, r_mph_ph)
                    else:
                        (model_x_mph, model_p_mph, model_t_mph) = (None, None,
                                                                   None)

                    plot_kwargs['model_p_mph'] = model_p_mph
                    plot_kwargs['model_x_mph'] = model_x_mph

                # Process LLE points
                if not VLE_only:
                    if not data_only:
                        model_x_ph, model_p_ph, model_t_ph = \
                            self.process_LLE_range(p, P_range, T_range,
                                                   r_ph_eq)
                    else:
                        (model_x_ph, model_p_ph, model_t_ph) = (None, None,
                                                                   None)
                    plot_kwargs['model_p_ph'] = model_p_ph
                    plot_kwargs['model_x_ph'] = model_x_ph

                # Plot each isotherm
                if p.m['n'] == 2:
                    self.plot_iso_t_bin(t, p,
                                        data_p=data_p,
                                        data_x_mph=data_x_mph,
                                        data_x_ph = data_x_ph,
                                    #    model_p_mph=model_p_mph,
                                        #model_t_mph=model_t_ph,
                                   #     model_x_mph=model_x_mph,
                                    #    model_p_ph=model_p_ph,
                                        #model_t_ph=model_t_ph,
                                   #     model_x_ph=model_x_ph,
                                        LLE_only=LLE_only,
                                        VLE_only=VLE_only,
                                        **plot_kwargs)
                elif p.m['n'] == 3:
                    pass
                else:
                    import logging
                    logging.warn('Dimensionality too high, ignoring plot'
                                 'request')

        if P is not None:
            for Pre in P:
                # Find model results and data points
                (P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph,
                 data_x_mph,
                 data_x_ph, data_t, data_p) = self.iso_range(s, p,
                                                           g_x_func,
                                                           T=None,
                                                           P=Pre,
                                                           res=res,
                                                           n=n,
                                                           tol=tol,
                                                           gtol=gtol,
                                                           n_dual=n_dual,
                                                           phase_tol=phase_tol,
                                                           LLE_only=LLE_only,
                                                           VLE_only=VLE_only,
                                                           Plot_Results=True,
                                                           data_only=data_only)

                # Process VLE points
                plot_kwargs = {}
                if not LLE_only:
                    model_x_mph, model_p_mph, model_t_mph = \
                        self.process_VLE_range(p, P_range, T_range,
                                               r_mph_eq, r_mph_ph)

                plot_kwargs['model_t_mph'] = model_t_mph
                plot_kwargs['model_x_mph'] = model_x_mph

                # Process LLE points
                if not VLE_only:
                    # Process results
                    model_x_ph, model_p_ph, model_t_ph = \
                        self.process_LLE_range(p, P_range, T_range,
                                               r_ph_eq)

                    plot_kwargs['model_t_ph'] = model_t_ph
                    plot_kwargs['model_x_ph'] = model_x_ph

                # Plot each isotherm
                if p.m['n'] == 2:
                    self.plot_iso_p_bin(Pre, p,
                                        data_p=data_p,
                                        data_x_mph=data_x_mph,
                                        data_x_ph=data_x_ph,
                                        # model_p_mph=model_p_mph,
                                   #     model_t_mph=model_t_ph,
                                   #     model_x_mph=model_x_mph,
                                        # model_p_ph=model_p_ph,
                                   #     model_t_ph=model_t_ph,
                                   #     model_x_ph=model_x_ph,
                                        LLE_only=LLE_only,
                                        VLE_only=VLE_only,
                                        **plot_kwargs)
                elif p.m['n'] == 3:
                    pass
                else:
                    import logging
                    logging.warn(
                        'Dimensionality too high, ignoring plot'
                        'request')
        return



    def iso_range(self, s, p, g_x_func, T=None, P=None, res=30, n=1000,
                  tol=1e-9, gtol=1e-2, n_dual=300, phase_tol=1e-3,
                  LLE_only=False, VLE_only=False, Plot_Results=False,
                  data_only=False):
        """
        Function used to find model ranges over isotherms/bars and organize
        the results into data containers that can be used plot functions.

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
                   Returns a class containing scalar value .m['g_mix']['t']

        T : float
            Isotherm to simulate

        P: float
           Isobar to simulate

        res : integer
              Specifies the number of data points to be simulated within the
              specified range.

        tol : scalar, optional
              Tolerance used in ``dual_equal``, if epsilon >= UBD - LBD that will
              terminate the routine.

        gtol : scalar, optional
              Minimum tolerance between hyperplane solution
              Note: The Dual solution is not perfect so a low tolerance is
              required, but a too low tol could potentially include points that do
              not truly lie on the equilibrium plane within the considered
              instability region.

        n_dual : scalar, optional
                Number of sampling points used in the tgo routine in solving LBD
                of the dual problem.
                Note: It is recommended to use at least ``100 + p.m['n'] * 100``

        phase_tol : scalar, optional
                    The minimum seperation between equilibrium planes required to
                    be considered a phase. Defaults to 0.001

        LLE_only : boolean, optional
                   If True then only phase seperation of same volume root
                   instability will be calculated.

        VLE_only : boolean, optional
                   If True then phase seperation of same volume root instability
                   will be ignored.

        Plot_Results : boolean, optional
                       If True the g_mix curve with tie lines will be plotted for
                       binary and ternary systems.

        Returns
        -------

        data_x_mph : dict containing vector containing of all the equilibrium
                     points in the isotherm/bar (VLE type only)

        data_x_ph : dict containing vector containing of all the equilibrium
                    points in the isotherm/bar  (VLE type and LLE type
                    (TODO separate))

        data_t : vector containing all the temperature points in the iso-
                     therm/bar

        data_p : vector containing all the pressure points in the isotherm/bar

        P_range : vector of size ``res``
                  Range of pressure points over the min/max of the isotherm/bar

        T_range : vector of size ``res``
                  Range of temperature points over the min/max of the
                   isotherm/bar

        r_ph_eq : list of size ``res`` containing ph_eq returns:
            ph_eq : dict containing keys for each phase in p.m['Valid phases'],
             ex:
                ph_eq[ph] : list containing composition vectors
                            Contains a list of equilibrium points of phase (ph)
                            seperations in the same volume root of the EOS
                            (ex. LLE type)

        r_mph_eq : list of size ``res`` containing mph_eq returns:
            mph_eq : list containing composition vectors
                     contains a list of equilibrium points of phase
                     seperations in different volume roots of the EOS (mph)
                     (ex. VLE type)

        r_mph_ph  : list of size ``res`` containing mph_ph returns:
            mph_ph : list containing strings
                     containts the phase string of the corresponding ``mph_eq``
                     equilibrium point
        """
        import numpy
        from ncomp import equilibrium_range as er

        if T is not None:
            iso_ind = numpy.where(numpy.array(p.m['T']) == T)
            data_p = numpy.array(p.m['P'])[iso_ind]
            PT_Range = [(min(data_p), max(data_p)),
                        (T, T)]
            data_t = None

        if P is not None:
            iso_ind = numpy.where(numpy.array(p.m['P']) == P)
            data_t = numpy.array(p.m['T'])[iso_ind]
            PT_Range = [(min(data_t), max(data_t)),
                        (P, P)]
            data_p = None

        # Organize data
        # VLE phases
        data_x_mph = {}
        for ph in p.m['Valid phases']:
            data_x_mph[ph] = []
            data_x_mph[ph].append([])  # Empty tuple for 0 index
            for comp_n in range(1, p.m['n']):
                data_x_mph[ph].append(numpy.array(p.m[ph][comp_n])[iso_ind])

        # LLE phases
        data_x_ph = {}
        for ph in p.m['Data phases']:
            data_x_ph[ph] = []
            data_x_ph[ph].append([])  # Empty tuple for 0 index
            for comp_n in range(1, p.m['n']):
                data_x_ph[ph].append(numpy.array(p.m[ph][comp_n])[iso_ind])

        # Find model outputs
        if not data_only:
            P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph = \
                er(g_x_func, s, p, PT_Range=PT_Range, n=n, res=res, tol=tol,
                   gtol=gtol, n_dual=n_dual, phase_tol=phase_tol,
                   LLE_only=LLE_only, VLE_only=VLE_only,
                   Plot_Results=Plot_Results)

        else:
            r_ph_eq, r_mph_eq, r_mph_ph = None, None, None

        return (P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph, data_x_mph,
                data_x_ph, data_t, data_p)


    def process_LLE_range(self, p, P_range, T_range, r_ph_eq):  # UNTESTED
        """
        Process the equilibrium points found in ncomp.equilibrium_range into
        plotable results for LLE type equilibrium points.

        Parameters
        ----------
        p : class
            Contains the dictionary describing the parameters.

        P_range: list
                 contains the pressure points at each model simulation

        T_range: list
                 contains the temperature points at each model simulation

        r_ph_eq : list containing ph_eq returns:
            ph_eq : dict containing keys for each phase in p.m['Valid phases'], ex:
                ph_eq[ph] : list containing composition vectors
                            Contains a list of equilibrium points of phase (ph)
                            seperations in the same volume root of the EOS
                            (ex. LLE type)

        Returns
        -------
        model_x_ph: dict containing equilibrium tie line vectors for each phase

        model_p_ph: dict containing pressure vectors at each tie line

        model_t_ph: dict containing temperature vectors at each tie line
        """

        model_x_ph = {}  # LLE type equilibrium tie lines
        model_p_ph = {}
        model_t_ph = {}
        for ph in p.m['Data phases']:
            model_x_ph[ph] = []
            model_p_ph[ph] = []
            model_t_ph[ph] = []

        for i in range(len(P_range)):
            for ph in p.m['Valid phases']:
                if len(r_ph_eq[i][ph]) > 0:  # Equilibrium point found
                    for j in range(len(r_ph_eq[i][ph])):
                        if len(r_ph_eq[i][ph][j]) > 1:  # discard single points
                            model_x_ph[ph].append(r_ph_eq[i][ph][j][0])
                            model_p_ph[ph].append(P_range[i])
                            model_t_ph[ph].append(T_range[i])
                            l = 0
                            for ph2 in p.m['Data phases']:
                                l += 1
                                if ph2 is not ph:
                                    try:
                                        model_x_ph[ph].append(
                                            r_ph_eq[i][ph][j][l])
                                        model_p_ph[ph].append(P_range[i])
                                        model_t_ph[ph].append(T_range[i])
                                    except(IndexError):
                                        model_x_ph[ph].append(None)
                                        model_p_ph[ph].append(None)
                                        model_t_ph[ph].append(None)
                            #for l in range(1, len(r_ph_eq[i][ph][j])):
                            #    for

                                # Attach a pressure and temperature
                                # point for each of these to keep dims

        return model_x_ph, model_p_ph, model_t_ph


    def process_VLE_range(self, p, P_range, T_range, r_mph_eq, r_mph_ph):
        """
        Process the equilibrium points found in ncomp.equilibrium_range into
        plotable results for VLE type equilibrium points.

        Parameters
        ----------

        p : class
            Contains the dictionary describing the parameters.

        P_range: list
                 contains the pressure points at each model simulation

        T_range: list
                 contains the temperature points at each model simulation

        r_mph_eq : list containing mph_eq returns:
            mph_eq : list containing composition vectors
                     contains a list of equilibrium points of phase
                     seperations in different volume roots of the EOS (mph)
                     (ex. VLE type)

        r_mph_ph  : list containing mph_ph returns:
            mph_ph : list containing strings
                     containts the phase string of the corresponding ``mph_eq``
                     equilibrium point

        Returns
        -------
        model_x_mph: dict containing equilibrium tie line vectors for each
                     phase

        model_p_mph: dict containing pressure vectors at each tie line

        model_t_mph: dict containing temperature vectors at each tie line
        """

        # Set empty containers for all equilibrium points
        model_x_mph = {}  # VLE type equilibrium container
        model_p_mph = {}
        model_t_mph = {}
        for ph in p.m['Valid phases']:
            model_x_mph[ph] = []
            model_p_mph[ph] = []
            model_t_mph[ph] = []

        for i in range(len(P_range)):
            print('r_mph_eq[i] = {}'.format(r_mph_eq[i]))
            print('r_mph_ph[i] = {}'.format(r_mph_ph[i]))
            print('model_p_mph = {}'.format(model_p_mph))
            if len(r_mph_eq[i]) > 0:  # Equilibrium point found
                for j in range(len(r_mph_eq[i])):
                    #if len(r_mph_eq[i][j]) > 1:  # discard single
                        print'r_mph_eq[i][j] ={}'.format(r_mph_eq[i][j])
                        print'r_mph_ph[i][j] ={}'.format(r_mph_ph[i][j])
                        # points
                        for l in range(len(r_mph_eq[i][j])):
                            print'r_mph_eq[i][j] ={}'.format(r_mph_eq[i][j][l])
                            print'r_mph_ph[i][j] ={}'.format(r_mph_ph[i][j][l])
                            if len(r_mph_eq[i][j][l]) > 1:  # discard single
                                for q in range(len(r_mph_eq[i][j][l])):
                                    model_x_mph[r_mph_ph[i][j][l][q]].append(
                                        r_mph_eq[i][j][l][q])
                                    model_p_mph[r_mph_ph[i][j][l][q]].append(P_range[i])
                                    model_t_mph[r_mph_ph[i][j][l][q]].append(T_range[i])
                            # Attach a pressure and temperature
                            # point for each of these to keep dims

        return model_x_mph, model_p_mph, model_t_mph


    def plot_iso_t_bin(self, T, p, data_p=None, data_x_mph=None,
                       data_x_ph=None, model_p_mph=None, model_x_mph=None,
                       model_p_ph=None, model_x_ph=None,
                       k=['All'], FigNo=None, plot_options=None,
                       plot_tie_lines=True, LLE_only=False, VLE_only=False):
        """
        Plot binary isotherms for the specified data and model ranges

        Parameters
        ----------

        T : float
            Temperature of isotherm to plot.

        p : class
            Contains the dictionary describing the parameters.

        data_p : vector
                 Pressure data values at each point

        data_x_mph : dict containing vector containing of all the equilibrium
                     points in the isotherm/bar (VLE type only)

                 Contains the composition data points at each data_p for
                 every valid phase

                 ex. data_x = {'x': [p.m['x'][1][25:36],  # x_1
                                     p.m['x'][2][25:36]], # x_2

                               'y': [p.m['y'][1][25:36],  # y_1
                                     p.m['y'][2][25:36]]  # y_2
                               }

        data_x_ph : dict containing vector containing of all the equilibrium
                    points in the isotherm/bar  (VLE type and LLE type
                    (TODO separate))

        model_x_mph: dict containing equilibrium tie line vectors for each
                     phase

        model_p_mph: dict containing pressure vectors at each tie line

        model_x_ph: dict containing equilibrium tie line vectors for each phase

        model_p_ph: dict containing pressure vectors at each tie line

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


        # VLE envelopes
        if not LLE_only:
            # Plot data points
            if (data_x_mph is not None) and (len(data_x_mph) > 0):
                for ph in k:
                    plot.plot(data_x_mph[ph][1], data_p, 'x',
                              label='{} data'.format(ph))
            # Plot model points
            if model_p_mph is not None and (len(model_p_mph) > 0):
                for ph in k:
                    print('model_x_mph[ph] = {}'.format(model_x_mph[ph]))
                    print('model_p_mph[ph] = {}'.format(model_p_mph[ph]))
                    # plot.plot(model_x[ph][1], model_p, '-',
                    #           label='{} model'.format(ph))
                    plot.plot(model_x_mph[ph], model_p_mph[ph], '-',
                              label='{} model VLE'.format(ph))

        # LLE envelopes
        if not VLE_only:
            # Plot data points
            if (data_x_ph is not None) and (len(data_x_ph) > 0):
                for ph in p.m['Data phases']:
                    plot.plot(data_x_ph[ph][1], data_p, 'x',
                              label='{} data'.format(ph))
            # Plot model points
            if model_p_ph is not None and (len(model_p_ph) > 0):
                for ph in k:
                    print('model_x_ph[ph] = {}'.format(model_x_ph[ph]))
                    print('model_p_ph[ph] = {}'.format(model_p_ph[ph]))
                    # plot.plot(model_x[ph][1], model_p, '-',
                    #           label='{} model'.format(ph))
                    plot.plot(model_x_ph[ph], model_p_ph[ph], '-',
                              label='{} model'.format(ph))

        plot.xlabel(r"$z_1$", fontsize=14)
        plot.ylabel(r"P (Pa)", fontsize=14)
        plot.title("{}-{} isotherm at T = {}".format(p.c[1]['name'][0],
                                                     p.c[2]['name'][0],
                                                     T))
        plot.legend()
        return

    def plot_iso_p_bin(self, P, p, data_t=None, data_x_mph=None,
                       data_x_ph=None, model_t_mph=None, model_x_mph=None,
                       model_t_ph=None, model_x_ph=None,
                       k=['All'], FigNo=None, plot_options=None,
                       plot_tie_lines=True, LLE_only=False, VLE_only=False):
        """
        Plot binary isobars for the specified data and model ranges

        Parameters
        ----------

        P : float
            Pressure of isobar to plot.

        p : class
            Contains the dictionary describing the parameters.

        data_t : vector
                 Temperature data values at each point

        data_x_mph : dict containing vector containing of all the equilibrium
                     points in the isotherm/bar (VLE type only)

                 Contains the composition data points at each data_p for
                 every valid phase

                 ex. data_x = {'x': [p.m['x'][1][25:36],  # x_1
                                     p.m['x'][2][25:36]], # x_2

                               'y': [p.m['y'][1][25:36],  # y_1
                                     p.m['y'][2][25:36]]  # y_2
                               }

        data_x_ph : dict containing vector containing of all the equilibrium
                    points in the isotherm/bar  (VLE type and LLE type
                    (TODO separate))rm/bar  (VLE type and LLE type
                    (TODO separate))

        model_x_mph: dict containing equilibrium tie line vectors for each
                     phase

        model_t_mph: dict containing temperature vectors at each tie line

        model_x_ph: dict containing equilibrium tie line vectors for each phase

        model_t_ph: dict containing temperature vectors at each tie line

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

        # VLE envelopes
        if not LLE_only:
            # Plot data points
            if (data_x_mph is not None) and (len(data_x_mph) > 0):
                for ph in k:
                    plot.plot(data_x_mph[ph][1], data_t, 'x',
                              label='{} data'.format(ph))
            # Plot model points
            if (model_t_mph is not None) and (len(model_x_mph) > 0):
                for ph in k:
                    # plot.plot(model_x[ph][1], model_p, '-',
                    #           label='{} model'.format(ph))
                    plot.plot(model_x_mph[ph], model_t_mph[ph], '-',
                              label='{} model'.format(ph))

        # LLE envelopes
        if not VLE_only:
            # Plot data points
            if (data_x_ph is not None) and (len(data_x_ph) > 0):
                for ph in p.m['Data phases']:
                    plot.plot(data_x_ph[ph][1], data_t, 'x',
                              label='{} data'.format(ph))
            # Plot model points
            if model_t_ph is not None and (len(model_t_ph) > 0):
                for ph in k:
                    # plot.plot(model_x[ph][1], model_p, '-',
                    #           label='{} model'.format(ph))
                    plot.plot(model_x_ph[ph], model_t_ph[ph], '-',
                              label='{} model'.format(ph))

        plot.xlabel(r"$z_1$", fontsize=14)
        plot.ylabel(r"P (Pa)", fontsize=14)
        plot.title("{}-{} isotherm at T = {}".format(p.c[1]['name'][0],
                                                     p.c[2]['name'][0],
                                                     P))
        plot.legend()
        return

    def plot_iso_tern(self, p, data_x_ph=None, data_x_mph=None,
                      model_x_ph=None, model_x_mph=None,
                      plot_phase_envelopes=False,
                      LLE_only=False, VLE_only=False):

        import ternary

        scale = 1.0
        figure, tax = ternary.figure(scale=scale)

        # Draw Boundary and Gridlines
        tax.boundary(linewidth=1.5)
        tax.gridlines(color="black", multiple=0.1)

        # Set Axis labels and Title
        fontsize = 16
        tax.set_title("Various Lines", fontsize=fontsize)
        tax.left_axis_label("$x_3$", fontsize=fontsize)
        tax.right_axis_label("$x_2$", fontsize=fontsize)
        tax.bottom_axis_label("$x_1$", fontsize=fontsize)

        # Plot data

        # VLE
        if not LLE_only:
            for ph in p.m['Valid phases']:
                for i in range(len(data_x_mph)):
                    if i < (len(data_x_mph)-1):
                        tax.line(data_x_mph[ph][i],  # ex. x
                                 data_x_mph[ph][i + 1],  # ex. xII
                                 linewidth=1.0, marker='x',
                                 # color='green',
                                 linestyle="--", label='Data')

        # LLE
        if not VLE_only:
            for ph in p.m['Valid phases']:
                for lph in p.m['Data phases']:
                    if lph not in p.m['Valid phases']:
                        for i in range(len(data_x_ph)):
                            tax.line(data_x_ph[ph][i],  # ex. x
                                     p.m[lph][i],  # ex. xII
                                     linewidth=1.0, marker='x',
                                     # color='green',
                                     linestyle="--", label='Data')



        # Plot model results

        # LLE
        if not VLE_only:
            for i in range(len(p.m['Valid phases'])):
                for j in range(i + 1, len(p.m['Valid phases'])):
                    ph1 = p.m['Valid phases'][i]
                    ph2 = p.m['Valid phases'][j]
                    tax.line(data_x_ph[ph1][i],  # ex. x
                             data_x_ph[ph2][i],  # ex. xII
                             linewidth=1.0, marker='.',
                             linestyle="-", label='Model {}-{} LLE'.format(ph1,
                                                                         ph2))
        # VLE
        if not LLE_only:
            for i in range(len(data_x_mph)):
                if i < len(data_x_mph):
                    for ph in p.m['Valid phases']:
                        tax.line(data_x_ph[ph][i],  # ex. x
                                 data_x_ph[ph][i + 1],  # ex. xII
                                 linewidth=1.0, marker='.',
                                 linestyle="-", label='Model {} VLE'.format(ph))

        tax.ticks(axis='lbr', multiple=0.1, linewidth=1)

        tax.legend()
        tax.show()
        return


def plot_ep(func, x_r, s, p, args=()):
    """
    Plot the speficied single var input error function over a range size x_r
    """
    from matplotlib import pyplot as plot
    #% Binary
    if p.m['n'] == 2:
        from scipy import linspace
        X_r = linspace(1e-5, 0.9999, x_r)
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
        #plot.show()

    if p.m['n'] == 3:
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plot
        from matplotlib import cm
        import numpy
        x_range = numpy.linspace(1e-15, 1.0, x_r)
        y_range = numpy.linspace(1e-15, 1.0, x_r)
        xg, yg = numpy.meshgrid(x_range, y_range)
        func_r = numpy.zeros((x_r, x_r))
        for i in range(xg.shape[0]):
            for j in range(yg.shape[0]):
                X = [xg[i, j], yg[i, j]]  # [x_1, x_2]
                if sum(X) > 1.0:  # 1.0:
                    func_r[i, j]  = None
                else:
                    f_out = func(X, *args)  # Scalar outputs
                    func_r[i, j] = numpy.float64(f_out)


        # Plots
        fig = plot.figure()
        ax = fig.gca(projection='3d')
        X, Y = xg, yg

        # Gibbs phase surfaces
        Z = func_r

        if True:
            print 'numpy.min(Z) = {}'.format(numpy.nanmin(Z))

            cset = ax.contourf(X, Y, Z, zdir='z', offset=numpy.nanmin(Z)-0.05,
                               cmap=cm.coolwarm)
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.3,
                            cmap=cm.coolwarm)
        if False:
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                   cmap=cm.coolwarm, linewidth=0,
                                   antialiased=True, alpha=0.5, shade = True)

            fig.colorbar(surf, shrink=0.5, aspect=5)

        ax.set_xlabel('$x_1$')
        ax.set_xlim(0, 1)
        ax.set_ylabel('$x_2$')
        ax.set_ylim(0, 1)
        ax.set_zlabel('$\Delta g$', rotation=90)
        plot.show()

    return


if __name__ == '__main__':
     plot_options = {'text.usetex' : True, # Options for all plots
                     'font.size' : 12,
                     'font.family' : 'lmodern',
                     'text.latex.unicode': True
                     }

     import ternary

     scale = 1.0
     figure, tax = ternary.figure(scale=scale)

     # Draw Boundary and Gridlines
     tax.boundary(linewidth=1.5)
     tax.gridlines(color="black", multiple=0.1)

     # Set Axis labels and Title
     fontsize = 16
     tax.set_title("Various Lines", fontsize=fontsize)
     tax.left_axis_label("$x_3$", fontsize=fontsize)
     tax.right_axis_label("$x_2$", fontsize=fontsize)
     tax.bottom_axis_label("$x_1$", fontsize=fontsize)

     # Draw an arbitrary line, ternary will project the points for you
     p1 = (0.5, 0.5, 0.1)
     #    (x_1, x_2, x_3)
     p1 = (0.3, 0.5, 0.2)
     p2 = (0.1, 0.1, 0.8)

     tax.line(p1, p2, linewidth=1.0, marker='x',  # color='green',
              linestyle="--", label='Data')

     p1 = (0.4, 0.4, 0.2)
     p2 = (0.15, 0.15, 0.7)
     tax.line(p1, p2, linewidth=1.0, marker='.',  # color='green',
              linestyle="-", label='Model')

     tax.ticks(axis='lbr', multiple=0.1, linewidth=1)

     tax.legend()
     tax.show()