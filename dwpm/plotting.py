import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from mpl_toolkits.mplot3d import Axes3D

#%% Binary plots
def plot_g_binary(
    g_func,
    g_func_args=None,
    xlabel='$x_1$',
    ylabel='Gibbs Free Energy, $g(x)$',
    sol_label='Solution tangent plane',
    x_min=1e-4,
    x_max=1 - 1e-4,
    num_points=200,
    eq_points=None,
    # Default colors:
    color_surface=(1.0, 0.498, 0.0),  # Approx "tab:orange"
    color_planes=(0.1216, 0.4667, 0.7059),  # Approx "tab:blue"
    color_vertical=(0.5, 0.5, 0.5)         # Gray for vertical lines
):
    """
    Plots a one-dimensional Gibbs free energy function G(x) for a binary mixture,
    with optional equilibrium highlights.

    Parameters
    ----------
    g_func : callable
        A Gibbs free energy function of one variable, G(x),
        with domain x in (0,1).
    x_min : float, optional
        Minimum x-value for plotting (to avoid singularities at x=0).
    x_max : float, optional
        Maximum x-value for plotting (to avoid singularities at x=1).
    num_points : int, optional
        Number of points in the plotting grid from x_min to x_max.
    eq_points : list or tuple of floats, optional
        Typically two solution points [x_alpha, x_beta] for an equilibrium split.
        If provided:
          - Draws gray dashed vertical lines at each solution
          - Attempts to label them below the x-axis as x^α, x^β
          - Draws the chord (straight line) between (x_alpha, G(x_alpha)) and
            (x_beta, G(x_beta)) both across the entire domain and restricted to
            the interval [x_alpha, x_beta].
    color_surface : tuple of 3 floats
        RGB color for the G(x) surface curve. Default is approximately tab:orange.
    color_planes : tuple of 3 floats
        RGB color for the chords. Default is approximately tab:blue.
    color_vertical : tuple of 3 floats
        RGB color for the vertical lines. Default is gray.

    Returns
    -------
    None
        Displays a matplotlib figure.

    Example
    -------
    >>> def example_g(x):
    ...     # simple polynomial curve, x is as vector of size 2
    ...     return 0.01*x[0] + (x[0]-0.5)**4 - 0.2*(x[0]-0.5)**2
    >>> plot_g_binary(example_g, eq_points=[0.2, 0.8])
    """
    x_vals = np.linspace(x_min, x_max, num_points)
    g_vals = [g_func([x, 1-x], *g_func_args) for x in x_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(x_vals, g_vals, color=color_surface, label='G(x)')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.title('Change in reduced Gibbs Free Energy')
    plt.grid(True)

    # If eq_points provided (often two), mark them
    if eq_points is not None and len(eq_points) > 0:
        # Sort for consistent left->right
        eq_points_sorted = sorted(eq_points)

        # Vertical lines & chord(s)
        for i, x_star in enumerate(eq_points_sorted):
            # Draw a dashed vertical line
            plt.axvline(x_star, color=color_vertical, linestyle='--', alpha=0.8)

            # Attempt to label below the x-axis
            # If too clipped in your environment, remove or adjust the offset
            label_text = r"$x^{\alpha}$" if i == 0 else r"$x^{\beta}$"
            # We'll place the annotation at (x_star, 'bottom of axes'):
            ax = plt.gca()
            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
            plt.annotate(
                label_text,
                xy=(x_star, 0),  # x in data coords, y=0 in axes coords
                xycoords=trans,
                xytext=(0, -25),  # shift downward by 25 points
                textcoords='offset points',
                ha='center',
                va='top',
                color=color_vertical,
                annotation_clip=False  # allow drawing outside axes
            )

        # If we have at least two points, draw the chord
        if len(eq_points_sorted) >= 2:
            x_alpha = eq_points_sorted[0]
            x_beta = eq_points_sorted[1]
            G_alpha = g_func([x_alpha, 1 - x_alpha], *g_func_args)
            G_beta = g_func([x_beta, 1 - x_beta], *g_func_args)

            if abs(x_beta - x_alpha) > 1e-14:
                # Slope of the chord
                m = (G_beta - G_alpha) / (x_beta - x_alpha)
                print(f'm = {m}')
                # Function for the chord across domain
                def chord(x):
                    return G_alpha + m*(x - x_alpha)

                # (a) Full-domain chord
                chord_full = [chord(x) for x in x_vals]
                plt.plot(x_vals, chord_full,
                         color=color_planes, linestyle='--', alpha=0.6,
                         #label='Chord (full domain)'
                         )

                # (b) Restricted chord from x_alpha to x_beta
                x_sub = np.linspace(x_alpha, x_beta, 50)
                chord_sub = [chord(x) for x in x_sub]
                plt.plot(x_sub, chord_sub,
                         color=color_planes, linestyle='-', linewidth=2.0,
                         label=sol_label)

    plt.legend()
    plt.show()


#%% Ternary plots

def plot_g_ternary(
    g_func,
    eq_planes=None,
    # Domain sampling
    num_points=100,
    # Colors & style
    color_surface=(1.0, 0.498, 0.0),
    alpha_surface=0.6,
    color_planes=(0.1216, 0.4667, 0.7059),
    alpha_planes=0.4,
    # Labels
    x_label=r"$x_1$",
    y_label=r"$x_2$",
    z_label=r"$G(x_1,x_2)$",
    title=None
):
    """
    Plots a ternary Gibbs free energy surface G(x1, x2, 1 - x1 - x2) over the
    simplex region { 0 <= x1, x2, x1+x2 <= 1 }, plus optional linear planes
    (tie-lines) given by eq_planes.

    Parameters
    ----------
    g_func : callable
        A function g_func(x) -> float, where x is a 3-component composition
        [x1, x2, x3], with x3 = 1 - x1 - x2. Must be valid for 0 <= x1, x2,
        x1+x2 <= 1.
    eq_planes : list of lists, optional
        Each element is a list [G_p, x1, lambda1, x2, lambda2], representing
        a plane T(x,y) = G_p + lambda1*(x - x1) + lambda2*(y - x2).
        The plane is plotted only in the valid ternary domain (x + y <= 1).
    num_points : int, optional
        The number of points for each dimension in the mesh (x1, x2).
        The total grid size is num_points^2.
    color_surface : str or tuple
        The color for the main Gibbs surface. Default is 'red'.
    alpha_surface : float
        Transparency alpha for the main surface. Default is 0.3.
    color_planes : str or tuple
        The color for the tie-plane surfaces. Default is 'blue'.
    alpha_planes : float
        Transparency alpha for the planes. Default is 0.2.
    x_label : str, optional
        Label for the x-axis (x1).
    y_label : str, optional
        Label for the y-axis (x2).
    z_label : str, optional
        Label for the z-axis (G).
    title : str, optional
        Title for the 3D plot.

    Returns
    -------
    None
        Displays a matplotlib 3D figure with the main G surface and any planes.

    Example
    -------
    >>> def example_g(x):
    ...     x1, x2, x3 = x
    ...     # A toy function, say ideal mixture + random
    ...     import math
    ...     return x1*math.log(x1 + 1e-12) + x2*math.log(x2 + 1e-12) \
    ...            + x3*math.log(x3 + 1e-12)
    ...
    >>> planes = [
    ...     [0.5, 0.3, 1.2, 0.4, -0.9]  # G_p, x1, lam1, x2, lam2
    ... ]
    >>> plot_g_ternary(example_g, eq_planes=planes)
    """

    # 1) Build mesh of x1, x2 from 0..1, skipping values where x1+x2>1
    x_vals = np.linspace(1e-12, 1 - 1e-12, num_points)
    y_vals = np.linspace(1e-12, 1 - 1e-12, num_points)
    X, Y = np.meshgrid(x_vals, y_vals)

    # Prepare Z array for the main Gibbs surface
    Z = np.zeros_like(X, dtype=float)
    Z[:] = np.nan  # default to NaN outside domain

    # For each grid point, if x1+x2 <= 1 => x3=1-(x1+x2), evaluate G
    for i in range(num_points):
        for j in range(num_points):
            x1 = X[i,j]
            x2 = Y[i,j]
            if x1 + x2 <= 1.0 - 1e-12:
                x3 = 1.0 - x1 - x2
                Z[i,j] = g_func([x1, x2, x3])
            else:
                Z[i,j] = np.nan  # outside simplex

    # 2) Create a figure + 3D axis
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # 3) Plot the main surface
    ax.plot_surface(
        X, Y, Z,
        color=color_surface,
        alpha=alpha_surface,
        linewidth=0,
        rcount=100, ccount=100
    )

    # 4) If eq_planes is given, plot each plane
    if eq_planes is not None:
        for plane in eq_planes:
            # plane = [G_p, x1, lambda1, x2, lambda2]
            G_p, refx1, lam1, refx2, lam2 = plane

            # We'll build a Z_plane array
            Zp = np.zeros_like(X, dtype=float)
            Zp[:] = np.nan
            for i in range(num_points):
                for j in range(num_points):
                    x1 = X[i,j]
                    x2 = Y[i,j]
                    if x1 + x2 <= 1.0:
                        # T(x,y) = G_p + lam1*(x - refx1) + lam2*(y - refx2)
                        Zp[i,j] = G_p + lam1*(x1 - refx1) + lam2*(x2 - refx2)
                    else:
                        Zp[i,j] = np.nan

            # Plot the plane surface
            ax.plot_surface(
                X, Y, Zp,
                color=color_planes,
                alpha=alpha_planes,
                linewidth=0
            )

    # 5) Add axis labels & title
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(np.nanmin(Z), np.nanmax(Z))  # or you can let MPL autoscale
    ax.set_title(title)

    plt.show()