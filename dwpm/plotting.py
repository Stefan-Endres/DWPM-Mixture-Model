import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms


#%% Binary plots
def plot_g_binary(
    g_func,
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
    ...     # simple polynomial curve
    ...     return 0.01*x + (x-0.5)**4 - 0.2*(x-0.5)**2
    >>> plot_g_binary(example_g, eq_points=[0.2, 0.8])
    """
    x_vals = np.linspace(x_min, x_max, num_points)
    g_vals = [g_func(x) for x in x_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(x_vals, g_vals, color=color_surface, label='G(x)')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Change in reduced Gibbs Free Energy')
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
            G_alpha = g_func(x_alpha)
            G_beta = g_func(x_beta)

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