#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
#%%
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D

fig = pl.figure()
ax = Axes3D(fig)
X = np.arange(-4, 4, 0.25)
Y = np.arange(-4, 4, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X ** 2 + Y ** 2)
Z = np.sin(R)

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=pl.cm.hot)
ax.contourf(X, Y, Z, zdir='z', offset=-2, cmap=pl.cm.hot)
ax.set_zlim(-2, 2)

pl.show()
#%%

import numpy as np
r, theta = np.mgrid[0:10, -np.pi:np.pi:10j]
x = r * np.cos(theta)
y = r * np.sin(theta)
z = np.sin(r)/r

from mayavi import mlab
mlab.mesh(x, y, z, colormap='gist_earth', extent=[0, 1, 0, 1, 0, 1])

mlab.mesh(x, y, z, extent=[0, 1, 0, 1, 0, 1], representation='wireframe', 
          line_width=1, color=(0.5, 0.5, 0.5))
          
#%%
import numpy as np
import matplotlib.pyplot as pl

def f(x,y):
    return (1 - x / 2 + x**5 + y**3) * np.exp(-x**2 -y**2)

n = 256
x = np.linspace(-3, 3, n)
y = np.linspace(-3, 3, n)
X,Y = np.meshgrid(x, y)

pl.axes([0.025, 0.025, 0.95, 0.95])

pl.contourf(X, Y, f(X, Y), 8, alpha=.75, cmap=pl.cm.hot)
C = pl.contour(X, Y, f(X, Y), 8, colors='black', linewidth=.5)
pl.clabel(C, inline=1, fontsize=10)

pl.xticks(())
pl.yticks(())
pl.show()