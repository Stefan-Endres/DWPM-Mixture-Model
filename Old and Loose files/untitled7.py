#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""


import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sixhump(x):
    return (4 - 2.1*x[0]**2 + x[0]**4 / 3.) * x[0]**2 + x[0] * x[1] + (-4 + \
        4*x[1]**2) * x[1] **2

x = np.linspace(-2, 2)
y = np.linspace(-1, 1)
xg, yg = np.meshgrid(x, y)

zpoints = sixhump([xg, yg])

#plt.figure()  # simple visualization for use in tutorial
#plt.imshow(sixhump([xg, yg]))
#plt.colorbar()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(xg, yg, zpoints, rstride=1, cstride=1,
                       cmap=plt.cm.jet, linewidth=0, antialiased=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)')
ax.set_title('Six-hump Camelback function')

#%%

import numpy
from scipy import optimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = numpy.linspace(-5, 5)
y = numpy.linspace(-5, 5)
xg, yg = numpy.meshgrid(x, y)
z = numpy.zeros((50,50))

from random import randint as rnd
for i in range(xg.shape[0]):
    for j in range(xg.shape[0]):
        # Do stuff to find "z"
        z[i,j] = rnd(1,5)


#zpoints = sixhump([xg, yg])

#plt.figure()  # simple visualization for use in tutorial
#plt.imshow(sixhump([xg, yg]))
#plt.colorbar()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(xg, yg, z, rstride=1, cstride=1,
                       cmap=plt.cm.jet, linewidth=0, antialiased=False)

ax.set_xlabel('r')
ax.set_ylabel('s')
ax.set_zlabel('$\epsilon$',rotation='vertical')
#ax.set_title('Six-hump Camelback function')
