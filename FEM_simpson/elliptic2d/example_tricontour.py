# https://scicomp.stackexchange.com/questions/31463/efficiently-plot-a-finite-element-mesh-solution-with-matplotlib

import numpy as np
import matplotlib.pyplot as plt

nodes= np.array([
        [0.000, 0.000],
        [1.000, 0.000],
        [2.000, 0.500],
        [0.000, 1.000],
        [1.000, 1.000],
        [1.750, 1.300],
        [1.000, 1.700]])
eles = np.array([
        [1, 2, 5],
        [5, 4, 1],
        [2, 3, 6],
        [6, 5, 2],
        [4, 5, 7],
        [5, 6, 7]])
node_vals = [1, 2, 1, 2, 7, 4, 5]

x, y = nodes.T
plt.tricontourf(x, y, eles - 1, node_vals, 12)
plt.colorbar()
plt.show()