import numpy as np
from optim_golden_ratio import *

def ladder_length(x):
    h = 4.0
    d = 4.0
    y = h*x/(x - d)
    return np.sqrt(x**2 + y**2)

# x should be larger than d
import matplotlib.pyplot as plt
xgrid = np.linspace(4.1, 20.0, 1000)
ygrid = ladder_length(xgrid)
plt.clf()
plt.plot(xgrid, ygrid)
plt.grid(True)
plt.tight_layout()
plt.savefig("IMG_chapra_exercise_13_22.pdf")

xopt, _ = optim_golden_ratio(ladder_length, 4.1, 20.0, sign=-1)
print("Minimum ladder length: ", ladder_length(xopt))
# we don't use minus here