# Execute some setup in another file
exec(open("setup_problem_01.py", encoding="utf-8").read())
filesave = "IMG_basis_lin_01.pdf"

#exec(open("setup_problem_02.py", encoding="utf-8").read())
#filesave = "IMG_basis_lin_02.pdf"

Lnum = 2.0

import numpy as np
import matplotlib.pyplot as plt
NptsPlot = 200
xplt = np.linspace(0.0, Lnum, 200)
yplt = np.zeros(NptsPlot)
plt.clf()
for ibasis in range(Nnodes):
    for i in range(NptsPlot):
        yplt[i] = Nfuncs[ibasis].subs({x: xplt[i], L: Lnum})
    plt.plot(xplt, yplt, label="Basis-" + str(ibasis+1))
plt.legend()
plt.grid(True)
plt.savefig(filesave, dpi=150)
