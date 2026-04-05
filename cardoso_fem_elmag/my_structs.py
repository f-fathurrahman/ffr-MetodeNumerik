import numpy as np

class Element:
    def __init__(self):
        self.x = self.y = None
        self.a = self.b = self.c = np.zeros(3)
        self.Delta = 0

    def setNodes(self, x, y):
        self.x, self.y = x, y
        self.xc = (x[0] + x[1] + x[2])/3
        self.yc = (y[0] + y[1] + y[2])/3
        self.a = [ x[1]*y[2] - x[2]*y[1],
                   x[2]*y[0] - x[0]*y[2],
                   x[0]*y[1] - x[1]*y[0] ]
        self.b = [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
        self.c = [x[2] - x[1], x[0] - x[2], x[1] - x[0]]
        self.Delta = abs(self.b[1]*self.c[2] - self.b[2]*self.c[1]) / 2


class ElectrostaticElement(Element):
    def __init__(self, eps=8.854187817e-12, rho=0):
        super().__init__()
        self.jpg = eps
        self.rho = rho
        self.C = np.zeros((3, 3))
        self.Q = np.zeros(3)
    
    def setNodes(self, x, y):
        super().setNodes(x, y)
        fac = self.jpg / (4 * self.Delta)
        self.C = fac * (np.outer(self.b, self.b) + np.outer(self.c, self.c))
        self.Q = (self.rho * self.Delta / 3) * np.ones(3)
    
    def setNodePotentials(self, V):
        self.Ex = -np.dot(self.b, V) / (2 * self.Delta)
        self.Ey = -np.dot(self.c, V) / (2 * self.Delta)
        self.modE = np.hypot(self.Ex, self.Ey)

    def setProperties(self, eps, rho):
        self.jpg = eps # why jpg ???
        self.rho = rho

