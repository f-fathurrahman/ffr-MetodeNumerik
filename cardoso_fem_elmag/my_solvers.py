import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:
    def __init__(self, nodes, triElements, elements, boundaryConditions):
        self.nodes = nodes
        self.triElements = triElements
        self.elements = elements
        self.boundaryConditions = boundaryConditions
        self.numNodes = len(nodes)
        self.globalMatrix = lil_matrix((self.numNodes, self.numNodes))
        self.globalVector = np.zeros(self.numNodes)
        self.potential = np.zeros(self.numNodes)

    def apply_boundary_conditions(self):
        self.fixedNodes = []
        self.fixedValues = np.full(self.numNodes, np.nan)
        for group, props in self.boundaryConditions.items():
            self.fixedNodes.extend(props["nodes"])
            self.fixedValues[props["nodes"]] = props["potential"]
        #
        self.fixedNodes = np.array(self.fixedNodes)
        self.freeNodes = np.where(np.isnan(self.fixedValues))[0]
        self.fixedValues[np.isnan(self.fixedValues)] = 0
    
    def solve(self):
        # Eliminate fixed nodes from the system
        matrix_reduced = self.globalMatrix[self.freeNodes][:,self.freeNodes]
        vector_reduced = self.globalVector[self.freeNodes] - \
                         self.globalMatrix[self.freeNodes][:, self.fixedNodes] @ self.fixedValues[self.fixedNodes]
        self.potential[self.freeNodes] = spsolve(matrix_reduced,vector_reduced)
        self.potential[self.fixedNodes] = self.fixedValues[self.fixedNodes]


class ElectrostaticSolver(Solver):
    def assemble_global_matrix_and_vector(self):
        for e, indices in enumerate(self.triElements):
            for i, iGlobal in enumerate(indices):
                self.globalVector[iGlobal] += self.elements[e].Q[i]
                for j, jGlobal in enumerate(indices):
                    self.globalMatrix[iGlobal, jGlobal] += self.elements[e].C[i, j]
    
    def get_potential(self):
        return self.potential
    
    def get_electric_field(self):
        Ex, Ey, modE, xc, yc = [], [], [], [], []
        for element, indices in zip(self.elements, self.triElements):
            element.setNodePotentials(self.potential[indices])
            Ex.append(element.Ex)
            Ey.append(element.Ey)
            modE.append(element.modE)
            xc.append(element.xc)
            yc.append(element.yc)
        return np.array(Ex), np.array(Ey), np.array(modE), np.array(xc), np.array(yc)


