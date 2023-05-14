from safep.solve_deformation import SolveDeformation
from safep.plotter import Plotter
from safep.element import Element
from safep.node import Node

class StructuralAnalysis:
    def __init__(self, EA=1.0, EI=1.0):
        self.nodes=[]
        self.elements=[]
        self.solver = SolveDeformation(self.nodes,self.elements)
        self.plotter = Plotter(self.nodes,self.elements)
