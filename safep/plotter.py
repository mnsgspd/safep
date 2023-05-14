import matplotlib.pyplot as plt

class Plotter:
    def __init__(self,nodes,elements):
        self.nodes=nodes
        self.elements=elements

        self.support_size=None
        self.hinge_size=None
        self.hinge_shift=False
        self.hinge_offset=False
        self.font_size=20.0
        self.member_width=3.0
        self.int_format="{:.8g}"
        self.float_format="{:d}"
    
