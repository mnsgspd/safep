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
        """Assign a material property value of 1 if no material property value is set."""
        self.EA=EA
        self.EI=EI

    def add_node(self,x,y):
        print(x,y)
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            print(f"nodes[{len(self.nodes)}] = [{x},{y}]")
            node_num=len(self.nodes)
            self.nodes.append(Node(node_num,x,y));
        elif isinstance(x, list) and isinstance(y, list) and len(x) == len(y):
            for i in range(len(x)):
                if isinstance(x[i],(int,float)) and isinstance(y[i],(int,float)):
                    print(f"nodes[{len(self.nodes)}] = [{x[i]},{y[i]}]")
                    node_num=len(self.nodes)
                    self.nodes.append(Node(node_num,x[i],y[i]));
                else:
                    raise ValueError(f"Invalid node coordinates: {x[i]}, {y[i]}")
        else:
            raise ValueError("Invalid arguments(add_node).\
            Expected two numbers or two lists of same length.")
        
    def add_beam(self,n0,n1,EA=None, EI=None, hinge0=False,hinge1=False,\
                 element_type="beam"):
        if EA is None:
            EA = self.EA
        if EI is None:
            EI = self.EI
        if isinstance(n0,int) and isinstance(n1,int):
            print(f"elements[{len(self.elements)}] =[{n0},{n1}] : {element_type}")
            self.elements.append(Element(self.nodes[n0],self.nodes[n1],EA,EI,hinge0,hinge1))
            self.elements[-1].calc_geometry()
        elif isinstance(n0, list) and isinstance(n1, list) and len(n0) == len(n1):
            for i in range(len(n0)):
                if(n0[i]<0 or len(self.nodes) <= n0[i] or \
                   n1[i]<0 or len(self.nodes) <= n1[i]):
                    raise ValueError(f"Invalid arguments(add_{element_type}).\
                    The node number must be less than {len(self.nodes)}.")
                if isinstance(n0[i],int) and isinstance(n1[i],int):
                    print(f"node[{len(self.nodes)}] = [{x[i]},{y[i]}]")
                    self.elements.append(Element(n0[i],n1[i],EA,EI,hinge0,hinge1))
                    self.elements[-1].calc_geometry()
                else:
                    raise ValueError(f"Invalid element number: {x[i]}, {y[i]}")
        else:
            raise ValueError("Invalid arguments(add_{element_type}).\
            Expected two numbers or two lists of same length.")

    def add_truss_member(self,n0,n1,EA=None):
        self.add_beam(n0,n1,EA,hinge0=True,hinge1=True,element_type="truss_member")

    def set_node_boundary_condition(self,node_num,x_specified=None,y_specified=None,
                                    theta_specified=None,u=None,v=None,theta=None):
        if(len(self.nodes)<=node_num):
            print("set_node_boundary_condtion error")
            print("Check node number.")
            sys.exit()
        self.nodes[node_num].boundary_condition(x_specified,y_specified,theta_specified,\
                                                u,v,theta)        

    def plot_structure(self,show_node_number=True,show_element_number=True):
        self.plotter.plot_structure(show_node_number=True,show_element_number=True)

    def plot_deformation(self,deformation_scale=None):
        self.plotter.plot_deformation(deformation_scale)

    def solve_deformation(self):
        self.update_hinge()
        self.solver.solve()

    def apply_load(self,node_num,x_force=None,y_force=None,moment_force=None):
        self.nodes[node_num].set_load(x_force,y_force,moment_force)

    def convert_hinge(self,node_num):
        self.nodes[node_num].is_all_hinge=True
        
    def update_hinge(self):
        for i,e in enumerate(self.elements):
            n0=e.node[0]
            n1=e.node[1]
            n0.elements.append(i)
            n1.elements.append(i)
            if(e.hinge[0]==True):
                n0.num_hinge+=1
            if(e.hinge[1]==True):
                n1.num_hinge+=1
        for n in self.nodes:
            if len(n.elements)-n.num_hinge==1:
                n.is_all_hinge=True
            elif len(n.elements)-n.num_hinge==0:
                e=self.elements[n.elements[0]]
                if e.node[0]==n:
                    e.hinge[0]=False
                else:
                    e.hinge[1]=False
            elif len(n.elements)-n.num_hinge>1 and n.is_all_hinge==True:
                for e in n.elements:
                    for j in range(2):
                        if e.node[j]==n and e.hinge[j]==False:
                            e.hinge[j]==True
                            num_hinge+=1
                            break
                    if(len(n.e.ements)-n.num_hinge<=1):
                        break
        for e in self.elements:
            e.update_hinge_type()
            
    def apply_element_force(self,element_num,force_type,force):
        if(force_type=="concentrated"):
            self.elements[element_num].calc_concentrated_force(force)
        elif(force_type=="distributed"):
            self.elements[element_num].calc_distributed_force(force)
        else:
            print("calc_element_force error")
            print(force_type,"is unknown force type.")
            sys.exit()
