from safep.support import Support

class Node:
    def __init__(self,node_num,x,y,
                 x_force=0.0,y_force=0.0,moment_force=0.0,\
                 u_specified=False,v_specified=False,theta_specified=False,\
                 u=0.0,v=0.0,theta=0.0,is_all_hinge=False):
        self.node_num=node_num
        self.pos = [x,y]
        self.force=[x_force,y_force,moment_force]
        self.disp_specified=[u_specified,v_specified,theta_specified] #True or False
        self.disp=[u,v,theta]
        
        """Offset between node position and text input position        
        """
        self.offset=[None,None]
        """Indicates whether node connection is all hinged connections or not."
        """
        self.is_all_hinge=is_all_hinge
        self.elements=[]
        self.num_hinge=0

        self.support=Support(self.disp_specified)

    def print_info(self):
        print("--------")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")
        print("--------")

    def boundary_condition(self,x_specified,y_specified,theta_specified,\
                           u,v,theta):
        if(x_specified is not None):
            self.disp_specified[0] = x_specified
        if(y_specified is not None):
            self.disp_specified[1] = y_specified
        if(theta_specified is not None):
            self.disp_specified[2] = theta_specified
        if(u is not None):
            if(x_specified==True):
                self.disp[0]=u
            else:
                print(f"nodes[{self.node_num}]: x_specified!=True, but u is prescribed.")
                sys.exit()
        if(v is not None):
            if(y_specified==True):
                self.disp[1]=v
            else:
                print(f"nodes[{self.node_num}]: y_specified!=True, but v is prescribed.")
                sys.exit()
        if(theta is not None):
            if(theta_specified==True):
                self.disp[2]=theta
            else:
                print(f"nodes[{self.node_num}]: theta_specified!=True, but theta is prescribed.")
                sys.exit()                
        self.support.update(self.disp_specified)
        
    def set_load(self,x_force,y_force,moment_force):
        if x_force is not None:
            self.force[0] = x_force
        if y_force is not None:
            self.force[1] = y_force
        if moment_force is not None:
            self.force[2] = moment_force

