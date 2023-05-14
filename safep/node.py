class Node:
    def __init__(self,x,y,
                 x_force=0.0,y_force=0.0,moment_force=0.0,\
                 u_given=False,v_given=False,theta_given=False,\
                 u=0.0,v=0.0,theta=0.0,is_all_hinge=False):
        self.x = x
        self.y = y
        self.x_force=x_force
        self.y_force=y_force
        self.moment_force=moment_force
        self.u_given=u_given
        self.v_given=v_given
        self.theta_given=theta_given
        self.u=u
        self.v=v
        self.theta=theta
        
        """Offset between node position and text input position        
        """
        self.offset_x=None
        self.offset_y=None

        """Indicates whether node connection is all hinged connections or not."
        """
        self.is_all_hinge=is_all_hinge
        self.element=[]
        self.num_hinge=0
