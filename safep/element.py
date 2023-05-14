class Element:
    def __init__(self,node0,node1,EA,EI,hinge0,hinge1):
        self.node = [node0,node1]
        self.hinge = [hinge0,hinge1]
        self.EA = EA
        self.EI = EI
        self.alpha=0.0
        self.T0=0.0
        self.T1=0.0
        self.fixedend_force = np.zeros(6)
        self.temp_force = np.zeros(6)

        self.offset_x=None
        self.offset_y=None
        
        self.hinge_type=0
        if(hinge0==True):
            self.hinge_type+=1
        if(hinge1==True):
            self.hinge_type+=2
        # 0 : member with no hinges
        # 1 : member with a hinge at its begining
        # 2 : member with a hinge at its end
        # 3 : member with hinges at both ends
    
