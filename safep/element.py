import numpy as np
class Element:
    def __init__(self,node0,node1,EA,EI,hinge0,hinge1):
        self.node = [node0,node1]
        self.hinge = [hinge0,hinge1]
        self.theta = [0.0,0.0]
        self.EA = EA
        self.EI = EI
        self.alpha=0.0
        self.T0=0.0
        self.T1=0.0
        self.fixedend_force = np.zeros(6)
        self.temp_force = np.zeros(6)

        self.offset=[None,None]
        
        self.hinge_type=0
        if(self.hinge[0]==True):
            self.hinge_type+=1
        if(self.hinge[1]==True):
            self.hinge_type+=2
        # 0 : member with no hinges
        # 1 : member with a hinge at its begining
        # 2 : member with a hinge at its end
        # 3 : member with hinges at both ends

        #parameter for plot
        #self.member_width=None
        #self.member_color=None
        self.data = {}

    def print_info(self):
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")        

    
    def calc_geometry(self):
        dx=self.node[1].pos[0]-self.node[0].pos[0]
        dy=self.node[1].pos[1]-self.node[0].pos[1]
        self.length=np.sqrt(dx**2+dy**2)
        self.cs=dx/self.length
        self.sn=dy/self.length

    def calc_temp_force(self,alpha,T0,T1):
        self.alpha=alpha
        if isinstance(T0, float) or isinstance(T0, int):
            self.T0=T0;
            self.temp_force[0]=-self.EA*self.alpha*self.T0
            self.temp_force[3]=self.EA*self.alpha*self.T0
        if isinstance(T1, float) or isinstance(T1, int):
            self.T1=T1
            self.temp_force[2]=self.EI*self.alpha*self.T1
            self.temp_force[5]=-self.EI*self.alpha*self.T1        

    # convert a six-element vector in the local coordinate system to the global coordinate system.
    def l2g_vec6(self,vec):
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=self.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=self.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-self.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
        return np.dot(rotation_matrix.T,vec)    
        
    def g2l_vec6(self,vec):
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=self.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=self.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-self.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
        return np.dot(rotation_matrix,vec)    

    def calc_distributed_force(self,q):
        self.fixedend_force = np.array([0.0, 0.5, self.length/12.0, 0.0, 0.5, -self.length/12.0])*q*self.length
        
    def calc_concentrated_force(self,p):
        self.fixedend_force = np.array([0.0, 0.5, self.length/8.0, 0.0, 0.5, -self.length/8.0])*p
    
    def update_hinge_type(self):
        self.hinge_type=0
        if(self.hinge[0]==True):
            self.hinge_type+=1
        if(self.hinge[1]==True):
            self.hinge_type+=2        

    def l2g_vec2(self,x, y):
        rotation_matrix = np.array([[self.cs,-self.sn],\
                                    [self.sn, self.cs]])
        rotated_points = np.dot(rotation_matrix, np.array([x, y]))
        return rotated_points[0], rotated_points[1]            
