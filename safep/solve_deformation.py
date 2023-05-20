import numpy as np

class SolveDeformation:
    def __init__(self,nodes,elements):
        self.nodes=nodes
        self.elements=elements
        
    def solve(self):
        num = 3*len(self.nodes)
        A=np.zeros((num,num))
        b = np.zeros(num)        
        for e in self.elements:
            local_kmatrix=self.calc_local_kmatrix(e)   
            for iendpoint in range(2):
                inode=e.node[iendpoint].node_num
                for jendpoint in range(2):
                    jnode=e.node[jendpoint].node_num
                    for idir in range(3):
                        for jdir in range(3):
                            A[3*inode+idir,3*jnode+jdir]+=local_kmatrix[3*iendpoint+idir,3*jendpoint+jdir]

        for i,n in enumerate(self.nodes):
            if(n.is_all_hinge==True):
                A[3*i+2,:]=0.0
                A[:,3*i+2]=0.0
                A[3*i+2,3*i+2]=1.0
        for i,n in enumerate(self.nodes):
            if(n.disp_specified[0]==True):
                b[:]-=A[:,3*i]*n.disp[0]
                A[:,3*i]=0.0
                A[3*i,3*i]=-1.0
            else:
                b[3*i]+=n.force[0]
            if(n.disp_specified[1]==True):
                b[:]-=A[:,3*i+1]*n.disp[1]
                A[:,3*i+1]=0.0
                A[3*i+1,3*i+1]=-1.0
            else:
                b[3*i+1]+=n.force[1]
            if(n.disp_specified[2]==True):
                b[:]-=A[:,3*i+2]*n.disp[2]
                A[:,3*i+2]=0.0
                A[3*i+2,3*i+2]=-1.0
            else:
                b[3*i+2]+=n.force[2]
        for e in self.elements:
            print(e.fixedend_force)
            fixedend_force=e.l2g_vec6(e.fixedend_force)
            temp_force=e.l2g_vec6(e.temp_force)
            for iendpoint in range(2):
                inode=e.node[iendpoint].node_num
                for idir in range(3):
                    b[3*inode+idir]+=fixedend_force[3*iendpoint+idir]\
                        +temp_force[3*iendpoint+idir]   

        print(A,b)
                    
        x = np.linalg.solve(A, b)
        
        for i,n in enumerate(self.nodes):
            if(n.disp_specified[0]==True):
                n.force[0]=x[3*i]
            else:
                n.disp[0]=x[3*i]
            if(n.disp_specified[1]==True):
                n.force[1]=x[3*i+1]
            else:
                n.disp[1]=x[3*i+1]
            if(n.disp_specified[2]==True):
                n.force[2]=x[3*i+2]
            else:
                n.disp[2]=x[3*i+2]

        for e in self.elements:
            n0=e.node[0]
            n1=e.node[1]
            global_disp=np.array([n0.disp[0],n0.disp[1],n0.disp[2],\
                                  n1.disp[0],n1.disp[1],n1.disp[2]])
            local_disp=e.l2g_vec6(global_disp)
            if e.hinge_type==0:
                e.theta[0]=local_disp[2]
                e.theta[1]=local_disp[5]
            elif e.hinge_type==1:
                e.theta[0]=1.5/e.length*(-local_disp[1]+local_disp[4])\
                    -0.5*local_disp[5]+0.25*e.length/e.EI+e.fixedend_force[2]
                e.theta[1]=local_disp[5]
            elif e.hinge_type==2:
                e.theta[0]=e.n0.disp[2]
                e.theta[1]=1.5/e.length*(-local_disp[1]+local_disp[4])\
                    -0.5*local_disp[2]+0.25*e.length/e.EI*e.fixedend_force[5]
            elif e.hinge_type==3:
                e.theta[0]=(-local_disp[1]+local_disp[4])/e.length\
                    +e.length/(6.0*e.EI)\
                    *(2*e.fixedend_force[2]-e.fixedend_force[5])
                e.theta[1]=(-local_disp[1]+local_disp[4])/e.length\
                    +e.length/(6.0*e.EI)\
                    *(2*e.fixedend_force[5]-e.fixedend_force[2])
            
        self.print_result()
        
    def print_result(self):
        print("#node : x-disp, y-disp, angle")
        for i,n in enumerate(self.nodes):
            print(i,":",n.disp)
        print("#node : x-force, y-force, moment")
        for i,n in enumerate(self.nodes):
            print(i,":",n.force)

    def calc_local_kmatrix(self,e):
        local_kmatrix=np.zeros((6,6))

        
        print(e.EA,e.length,"hoge")

        
        local_kmatrix[0,0]=local_kmatrix[3,3]=e.EA/e.length
        local_kmatrix[0,3]=local_kmatrix[3,0]=-e.EA/e.length
        if(e.hinge_type==0):
            local_kmatrix[1,1]=local_kmatrix[4,4]=12.0*e.EI/e.length**3
            local_kmatrix[1,4]=local_kmatrix[4,1]=-12.0*e.EI/e.length**3
            local_kmatrix[1,2]=local_kmatrix[2,1]=local_kmatrix[1,5]=local_kmatrix[5,1]=6.0*e.EI/e.length**2
            local_kmatrix[2,4]=local_kmatrix[4,2]=local_kmatrix[4,5]=local_kmatrix[5,4]=-6.0*e.EI/e.length**2
            local_kmatrix[2,2]=local_kmatrix[5,5]=4.0*e.EI/e.length
            local_kmatrix[2,5]=local_kmatrix[5,2]=2.0*e.EI/e.length
        if(e.hinge_type==1):
            local_kmatrix[1,1]=local_kmatrix[4,4]=3.0*e.EI/e.length**3
            local_kmatrix[1,4]=local_kmatrix[4,1]=-3.0*e.EI/e.length**3
            local_kmatrix[1,5]=local_kmatrix[5,1]=3.0*e.EI/e.length**2
            local_kmatrix[4,5]=local_kmatrix[5,4]=-3.0*e.EI/e.length**2
            local_kmatrix[5,5]=3.0*e.EI/e.length
            f=e.fixedend_force
            e.fixedend_force[1]=f[1]-1.5*f[2]/e.length
            e.fixedend_force[2]=0.0
            e.fixedend_force[4]=f[4]+1.5*f[2]/e.length
            e.fixedend_force[5]=f[5]-0.5*f[2]
        if(e.hinge_type==2):
            local_kmatrix[1,1]=local_kmatrix[4,4]=3.0*e.EI/e.length**3
            local_kmatrix[1,4]=local_kmatrix[4,1]=-3.0*e.EI/e.length**3
            local_kmatrix[1,2]=local_kmatrix[2,1]=3.0*e.EI/e.length**2
            local_kmatrix[2,4]=local_kmatrix[4,2]=-3.0*e.EI/e.length**2
            local_kmatrix[2,2]=3.0*e.EI/e.length
            f=e.fixedend_force.copy()
            e.fixedend_force[1]=f[1]-1.5*f[5]/e.length
            e.fixedend_force[2]=f[2]-0.5*f[5]
            e.fixedend_force[4]=f[4]+1.5*f[2]/e.length
            e.fixedend_force[5]=0.0
        if(e.hinge_type==3):
            f=e.fixedend_force.copy()
            print(e.fixedend_force)
            e.fixedend_force[1]=f[1]-(f[2]+f[5])/e.length
            e.fixedend_force[2]=0.0
            e.fixedend_force[4]=f[4]+(f[2]+f[5])/e.length
            e.fixedend_force[5]=0.0
            print(e.fixedend_force,"hoge")
            
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=e.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=e.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-e.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
           
        return np.dot(np.dot(rotation_matrix.T,local_kmatrix),rotation_matrix)

