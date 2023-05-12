#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math
import datetime
import sys
from matplotlib.patches import Rectangle
# In[2]:


class Node:
    def __init__(self,x,y,x_force=0.0,y_force=0.0,moment_force=0.0,\
                 x_fix=False,y_fix=False,rotation_fix=False,\
                 u=0.0,v=0.0,theta=0.0,is_hinge=False):
        self.x = x
        self.y = y
        self.x_force=x_force
        self.y_force=y_force
        self.moment_force=moment_force
        self.x_fix=x_fix
        self.y_fix=y_fix
        self.rotation_fix=rotation_fix
        self.u=u
        self.v=v
        self.theta=theta
        
        self.is_hinge=is_hinge
        self.element=[]
        self.num_hinge=0
        
    def print_info(self):
        print(self.x, self.y,self.x_fix,self.y_fix,self.rotation_fix)
        
    def boundary_condition(self,x_fix,y_fix,rotation_fix,u,v,theta):
        if(x_fix is not None):
            self.x_fix = x_fix
        if(y_fix is not None):
            self.y_fix = y_fix
        if(rotation_fix is not None):
            self.rotation_fix = rotation_fix
        if(u is not None):
            if(x_fix==True):
                self.u=u
            else:
                print("x_fix==True, but u is prescribed.")
                sys.exit()
        
        if(v is not None):
            if(y_fix==True):
                self.v=v
            else:
                print("y_fix==True, but v is prescribed.")
                sys.exit()
        
        if(theta is not None):
            if(rotation_fix==True):
                self.theta=theta
            else:
                print("rotation_fix==True, but theta is prescribed.")
                sys.exit()
                
    def load(self,x_force,y_force,moment_force):
        self.x_force = x_force
        self.y_force = y_force
        self.moment_force = moment_force


# In[3]:


class Element:
    def __init__(self,node0,node1,EA,EI,hinge0,hinge1):
        self.node = [node0,node1]
        self.hinge = [hinge0,hinge1]
        self.EA = EA
        self.EI = EI
        self.fixedend_force = np.zeros(6)
        self.temp_force = np.zeros(6)
        
        self.hinge_type=0
        if(hinge0==True):
            self.hinge_type+=1
        if(hinge1==True):
            self.hinge_type+=2
        # 0 : member with no hinges
        # 1 : member with a hinge at its begining
        # 2 : member with a hinge at its end
        # 3 : member with hinges at both ends
    
    def print_info(self):
        print(node,hinge,EA,EI)
        
    def calc_geometry(self,node):
        dx=node[self.node[1]].x-node[self.node[0]].x
        dy=node[self.node[1]].y-node[self.node[0]].y
        self.length=np.sqrt(dx*dx+dy*dy)
        self.cs=dx/self.length
        self.sn=dy/self.length
    
    def calc_temp_force(self,alpha,T0,T1):
        if isinstance(T0, float) or isinstance(T0, int):
            self.T0=T0
            self.temp_force[0]=-self.EA*alpha*T0/self.length
            self.temp_force[3]=self.EA*alpha*T0/self.length
        if isinstance(T1, float) or isinstance(T1, int):
            self.temp_force[2]=self.EI*alpha*T1
            self.temp_force[5]=-self.EI*alpha*T1
        
    def calc_distributed_force(self,q):
        self.fixedend_force = np.array([0.0, 0.5, self.length/12.0, 0.0, 0.5, -self.length/12.0])*q*self.length
        
    def calc_concentrated_force(self,p):
        self.fixedend_force = np.array([0.0, 0.5, self.length/8.0, 0.0, 0.5, -self.length/8.0])*p
    
    def rotate_fixedend_force(self):
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=self.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=self.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-self.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
        return np.dot(rotation_matrix.T,self.fixedend_force)
    
    def rotate_temp_force(self):
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=self.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=self.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-self.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
        return np.dot(rotation_matrix.T,self.temp_force)
    
    def update_hinge_type(self):
        self.hinge_type=0
        if(self.hinge[0]==True):
            self.hinge_type+=1
        if(self.hinge[1]==True):
            self.hinge_type+=2


# In[4]:


class Solve2Dframe:
    def __init__(self,node,element):
        self.node=node
        self.element=element
        
    def calc(self):
        num = 3*len(self.node)
        A=np.zeros((num,num))
        b = np.zeros(num)        
        for i,e in enumerate(self.element):
            local_kmatrix=self.calc_local_kmatrix(e)   
            #print(i,local_kmatrix)
            for iendpoint in range(2):
                inode=e.node[iendpoint]
                for jendpoint in range(2):
                    jnode=e.node[jendpoint]
                    for idir in range(3):
                        for jdir in range(3):
                            A[3*inode+idir,3*jnode+jdir]+=local_kmatrix[3*iendpoint+idir,3*jendpoint+jdir]
            
        for i,n in enumerate(self.node):
            if(n.x_fix==True):
                b[:]-=A[:,3*i]*n.u
                A[:,3*i]=0.0
                A[3*i,3*i]=-1.0
            else:
                b[3*i]+=n.x_force
            if(n.y_fix==True):
                b[:]-=A[:,3*i+1]*n.v
                A[:,3*i+1]=0.0
                A[3*i+1,3*i+1]=-1.0
            else:
                b[3*i+1]+=n.y_force
            if(n.rotation_fix==True):
                b[:]-=A[:,3*i+2]*n.theta
                A[:,3*i+2]=0.0
                A[3*i+2,3*i+2]=-1.0
            else:
                b[3*i+2]+=n.moment_force
                
        for i,e in enumerate(self.element):
            fixedend_force=e.rotate_fixedend_force()  
            temp_force=e.rotate_temp_force()    
            for iendpoint in range(2):
                inode=e.node[iendpoint]
                for idir in range(3):
                    b[3*inode+idir]+=fixedend_force[3*iendpoint+idir]+temp_force[3*iendpoint+idir]   
        #print(b)
        #print(A)
        x = np.linalg.solve(A, b)
        
        for i,n in enumerate(self.node):
            if(n.x_fix==True):
                n.x_force=x[3*i]
            else:
                n.u=x[3*i]
            if(n.y_fix==True):
                n.y_force=x[3*i+1]
            else:
                n.v=x[3*i+1]
            if(n.rotation_fix==True):
                n.moment_force=x[3*i+2]
            else:
                n.theta=x[3*i+2]
                
    def print_result(self):
        print("#node : x-disp, y-disp, angle")
        for i,n in enumerate(self.node):
            print(i,":",n.u,n.v,n.theta)
        print("#node : x-force, y-force, moment")
        for i,n in enumerate(self.node):
            print(i,":",n.x_force,n.y_force,n.moment_force)
        
    def calc_local_kmatrix(self,e):
        local_kmatrix=np.zeros((6,6))
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
            f=e.fixedend_force
            e.fixedend_force[1]=f[1]-1.5*f[5]/e.length
            e.fixedend_force[2]=f[2]-0.5*f[5]
            e.fixedend_force[4]=f[4]+1.5*f[2]/e.length
            e.fixedend_force[5]=0.0
        if(e.hinge_type==3):
            f=e.fixedend_force
            e.fixedend_force[1]=f[1]-(f[2]+f[5])/e.length
            e.fixedend_force[2]=0.0
            e.fixedend_force[4]=f[4]+(f[2]+f[5])/e.length
            e.fixedend_force[5]=0.0
            
        rotation_matrix=np.zeros((6,6))
        rotation_matrix[0,0]=rotation_matrix[1,1]=rotation_matrix[3,3]=rotation_matrix[4,4]=e.cs
        rotation_matrix[0,1]=rotation_matrix[3,4]=e.sn
        rotation_matrix[1,0]=rotation_matrix[4,3]=-e.sn
        rotation_matrix[2,2]=rotation_matrix[5,5]=1.0
           
        return np.dot(np.dot(rotation_matrix.T,local_kmatrix),rotation_matrix)
    


# In[5]:


def shape_function(x,L):
    return np.array([2.0*(x/L)**3-3.0*(x/L)**2+1,\
                     L*((x/L)**3-2.0*(x/L)**2+(x/L)),\
                     -2.0*(x/L)**3+3.0*(x/L)**2,\
                     L*((x/L)**3-(x/L)**2)])


# In[6]:


class DrawResult:
    def __init__(self,node,element):
        self.node=node
        self.element=element
        self.support_size=-1.0
        self.hinge_size=-1.0
        self.hinge_loc=0.95
        self.hinge_shift=False
        self.text_size=3.0
        self.shift_size=-1.0
        self.member_width=3.0
        self.cal_size_done=False
    
    def set_sizes(self):
        self.max_length=self.element[0].length
        for e in self.element:
            self.max_length=max(self.max_length,e.length)
        if self.support_size<0.0:
            self.support_size=0.06*self.max_length
        if self.hinge_size<0.0:
            self.hinge_size=0.03*self.max_length
        if self.shift_size<0:
            self.shift_size=0.01*self.max_length
       
    def structure(self):
        self.set_sizes()
 #       plt.xticks([])
 #       plt.yticks([])
        self.fig, self.ax = plt.subplots()
        self.ax.axis('equal')
        self.plot_node()
        self.plot_element()
        self.plot_node_number()
        self.plot_element_number()
        self.plot_hinge()
        self.plot_support()
        #self.plot_force()
        plt.show()
        plt.close()

    def calc_def_scale(self):
        self.max_def=np.sqrt(self.node[0].u**2+self.node[1]**2)
        for n in self.node:
            self.max_def=max(self.max_def,np.sqrt(n.x**2+n.y**2))
        return 0.15*self.max_length/self.max_def
        
    def deformation(self,def_scale):
        self.def_scale=def_scale
        self.set_sizes()
        if(def_scale<0):
            self.def_scale=self.calc_def_scale()
 #       plt.xticks([])
 #       plt.yticks([])
        self.fig, self.ax = plt.subplots()
        self.ax.axis('equal')
        self.plot_node()
        self.plot_element()
        self.plot_deformed_member()
        self.plot_hinge()
        self.plot_support()
        #self.plot_force()
        plt.show()
        plt.close()

    def axial_force(self):
        self.set_sizes()
        self.fig, self.ax = plt.subplots()
        self.ax.axis('equal')
        self.plot_node()
        self.plot_element()
        self.plot_axial_force()
        self.plot_support()
        
    def plot_axial_force(self):
        for i,e in enumerate(self.element):
            local_u0=e.cs*self.node[e.node[0]].u+e.sn*self.node[e.node[0]].v
            local_u1=e.cs*self.node[e.node[1]].u+e.sn*self.node[e.node[1]].v
            N=e.EA*(local_u1-local_u0)/e.length
            str_N = "{:.8g}".format(N)
            x=(self.node[e.node[0]].x+self.node[e.node[1]].x)*0.5
            y=(self.node[e.node[0]].y+self.node[e.node[1]].y)*0.5
            self.ax.text(x,y,str_N)
    def plot_node(self):
        for i,n in enumerate(self.node):
            circle = plt.Circle((n.x,n.y), 0.2*self.support_size, color='black', fill=True)
            self.ax.add_artist(circle)
    
    def plot_node_number(self):
        for i,n in enumerate(self.node):
            str_num = str(i)
            x=n.x+self.shift_size
            y=n.y+self.shift_size
            self.ax.text(x,y,str_num)
    
    def plot_element(self):
        for i,e in enumerate(self.element):
            x0=self.node[e.node[0]].x
            y0=self.node[e.node[0]].y
            x1=self.node[e.node[1]].x
            y1=self.node[e.node[1]].y
            self.ax.plot([x0,x1],[y0,y1],"b-",linewidth=self.member_width)
    
    def plot_element_number(self):
        for i,e in enumerate(self.element):
            x0=self.node[e.node[0]].x
            y0=self.node[e.node[0]].y
            x1=self.node[e.node[1]].x
            y1=self.node[e.node[1]].y
            str_num=str(i)
            self.ax.text((x0+x1)*0.5+self.shift_size,\
                         (y0+y1)*0.5+self.shift_size,str_num,color='blue')
        
    def plot_deformed_member(self):
        for i,e in enumerate(self.element):
            x0=self.node[e.node[0]].x+self.def_scale*self.node[e.node[0]].u
            y0=self.node[e.node[0]].y+self.def_scale*self.node[e.node[0]].v
            x1=self.node[e.node[1]].x+self.def_scale*self.node[e.node[1]].u
            y1=self.node[e.node[1]].y+self.def_scale*self.node[e.node[1]].v
            self.ax.plot([x0,x1],[y0,y1],"k-",linewidth=2)
    
    def plot_hinge(self):
        if(self.hinge_shift==True):
            for i,e in enumerate(self.element):
                node0=e.node[0]
                node1=e.node[1]
                hinge_scale=0.3
                if(e.hinge[0]==True):
                    xh=self.loc_hinge*self.node[node0].x+(1.0-self.loc_hinge)*self.node[node1].x
                    yh=self.loc_hinge*self.node[node0].y+(1.0-self.loc_hinge)*self.node[node1].y
                    circle = plt.Circle((xh,yh), self.hinge_size, color='black', fill=False)
                    self.ax.add_artist(circle)
                #circle = plt.Circle((xh,yh), 0.04, color='white', fill=True)
                #ax.add_artist(circle)
                if(e.hinge[1]==True):
                    xh=(1.0-self.loc_hinge)*self.node[node0].x+self.loc_hinge*self.node[node1].x
                    yh=(1.0-self.loc_hinge)*self.node[node0].y+self.loc_hinge*self.node[node1].y
                    circle = plt.Circle((xh,yh), self.hinge_size, color='black', fill=False)
                    self.ax.add_artist(circle)
                    #circle = plt.Circle((xh,yh), 0.04, color='white', fill=True)
                    #ax.add_artist(circle)
        if(self.hinge_shift==False):
            for i,e in enumerate(self.element):
                node0=e.node[0]
                node1=e.node[1]
                hinge_scale=0.3
                if(self.node[node0].is_hinge==False and e.hinge[0]==True):
                    xh=self.loc_hinge*self.node[node0].x+(1.0-self.loc_hinge)*self.node[node1].x
                    yh=self.loc_hinge*self.node[node0].y+(1.0-self.loc_hinge)*self.node[node1].y
                    circle = plt.Circle((xh,yh), self.hinge_size, color='black', fill=False)
                    self.ax.add_artist(circle)
                if(self.node[node1].is_hinge==False and e.hinge[1]==True):
                    xh=(1.0-self.loc_hinge)*self.node[node0].x+self.loc_hinge*self.node[node1].x
                    yh=(1.0-self.loc_hinge)*self.node[node0].y+self.loc_hinge*self.node[node1].y
                    circle = plt.Circle((xh,yh), self.hinge_size, color='black', fill=False)
                    self.ax.add_artist(circle)
            for i,n in enumerate(self.node):
                if(n.is_hinge==True):
                    circle = plt.Circle((n.x,n.y),self.hinge_size, color='black',fill=False)
                    self.ax.add_artist(circle)
                    

    def plot_support(self):
        for i,n in enumerate(self.node):
            if(n.x_fix==True and n.y_fix==True and n.rotation_fix==True):
                rect = Rectangle((n.x-self.support_size*0.5,n.y-self.support_size*0.5), self.support_size, self.support_size, fill=False)
                self.ax.add_patch(rect)
            if(n.x_fix==True and n.y_fix==True and n.rotation_fix==False):
                x=n.x
                y=n.y
                vertices = [(x, y), (x+0.6*self.support_size,y-self.support_size), (x-0.6*self.support_size,y-self.support_size), (x, y)]
                edges = [(0, 1), (1, 2), (2, 0)]
                triangle = plt.Polygon(vertices, fill=None, closed=True)
                self.ax.add_patch(triangle)
            if(n.x_fix==False and n.y_fix==True and n.rotation_fix==False):
                x=n.x
                y=n.y
                vertices = [(x, y), (x+0.6*self.support_size,y-self.support_size), (x-0.6*self.support_size,y-self.support_size), (x, y)]
                triangle = plt.Polygon(vertices, fill=None, closed=True)
                self.ax.plot([x+0.6*self.support_size,x-0.6*self.support_size],[y-1.2*self.support_size,y-1.2*self.support_size])
                self.ax.add_patch(triangle)
            if(n.x_fix==True and n.y_fix==False and n.rotation_fix==False):
                x=n.x
                y=n.y
                vertices = [(x, y), (x-self.support_size,y-0.6*self.support_size), (x-self.support_size,y+0.6*self.support_size), (x, y)]
                triangle = plt.Polygon(vertices, fill=None, closed=True)
                self.ax.plot([x-1.2*self.support_size,x-1.2*self.support_size],[y-0.6*self.support_size,y+0.6*self.support_size])
                self.ax.add_patch(triangle)
                


# In[7]:


class StructuralAnalysis:
    def __init__(self, EA=1, EI=1):
        self.EA = EA
        self.EI = EI
        self.element = []
        self.node = []
        self.sol = Solve2Dframe(self.node,self.element)
        self.drw = DrawResult(self.node, self.element)
    
    def make_node(self,x,y):
        
        print("node[",len(self.node),"] = [",x,",",y,"]")
        self.node.append(Node(x,y));
        
    def make_beam(self,node0,node1,EA=None, EI=None, hinge0=False,hinge1=False):
        if EA is None:
            EA = self.EA
        if EI is None:
            EI = self.EI
        print("element[",len(self.element),"] = [",node0,",",node1,"] : beam")
        self.element.append(Element(node0,node1,EA,EI,hinge0,hinge1))
        self.element[-1].calc_geometry(self.node)
    
    def make_truss_member(self,node0,node1,EA=None):
        if EA is None:
            EA = self.EA
        EI=self.EI
        hinge0=True
        hinge1=True
        print("element[",len(self.element),"] = [",node0,",",node1,"] : truss member")
        self.element.append(Element(node0,node1,EA,EI,hinge0,hinge1))
        self.element[-1].calc_geometry(self.node)
        
    def set_boundary_condition(self,node_num,x_fix=None,y_fix=None,
                               rotation_fix=None,u=None,v=None,theta=None):
        if(len(self.node)<=node_num):
            print("set_boundary_condtion error")
            print("Check node number.")
            sys.exit()
        self.node[node_num].boundary_condition(x_fix,y_fix,rotation_fix,u,v,theta)
        
    def calc_element_force(self,element_num,force_type,force):
        if(force_type=="concentrated"):
            self.element[element_num].calc_concentrated_force(force)
        elif(force_type=="distributed"):
            self.element[element_num].calc_distributed_force(force)
        else:
            print("calc_element_force error")
            print(force_type,"is unknown force type.")
            sys.exit()
            
    def apply_temperature(self,element_num,alpha,T0=None,T1=None):
        self.element[element_num].calc_temp_force(alpha,T0,T1)                
    
    def apply_load(self,node_num,x_force=0,y_force=0,moment_force=0):
        self.node[node_num].load(x_force,y_force,moment_force)
        
    def print_node(self):
        for n in self.node:
            n.print_info()
    
    def print_element(self):
        for e in self.element:
            e.print_info()
            
    def solve(self):
        self.update_hinge()
        self.sol.calc()
        
    def draw_structure(self):
        self.update_hinge()
        self.drw.structure()
    
    def draw_deformation(self,def_scale=-1.0):
        self.drw.deformation(def_scale)
        
    def draw_axial_force(self):
        self.drw.axial_force()
        
    def print_result(self):
        print("#node : x-disp, y-disp, angle")
        for i,n in enumerate(self.node):
            print(i,":",n.u,n.v,n.theta)
        print("#node : x-force, y-force, moment")
        for i,n in enumerate(self.node):
            print(i,":",n.x_force,n.y_force,n.moment_force)
    
    def printnow(self):
        print(datetime.datetime.now())
    
    def convert_hinge(self,node_num):
        self.node[node_num].is_hinge=True

    def update_hinge(self):
        for i,e in enumerate(self.element):
            node0=self.node[e.node[0]]
            if(node0.is_hinge==True):
                e.hinge[0]=True
            node1=self.node[e.node[1]]
            if(node1.is_hinge==True):
                e.hinge[1]=True
            node0.element.append(i)
            node1.element.append(i)
            if(e.hinge[0]==True):
                node0.num_hinge+=1
            if(e.hinge[1]==True):
                node1.num_hinge+=1
        for i,n in enumerate(self.node):
            if len(n.element)==n.num_hinge:
                n.is_hinge=True
        for i,n in enumerate(self.node):
            print(i,n.is_hinge,n.element)
            if n.is_hinge==True:
                for element_num in n.element:
                    e=self.element[element_num]
                    if(e.node[0]==i):
                        e.hinge[0]=False
                        break
                    elif(e.node[1]==i):
                        e.hinge[1]=False
                        break
        for i,e in enumerate(self.element):
            e.update_hinge_type()


# In[8]:

