import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from safep.shape_function import ShapeFunction as sf

class Plotter:
    def __init__(self,nodes,elements):
        self.nodes=nodes
        self.elements=elements

        self.support_size=None
        self.hinge_size=None
        self.hinge_shift=False
        self.hinge_offset=None
        self.text_offset=None
        self.deformation_scale=1.0
        self.node_size=None
        self.font_size=20.0
        self.member_width=3.0
        self.member_color="blue"
        self.deformed_member_color="black"
        self.float_format="{:.8g}"
        self.int_format="{:d}"


        self.max_length=0.0
    

    def set_sizes(self):
        self.max_length=0.0
        for e in self.elements:
            self.max_length=max(self.max_length,e.length)
            
        if self.support_size is None:
            self.support_size=0.06*self.max_length
            
        if self.hinge_size is None:
            self.hinge_size=0.03*self.max_length
            
        if self.text_offset is None:
            self.text_offset=0.05*self.max_length
            
        if self.hinge_offset is None:
            self.hinge_offset=0.05*self.max_length

        if self.node_size is None:
            self.node_size=0.01*self.max_length            

    def plot_structure(self,show_node_number,show_element_number):
        self.set_sizes()
        self.fig,self.ax=plt.subplots()
        self.ax.axis('equal')
        if(show_node_number==True):
            self.plot_node_number()
        if(show_element_number==True):
            self.plot_element_number()
        self.plot_member()
        self.plot_node()
        self.plot_support()
        plt.show()
        plt.close()
        
    def plot_deformation(self,deformation_scale):
        self.set_sizes()
        self.fig,self.ax=plt.subplots()
        self.ax.axis('equal')
        self.plot_member()
        self.plot_deformed_member(deformation_scale)
        self.plot_node()
        self.plot_support()
        plt.show()
        plt.close()

    def plot_moment_diagram(self):
        self.set_sizes()
        self.fig,self.ax=plt.subplots()
        self.ax.axis('equal')
        self.plot_member()
        self.plot_bending_moment()
        self.plot_node()
        self.plot_support()
        plt.show()
        plt.close()        

    def plot_member(self):
        for e in self.elements:
            x0=e.node[0].pos[0]
            y0=e.node[0].pos[1]
            x1=e.node[1].pos[0]
            y1=e.node[1].pos[1]
            if isinstance(e.data.get('member_width',None),(int, float)):
                mw=e.data['member_width']
            else:
                mw=self.member_width
            if isinstance(e.data.get('member_color',None),(int, float)):
                col=e.data['member_color']
            else:
                col=self.member_color

            self.ax.plot([x0,x1],[y0,y1],color=col,linewidth=self.member_width)

    def plot_node_number(self):
        offset=[None]*2
        for i,n in enumerate(self.nodes):
            str_num = str(i)
            for j in range(2):
                if(n.offset[j] is not None):
                    offset[j]=n.offset[j]
                else:
                    offset[j]=self.text_offset
            x=n.pos[0]+offset[0]
            y=n.pos[1]+offset[1]
            self.ax.text(x,y,str_num,fontsize=self.font_size)

    def plot_element_number(self):
        offset=[None]*2
        for i,e in enumerate(self.elements):
            x0=e.node[0].pos[0]
            y0=e.node[0].pos[1]
            x1=e.node[1].pos[0]
            y1=e.node[1].pos[1]
            str_num=str(i)
            for j in range(2):
                if(e.offset[j] is not None):
                    offset[j]=e.offset[j]
                else:
                    offset[j]=self.text_offset
            self.ax.text((x0+x1)*0.5+offset[0],\
                         (y0+y1)*0.5+offset[1],str_num,color='blue',fontsize=self.font_size)

    def plot_support(self):
        for n in self.nodes:
            if(n.support.support_type=="FixedEnd"):
                rect = Rectangle((n.pos[0]-self.support_size*0.5,n.pos[1]-self.support_size*0.5), self.support_size, self.support_size, fill=False)
                self.ax.add_patch(rect)
            if(n.support.support_type=="xRollerSupport"):
                self.plot_roller_support(n.pos[0],n.pos[1],0)
            if(n.support.support_type=="yRollerSupport"):
                self.plot_roller_support(n.pos[0],n.pos[1],90)
            if(n.support.support_type=="PinnedSupport"):
                self.plot_triangle(n.pos[0],n.pos[1],0)
        #            if(n.x_fix==True and n.y_fix==True and n.rotation_fix==False):
        #                x=n.x
        #                y=n.y
        #                vertices = [(x, y), (x+0.6*self.support_size,y-self.support_size), (x-0.6*self.support_size,y-self.support_size), (x, y)]
        #                edges = [(0, 1), (1, 2), (2, 0)]
        #                triangle = plt.Polygon(vertices, fill=None, closed=True)

    def rotate(self,x, y, theta_deg):
        theta_rad = np.radians(theta_deg)
        rotation_matrix = np.array([[np.cos(theta_rad), -np.sin(theta_rad)],
                                    [np.sin(theta_rad),  np.cos(theta_rad)]])
        rotated_points = np.dot(rotation_matrix, np.array([x, y]))
        return rotated_points[0], rotated_points[1]
        
    def translate(self,x, y, xt, yt):
        return [xi + xt for xi in x], [yi + yt for yi in y]

    def plot_roller_support(self,x,y,t=0):
        self.plot_bar(x,y,t)
        self.plot_triangle(x,y,t)
    
    
    def plot_bar(self,x, y, t=0.0):
        tx = [0.6*self.support_size, -0.6*self.support_size]
        ty = [-1.4*self.support_size, -1.4*self.support_size]
        x_rotated, y_rotated = self.rotate(tx,ty,t)
        x_final, y_final = self.translate(x_rotated, y_rotated, x, y)
        plt.plot(x_final, y_final,color="black")

    def plot_triangle(self,x, y, t=0.0):
        tx = [0, 0.6*self.support_size, -0.6*self.support_size]
        ty = [0, -self.support_size, -self.support_size]
        x_rotated, y_rotated = self.rotate(tx,ty,t)
        x_final, y_final = self.translate(x_rotated, y_rotated, x, y)
        plt.plot(x_final + [x_final[0]], y_final + [y_final[0]],color="black")

    def plot_node(self):
        for n in self.nodes:
            circle = plt.Circle((n.pos[0],n.pos[1]), self.node_size, color='black', fill=True)
            self.ax.add_artist(circle)        
    
    def calc_def_scale(self):
        self.max_def=0.0
        for n in self.nodes:
            self.max_def=max(self.max_def,np.sqrt(n.x**2+n.y**2))
        if self.max_def==0.0:
            return 0
        return 0.15*self.max_length/self.max_def

    def plot_deformed_member(self,deformation_scale):
        if deformation_scale is None:
            ds=self.deformation_scale
        else:
            ds=deformation_scale
        for e in self.elements:
            n0=e.node[0]
            n1=e.node[1]
            global_disp=np.array([n0.disp[0],n0.disp[1],e.theta[0],\
                                  n1.disp[0],n1.disp[1],e.theta[1]])
            local_disp=e.g2l_vec6(global_disp)
            bending_disp=np.array([local_disp[1],local_disp[2],\
                                   local_disp[4],local_disp[5]])
            print(bending_disp)
            x=np.linspace(0,e.length,21)
            y=np.zeros(21)
            print(len(x),len(y))
            v = np.empty_like(x)
            for i in range(len(x)):
                v[i] = np.dot(bending_disp, sf.cubic_shape_function(x[i], e.length))
            print(v)
            u=(local_disp[3]-local_disp[0])*x/e.length+local_disp[0]
            x+=ds*u
            y+=ds*v
            x_final,y_final=e.l2g_vec2(x,y)
            x_final,y_final=self.translate(x,y,n0.pos[0],n1.pos[1])
            if isinstance(e.data.get('member_width',None),(int, float)):
                mw=e.data['member_width']
            else:
                mw=self.member_width
            if isinstance(e.data.get('deformed_member_color',None),(int, float)):
                col=e.data['deformed_member_color']
            else:
                col=self.deformed_member_color
            self.ax.plot(x_final,y_final,color=col,linewidth=mw)

            
