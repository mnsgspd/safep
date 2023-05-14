class ShapeFunction:
    def shape_function(x,L):
        return np.array([2.0*(x/L)**3-3.0*(x/L)**2+1,\
                         L*((x/L)**3-2.0*(x/L)**2+(x/L)),\
                         -2.0*(x/L)**3+3.0*(x/L)**2,\
                         L*((x/L)**3-(x/L)**2)])
