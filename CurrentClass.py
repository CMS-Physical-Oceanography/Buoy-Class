import numpy as np
from VectorClass import Vector2d

class Currents:
    
    """
    Currents Class expects input data to have dimentions depth by time
    or depth by spatial dim. Angles are reported from the positive y axis.
    """
    
    def __init__(self,csflow,asflow,depths_ = np.nan):
        
        self.i = csflow
        self.j = asflow
        self.z = depths_
        
        
    def rot_angles(self,theta_prime):
        
        
        
        u_0 = self.i
        v_0 = self.j
        
        
        # calculates radial comp
        r_0 = np.sqrt(u_0**2 + v_0**2)
        # calculates angle from the positive y
        theta_0 = np.arctan2(u_0,v_0)
    
        # maps angles to new axis
        theta_f = theta_0 - np.deg2rad(theta_prime)
       
        # calculates components for new coordinate system
        self.i = r_0 * np.sin(theta_f) # cross-shelf
        self.j = r_0 * np.cos(theta_f) # along-shelf
        


    def transport(self):



        tx = range(len(self.i[0]))
        ty = range(len(self.j[0]))
        
        cross_shelf = lambda t : np.nansum(self.i[:,t])
        along_shelf = lambda t : np.nansum(self.j[:,t])

        cs_transports = np.array(list(map(cross_shelf,tx)))
        as_transports = np.array(list(map(along_shelf,ty)))

        return Vector2d(cs_transports,as_transports)



    def depth_average(self,real=True):



        # tx = range(len(self.i[0]))
        # ty = range(len(self.j[0]))
        
        # cross_shelf = lambda t : np.nanmean(self.i[:,t])
        # along_shelf = lambda t : np.nanmean(self.j[:,t])

        # cs_da = np.array(list(map(cross_shelf,tx)))
        # as_da = np.array(list(map(along_shelf,ty)))

        cs_da = np.nanmean(self.i,axis=0)
        as_da = np.nanmean(self.j,axis=0)


        if real == False:

            out = np.zeros(self.i.shape[1],dtype=np.complex_)
            out.real = cs_da
            out.imag = as_da
            return out 

        return Vector2d(cs_da,as_da) 

# class Ekman_Spiral:
#     def __init__(self,A,delta):
#         self.u = list(map(lambda z: A*np.exp(z/delta)*
#                      np.sin((np.pi/4)-(z/delta)),range(-h,1)))
#         self.v = list(map(lambda z: A*np.exp(z/delta)*
#                      np.cos((np.pi/4)-(z/delta)),range(-h,1)))