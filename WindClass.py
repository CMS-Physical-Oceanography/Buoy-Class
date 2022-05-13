import numpy as np
from VectorClass import Vector2d


class Wind(Vector2d):


    """
    Wind  inputs wind speed as self.i and direction as 
    self.j
 
    Like its parent class Vector2d, built in functions in 
    Wind assume self.i is the radial or horizontal 
    component and self.j is the angular or vertical component
    of the wind vector. To work with typical meterological data
    angles are taken to be clockwise from the positive y axis
    
    
                                                            
    __init__():

    self.i
    self.j
    self.cd # drag coeff.
    
    1) windstress(self, drag = 1.4E-3, airdens = 1.22, 
              return_vector = True,zero_axis = 'y',input_coordinates = 'polar)
              
       self.windstress() assumes that that inputed windspeeds and drag coefficients 
       are adjusted to the 10 meter wind speed equivalents. If data is collected at 
       a different height call self.cdnlp before

                                                               """

    def __init__(self,i,j):
        
        super().__init__(i,j)
        self.cd = np.ones(len(self.i))*1.15e-3 


    def cdnlp(self,heights):

        """
        =======================================================
        MODIFIED FROM:
        https://github.com/pyoceans/python-airsea/blob/master/airsea/windstress.py
        
        The python-airsea toolbox was translated from the 
        matlab version.  

        This functions calculates U10 and drag coefficient cd 
        at different heights and windspeeds for the values in 
        self.i as outlined in Large & Pond 1981. 

        inputs:
        
        self

        heights = vertical distance from the mesurement height
                  to the surface roughness(Large & Pond) in meters.

        returns: Nothing

        re-assigns: 
        
        self.i = u10
        
        self.cd = cd

        =======================================================
        """
        sp = self.i
        
        tol = 0.001  # minimum convergence 
                  
           
        z = np.array(heights)
        a = np.log(z / 10.) / .4  # Log-layer correction factor.
        u10o = np.zeros(sp.shape)
        cd = 1.15e-3 * np.ones(sp.shape)
        u10 =  sp / (1 + a * np.sqrt(cd))
        ii = np.abs(u10 - u10o) > tol
        
        while np.any(ii):
            u10o = u10
            cd = (4.9e-4 + 6.5e-5 * u10o)  # Compute cd(u10).
            cd[u10o < 10.15385] = 1.15e-3
            u10 = sp / (1 + a * np.sqrt(cd))  # Next iteration.
            # Keep going until iteration converges.
            ii = np.abs(u10 - u10o) >tol

        self.i = u10
        self.cd = cd
        

    def Ekmantransport(self,swdens = 1023,lat = 33.5):

        wstrs = self.windstress(return_vector= True)
        
        cor = coriolis(lat)
        Uekx = wstrs.j /(swdens*cor)
        Ueky = wstrs.i / (swdens*cor)

        return Vector2d(Uekx,Ueky)



    def windstress(self, airdens = 1.22,real=True):
        
        """

        MAKE SURE TO CALL SELF.CDNLP OR DRAG
        WILL BE SET TO A DEFULT VALUE OF 
        1.15e-3.
        """
 
        drag = self.cd
       
        rad = np.deg2rad(self.j)
        wndstress = drag*airdens*self.i**2
        j = wndstress*np.cos(rad)
        i = wndstress*np.sin(rad)

        if real==False:

            out = np.zeros(len(self.i),dtype=np.complex_)
            out.real = i
            out.imag = j
            return out

        else:

            return Vector2d(i,j)
