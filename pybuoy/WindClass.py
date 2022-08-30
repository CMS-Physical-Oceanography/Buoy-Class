import numpy as np
from .OceanBase import Vector2d

"""
===========================================================
-                      REFERENCES
-[1] Large & Pond, Open Ocean Momentum Flux Measurements in 
-    Moderate to Strong Winds, Journal of Physical Oceanography
-    vol. 11, (1981) pp. 324-336.
-[2] V. W. Ekman, On the influence of the Earths
-    rotation on ocean-currents, Arkiv Matematik,
-    Astron. Fysik, 2 (1905) pp. 1-53.
==========================================================="""

class Wind(Vector2d):
    """
    Wind  inputs wind speed as self.i and direction as 
    self.j, or as self.i = u and self.j = v. Like its parent 
    class Vector2d, built in functions in Wind assume self.i 
    is the radial or horizontal component and self.j is the 
    angular or vertical component of the wind vector. To work 
    with typical meterological data angles are taken to be clockwise 
    from the positive y axis.                               
    ATTRIBUTES:
        -self.i = Cartesian u or polar r component of wind
        -self.j = Cartesian v or polar theta component of wind
        -self.cd = drag coefficient. Defult value is 1.15-3
    METHODS:
        -Re-Assigns Attributes: self.new_coordsys(),self.cdnlp().
        -
        -Returns Wind() object: self.windstress(),self.ekmantransport().
    IMPROVEMENTS:
    Add different options for calculating the drag-coefficient and u10. This 
    could be done by adapting the full wind stress function in the Python 
    air-sea toolbox (GitHub link in self.windstress())."""

    def __init__(self,i=None,j=None):
        super().__init__(i,j)
        self.cd = None

#     def initialize(self,cd):
#         """
#         This function is called in self.__init__
#         to assign the drag coefficient based on the 
#         inputted cd. Defult is cd = 1.15e-3."""

#         if self.i[0]==None:
#             return None
#         elif cd==None:
#             return np.ones(len(self.i))*1.15e-3 
#         else: 
#             return cd

    def new_coordsys(self,y_displacement,cart=False):
        """
        This function applies a clockwise rotation 
        y_displacement degrees from the current 
        positive y axis using Vector2d.rot_angles(theta,cart=False).
        INPUTS:
            -self in polar coordinates.
        REASSIGNS:
            -self.j = self.j - y_displacement."""

        self.rot_angles(y_displacement,cart=cart)

    def cdnlp(self,heights):
        """
        MODIFIED FROM:
        https://github.com/pyoceans/python-airsea/blob/master/airsea/windstress.py

        This functions calculates U10 and drag coefficient cd 
        at different heights and windspeeds for the values in 
        self.i as outlined in NOTE Large & Pond 1981. 
        INPUTS:
            -self in polar coordinates.
            -heights = vertical distance from the mesurement height
            -           to the surface roughness(Large & Pond) in meters.
        RE-ASSIGNS: 
            -self.i = u10
            -self.cd = cd"""

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
            cd = (4.9e-4 + 6.5e-5 * u10o)  # cd(u10).
            cd[u10o < 10.15385] = 1.15e-3
            u10 = sp / (1 + a * np.sqrt(cd))  # Next iteration.
            # Keep going until iteration converges.
            ii = np.abs(u10 - u10o) > tol

        self.i = u10
        self.cd = cd

    def windstress(self,cart=False,real=True,airdens = 1.22):
        """
        This function calculates the 10 meter wind stress as 
        WS = drag*airdens*(u10**2). Defult values are drag=1.15e-3,
        and air density = 1.22 kg/m**3.
        INPUTS:
            -self in polar (defult) or Cartesian coordinates. To input
             Cartesian wind coordinates pass with real=False.
        RETURNS:
            -Vector2d object with the self.i,self.j = u,v components of 
             wind stress if real = True (defult). 
             If real = False, this returns an array of complex vectors 
             u+sqrt(-1)*v."""
        
        drag = self.cd

        if cart == True: # if self.i,self.j = u,v
            i = drag*airdens*self.i**2
            j = drag*airdens*self.j**2
        else: # self.i,self.j = r,theta
            rad = np.deg2rad(self.j)
            wndstress = drag*airdens*(self.i**2)
            j = wndstress*np.cos(rad)
            i = wndstress*np.sin(rad)

        if real==False: # if complex output is desired
            out = np.zeros(len(self.i),dtype=np.complex_)
            out.real = i
            out.imag = j
            return out
        else: # for defult Vector2d(i,j) output
            return Vector2d(i,j)

    def ekmantransport(self,lat,swdens = 1025):
        """
        This function calculates the predicted transport from 
        Ekman (1905)[2].
        IMPUTS:
            -self in polar coordinates.
        RETURNS:
            -Vector2d object with self.i,self.j = u transport,v transport."""
        wstrs = self.windstress()
        cor = coriolis(lat)
        Uekx = wstrs.j /(swdens*cor)
        Ueky = wstrs.i / (swdens*cor)       
        
        return Vector2d(Uekx,Ueky)
    
    def lpf_wnd(self,cut_freq,order):

        from .filters import filter_comps
        from scipy.signal import butter
        A,B = butter(order,cut_freq)

        raw_data = np.zeros((len(self.i),2))
        raw_data[:,0] = self.j
        raw_data[:,1] = self.i

        self.i,self.j = filter_comps(raw_data,A,B)