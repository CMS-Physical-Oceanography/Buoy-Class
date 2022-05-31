import numpy as np
from .OceanBase import Vector2d
"""
===========================================================
-                      REFERENCES
-[1] M. S. Longuet-Higgins, Mass transport in water waves, 
-    Philos. Trans. Roy. Soc. London, A245, (1953) pp. 535-581.
-[2] M. Fewings, Observations of Cross-Shelf Flow
-    Driven by Cross-Shelf Winds on the Inner Continental Shelf, 
-    Journal Of Physical Oceanography, 38 (2008) pp. 2358-2378. 
==========================================================="""

class Waves(Vector2d):
    """
    Waves takes the phase speed as self.i and direction as 
    self.j, or i,j = u,v components of the bulk phase speed. Angles 
    are assumed to be the clockwise displacement from the positive 
    y-axis in degrees.                          
    ATTRIBUTES:
        -self.i = Cartesian u or polar r component of phase speed.
        -self.j = Cartesian v or polar theta component of phase speed.
        -self.swh = Significant wave height.
        -self.T = Wave period.
        -self.k = Wave number.
        -self.cg = Bulk group speed. 
        -self.depth = water depth.
    METHODS:
        -Re-Assigns Attributes: self.new_coordsys(),self.getk(),self.getc(),self.getcg().
        -
        -Returns Vector2d object: self.stokesdrift(),self.energy(),self.energyflux(),
        -                         self.transport(),self.waveforcing().
    IMPROVEMENTS:
    Add support for the full wave spectrum and include more solutions to the 
    Stokes-Drift."""

    def __init__(self,i=None,j=None,swh=None, k=None,depth=None,cg=None,T=None):
        super().__init__(i,j)
        self.swh,self.k,self.T,self.cg,self.depth = swh,k,T,cg,depth

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
     
    def getc(self):
        omega = (2*np.pi)/np.array(self.T)
        self.i = omega/self.k

    def getcg(self):

        Cg_ = lambda c,k,h :(  c*.5 * (1 + ((2*k*h)/np.sinh(2*k*h)))  )
        self.cg = Cg_(self.i,self.k,self.h)

    def getk(self):
        """
        Uses the Newton-Raphson method to calculate the wavenumber
        k with an inital guess k0 being for shallow water waves where
        k0 = omega/sqrt(gh) and omega is the angular frequency 
        omega = 2pi/T. The function f is the dispersion relation for 
        surface gravity waves 0 = -omega + sqrt(gk*tanh(kh))
        This method is outlined in: 
        "A close approximation of wave dispersion relation for direct 
        calculation of wavelength in any coastal water depth" -Zai-Jin You

        This function does not return anything, instead it just assigns
        the calculated value to self.k."""
        
        # create omega and the initial guess of k defined in the doc string
        omega = (np.pi*2)/np.array(self.T)
        
        k = omega/np.sqrt(9.81*self.depth) # make sure depth is set if error is thrown here 
    
        # create initial values using k0
        f = (9.81 * k * np.tanh(k*self.depth)) - omega**2 
        dfdk = 9.81*self.depth*k*( 1/(np.cosh(k*self.depth)**2) ) + (9.81*np.tanh(k*self.depth))
        ii = np.any(abs(f) > 1e-10)
        
        while ii: # while ii == True
                                                                  
            k = k - (f/dfdk)
            f = (9.81*k*np.tanh(k*self.depth)) - omega**2
            dfdk = 9.81*self.depth*k*( 1/(np.cosh(k*self.depth)**2) ) + (9.81 * np.tanh(k*self.depth))
            ii = np.any(abs(f) > 1e-10)
            
        self.k = k
     
    def stokesdrift(self,res=None):

        if res  == None:
            res = int(self.depth)

        lendata = len(self.i)
        # initialize return and depth arrays
        xshelf = np.zeros((res,lendata))
        yshelf = np.zeros((res,lendata))
        z = np.linspace(-0.05,-self.depth,res)

        # data slices used in calculations
        c = self.i # calculated with self.getc()
        theta = np.deg2rad(self.j) # 2/7/22: theta is degrees clockwise of +y
        swh = self.swh
        k_ = self.k # calculated with self.getk()
            
        for i in range(lendata): # for each data point in [start:stop]

            k = k_[i] 
            A =   ( (( 9.81*k*swh[i]**2) / ( 8*c[i] )) )
            # Cross-Shelf Stokes Drift
            xshelf[:,i] =( (A/ np.sinh( 2*k*self.depth)) 
                                * np.cosh( 2*k*(z+self.depth) ))*np.sin(theta[i]) 
            # Along-Shelf Stokes Drift
            yshelf[:,i] =( (A/ np.sinh( 2*k*self.depth)) 
                                * np.cosh( 2*k*(z+self.depth) )*np.cos(theta[i]) )
                
        return Vector2d(xshelf,yshelf) # returns a vector of stokes drift components

    def energyflux(self,swdens=1025):
        h = self.depth
        Eflux = ( (((swdens/16)*9.81*swh**2*self.i)*
                (1 +((2* self.k*h) / (np.sinh(2*self.k*h))))) )

        rad = np.deg2rad(self.j)
    
        i = Eflux*np.sin(theta)
        j = Efluz*np.cos(theta)

        return Vector2d(i,j)

    def energy(self,swdens=1025):
        return (swdens/8) * 9.81 * swh**2 

    def transport(self):

        theta = np.deg2rad(self.j)
        Qw = (9.81*self.swh**2)/(16*self.i)

        i = Q * np.sin(theta)
        j = Q * np.cos(theta)

        return Vector2d(i,j)

    def waveforcing(self):
        """
        This function returns a measure of "wave forcing"
        used in [3] to bin current observations.
        INPUTS:
            -self in polar coordinates.
        RETURNS:
            -Vector2d object with self.i,self,j =
            -     Hsig^2*sin(theta),Hsig^2*cos(theta)"""     

        Hsig = self.swh
        traj = np.deg2rad(self.j)
        
        return(Vector2d((Hsig**2)*np.sin(traj),(Hsig**2)*np.cos(traj)))