import numpy as np

"""
Contains:
    -class: Vector2d
    -class: VectorPlots"""

class Vector2d:
    """
    DESCRIPTION:
    This class holds vector components self.i and self.j as
    1d numpy arrays and contains methods for basic rotations 
    and operations. Vector2d is the parent class of Wind(), 
    Waves(), and Currents().
    NOTE 
    -   ALL FUNCTIONS TAKE ANGLES TO BE REPORTED AS THE ANGLE 
    -   CLOCKWISE FROM THE POSITIVE Y-AXIS IN DEGREES FOR 
    -   WORKING WITH METEROLOGICAL DATA. 
    -       theta
    -       |^ /
    -    y  | /
    -       |/_________ x
    NOTE
    ATTRIBUTES:
        -self.i = cartesian x component or polar radial component
        -
        -self.j = cartesian y component or polar angular component
    METHODS:
        -Re-Assigns attribute values: self.rot_angles(theta), self.invert()
        -
        -Returns Vector2d objects: self.adjust(), self.mag(), self.polar2cart(),
        -                          self.cart2polar(), self.polar_rot(theta) 
    IMPROVEMENTS:
    It might be useful to be able to choose between degrees and radians and  which 
    axis angles are taken from.
    ============================================================================"""       
    
    def __init__(self,i,j):
        self.i = np.array(i,dtype=float)
        self.j = np.array(j,dtype=float)

    def adjust(self):
        """
        This function adds 360 degrees to 
        negative angles to keep values in 
        [0,360]. This is called consistently
        in Vector2d and should be made optional
        in the future.
        INPUTS:
            -self
        REASSIGNS:
            -self"""

        angles = self.j
        angles[angles<0] += 360 # 360 degrees

        return angles        

    def rot_angles(self,theta,cart=False):
        """
        NOTE Changes the value of self.j in Vector2d object its called with.

        This function rotates the angular component of a Vector2d in polar 
        coordinates (self.i=radial,self.j=angular) into a new coordinate 
        system where NOTE the positive y-axis is rotated clockwise theta degrees.
        INPUTS:
            -self with self.j in degrees or the vertical Cartesian component.
        REASSIGNS:
            -self.j = self.j-theta"""  
        if cart==False:
            self.j = self.j - theta # degrees
            self.j = self.adjust()
        elif cart == True:
            r_0 = np.sqrt(self.i**2 + self.j**2) # calculates radial comp
            theta_0 = np.arctan2(self.i,self.j) # calculates angle from the positive y
            theta_f = theta_0 - np.deg2rad(theta)

            # find new cartesian components
            self.i = r_0 * np.sin(theta_f)
            self.j = r_0 * np.cos(theta_f) 
    
    def invert(self,cart=False):
        """
        NOTE Changes the value of self.j in Vector2d object its called with.

        Subtracts 180 degrees from the current value of self.j to invert the 
        direction. Used to rotate wind and wave data reported in 
        direction of origin to direction of motion.
        INPUTS:
            -self
        REASSIGNS:
            -self.j"""                                                           
  
        self.rot_angles(180,cart)

    def mag(self):
        """
        This function returns the magnitude of 
        (self.i,self.j) where components are in 
        Cartesian coordiantes.
        INPUTS:
            -self
        RETURNS:
            -magnitude of self"""

        return np.sqrt(self.i**2+self.j**2)

    def polar2cart(self):
        """
        This function converts polar corrdinates 
        with radial and angular components i=r,j=theta
        into Cartesian coordiates i=x,j=y.
        INPUTS:
            -self
        RETURNS:
            -Vector2d object of self in Cartesian coordinates."""

        rad = np.deg2rad(self.j)
        i  = np.sin(rad)*self.i
        j  = np.cos(rad)*self.i
        
        return Vector2d(i,j)

    def cart2polar(self): 
        """
        This function converts Cartesian coordiantes
        with components i=horizontal,j=vertical into 
        polar components i=radial,j=angular.
        INPUTS:
            -self
        RETURNS:
            -Vector2d object of self in polar coordinates."""
        
        r = self.mag()
        self.j = np.arctan2(self.i,self.j)
        theta = self.adjust()

        return Vector2d(r,theta)

    def polar_rot(self,theta):
        """
        Maps polar coordinates speed and direction to a new
        coordinate system with the previous positive y-axis 
        NOTE rotated theta degrees in the clockwise direction.
        INPUTS:
            -self with components in POLAR COORDINATES as i = r 
             and j = theta.
            -theta = Angle of rotation in degrees taken clockwise from the current.
        RETURNS:
            -Vector2d object in Cartesian coordinates with the rotated components."""

        # maps to the new axis
        
        new_traj =self.j - theta
        new_traj = np.deg2rad(new_traj)
            
        # calculates the i and j components
        i = np.sin(new_traj)*self.i
        j = np.cos(new_traj)*self.i
                     
        return Vector2d(i,j)
    

