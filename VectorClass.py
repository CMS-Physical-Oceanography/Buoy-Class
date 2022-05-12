import numpy as np

class Vector2d:
    """
    Vector2d class functions are build assuming the i component is 
    the horizontal component in cartesian coords and the radial 
    component in polar coords.  Likewise j represents vertical and 
    angular components. Components can be changed with the 
    flipflop(self) class function. Angles are assumed to be taken 
    clockwise from the positive y-axis for dealing with ocean and 
    meteorlogical data
================================================================================"""       
    
    def __init__(self,i,j):
        self.i = i
        self.j = j
        
    """=============================================================================
                                                                           """
    
    def polar2cart(self,zero_axis ='y',return_vector=False):

        rad = np.deg2rad(self.j)
        
        if return_vector == False:

            if zero_axis == 'y':

                i = np.array(list(map(lambda x,y: x*y,np.sin(rad),self.i)))
                self.j = np.array(list(map(lambda x,y: x*y,np.cos(rad),self.i)))
                self.i = i

            elif zero_axis == 'x':
                
                i = np.array(list(map(lambda x,y: x*y,np.cos(rad),self.i)))
                self.j = np.array(list(map(lambda x,y: x*y,np.sin(rad),self.i)))
                self.i = i

            else:

                print('invalid zero_axis input')
                
        elif return_vector == True:
         

            
            i  = np.array(list(map(lambda x,y: x*y,np.cos(rad),self.i)))
            j  = np.array(list(map(lambda x,y: x*y,np.sin(rad),self.i)))
            
            if zero_axis == 'y':
                            
                return(Vector2d(j,i))
                
            elif zero_axis == 'x':
                
                return(Vector2d(i,j))

            else:

                print('invalid input to zero_axis')
                
        else:

            print('invalid input for return_vector')
            

            
        """
================================================================================
                                                                              """

    def timeslice(self,num = 1,period = 'year',segment=True):
        lengths = {'year':(24*365),'month':(24*30),'day':23,'all':len(self.i)}
        slice_ = (lengths[period]*num)-1
        
        if segment == True:
            self.i = self.i[slice_-lengths[period]:slic_]
            self.j = self.j[slice_-lengths[period]:slice_]
        else:
            self.i = self.i[slice_-lengths[period]:-1]
            self.j = self.j[slice_-lengths[period]:-1]

        
        """

=================================================================================
                                                                              """
            
    def cart2polar(self):
        # arctan2() takes the quadrant into account 
        # https://numpy.org/doc/stable/reference/generated/numpy.arctan2.html
        def adjust(angle):
            # Adjusts for negative values

            if angle < 0 :
                return angle + np.pi*2
            else:
                return angle 
        r = np.array(list(map(lambda x,y : np.sqrt(x**2+y**2),self.i,self.j)))
        theta = np.array(list(map(adjust,np.arctan2(self.j,self.i))))

        return Vector2d(r,theta)

    
    """
=================================================================================
                                                                                """      
    def mag(self):
        return np.sqrt(self.i**2+self.j**2)
    
    """
=================================================================================
                                                                               """

    def invert(self,units='deg'):
        """
    Maps a set of angles 180 degrees or pi radians from the current plane.
    Useful for directions that describe where something originates such as 
    wind direction.
                                                                       """
        def adjust(angle):
            # Adjusts for negative values

            if angle < 0 :
                return angle + 360
            else:
                return angle                                                                      
        
        if units == 'deg':
            a  = np.array(list(map(lambda x: x-180,self.j)))
            self.j = np.array(list(map(adjust,a)))
        elif units == 'rad':
            a = np.array(list(map(lambda x: x-np.pi,self.j)))
            self.j = np.array(list(map(adjust,a)))
        else:
            print('INVALID UNITS: keyword units can be either','\n'
           'deg or rad for degrees or radian')
   
        """
================================================================================
                                                                      """       
    
    def adjust(angle):
        # Adjusts for negative values
        if angle < 0 :
            return angle + 360
        else:
            return angle

                              
        """
================================================================================
                                                                              """              

    def rot_angles(self,theta):
        """
    Maps an angle of value {0,360} in degrees from a cardinal
    plane with N = 0 deg theta degrees clockwise of the
    North axis.
    
    Theta is the clockwise angle between the North axis
    (0 deg) and the new positive y axis in degrees
    
                                                    """    
        def adjust(angle):
            # Adjusts for negative values

            if angle < 0 :
                return angle + 360
            else:
                return angle

        #maps to the new plane 
        new_dir = np.array(list(map(lambda x: x-theta,self.j)))
        new_dir = np.array(list(map(adjust,new_dir)))

        self.j = new_dir

                              
        """
================================================================================

                                                                              """       

    def polar_rot(self,theta,deg = True):
        """
    Maps polar coordinates speed and direction {0,360}
    in degrees from a cardinal plane with North = 0 deg to a
    coordinate system with the positive y axis theta degrees
    clockwise from the North axis.
    
    Returns cartesian u and v components of the radial speed.  
    
    Theta is the clockwise angle between the North axis
    (0 deg) and the new positive y axis in degrees
    
                                                   """
       
        
        def adjust(angle):
            # Adjusts for negative values

            if angle < 0 :
                return angle + 360
            else:
                return angle

            
        # maps to the new plane 
        new_traj -= theta
        new_traj = np.array(list(map(adjust,new_traj)))


        # checks if units are in degrees
        if deg == True: 
            new_traj = np.deg2rad(new_traj)
            
        # calculates the u and v components
        u = np.sin(new_traj)*self.i
        v = np.cos(new_traj)*self.i
                     
        return Vector2d(u,v)