import numpy as np
from VectorClass import Vector2d

class Currents(Vector2d):
    """
    DESCRIPTION:
    This class was built to hold 2D current velocity data.
    Vector2d methods assume angles are the closkwise displacement 
    from the positive y-axis.
    NOTE
    CURRENTS ASSUMES self.i AND self.j ARE OF THE FORM:
    -
    -[[velocity_(t0,depth0), ... velocity_(tN,depth0)]
    - :
    - [velocity_(t0,depthN), ... velocity_(tN,depthN)]]
    -
    where each index holds a velocity measurement at the 
    cooresponding depth and time indicies. Time can be 
    replaced with horizontal space and the functions will 
    work the same with depthN in the first row.
    NOTE
    ATTRIBUTES:
        -self.i = cartesian u component of flow.
        -self.j = cartesian v component.
        -self.depth = depth(s) cooresponding to the columns of self.i,self.j 
        -             as float or 1D array. 
    METHODS:
        -Re-Assigns attribute values: self.new_coordsys(y_displacement)
        -
        -Returns Vector2d objects: self.depth_average(),self.transports()
    IMPROVEMENTS:
    Add a third flow dimension w. Allow current values to vary across two
    dimensions of space.
    ============================================================================"""
    def __init__(self,i,j,depth=None):
        super().__init__(i,j)
        self.depth = depth
    
    def new_coordsys(self,y_displacement,cart=True):
        """
        This function applies a clockwise rotation 
        y_displacement degrees from the current 
        positive y axis using Vector2d.rot_angles(theta,cart=True). 
        For typical ADCP data with j=North and i=East components, 
        this function re-definesthe positive y-axis to be 
        y_displacement degrees clockwise from North.
        INPUTS:
            -self in Cartesian coordinates if cart=True (defult)
        REASSIGNS:
            -self.i,self.j = u,v Cartesian components in new
            -                coordinate system if cart=True."""
        
        self.rot_angles(y_displacement,cart=cart)
        print(y_displacement)

    def depth_average(self,real=True):
        """
        This function calculates the depth-averaged
        flow by averaging each column for the i and j 
        components of flow.
        INPUTS:
            -self in Cartesian coordinates
        RETURNS:
            -Vector2d object containing the depth-averaged 
            -currents if real =True (defult). If real = False
            an array of complex vectors (i_da+sqrt(-1)*j_da) is 
            returned."""
        
        i_da = np.nanmean(self.i,axis=0)
        j_da = np.nanmean(self.j,axis=0)

        if real == False:
            out = np.zeros(self.i.shape[1],dtype=np.complex_)
            out.real = i_da
            out.imag = j_da
            return out
        return Vector2d(i_da,j_da) 
        
    def transports(self):
        """
        Calculates the i and j components of transport
        by summing each column of self.i and self.j.
        INPUTS:
            -self in Cartesian coordinates.
        RETURNS:
            -Vector2d object of calculated transports.""" 
        
        i_transport = np.nansum(self.i,axis=0)
        j_transport = np.nansum(self.j,axis=0)

        return Vector2d(i_transport,j_transport)