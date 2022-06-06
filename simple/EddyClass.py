import numpy as np
from .Coriolis import coriolis 
"""
===========================================================
-                     REFERENCES
-[1] S. J. Lentz,"Sensitivity of the Inner-Shelf Circulation 
-                 to the Form of Eddy Viscosity Profile", 
-                 Journal of Physical Oceanography, 
-                 Vol. 25 (1993) pp. 19-28
==========================================================="""

class EddyVisc:
    """This class calculates various forms of Eddy Viscosity 
    for use in the finite difference scheme in class: Model.
    Contains:
    -ATTRIBUTES:
        -self.Av = Nz x len(data) array of verticle viscosity                       
                    values with each column containing a profile.
        -self.dz = vertical spatial steps.
        -self.profiles = Dictionary of viscosity function used in 
                         creating a viscosity profile in class: model.
    -METHODS:
        -internally called: self.get_depth()
        -externally called: self.BilinearCutoff(),self.Constant()
    Last Edited: 5/11/22
        Changes: wrote docstring
    -TH"""

    def __init__(self):
        
        self.dz = []
        self.Av = []

    def get_depth(self,EkDepth,linear_ext=.1):
        """
        This function calculates the index of the
        vertical gridpoint where the discritization 
        is equivalent to a percentage: linear_ext of 
        the turbulant boundary layer: EkDepth. NOTE Used 
        to make viscosities that increase linearly from 
        a boundary.
        INPUTS:
            -EkDepth = Depth of the boundary layer
            -linear_ext= Percent of the boundary layer
                         wher the viscosity profile is 
                         linearly increasing.
        OUTPUT:
            -z_idx= Index of the desired depth given by 
                    linear_ext*Ekdepth. It will return 
                    the index of the grid point less than 
                    linear_ext*Ekdepth when it lies between 
                    two.
        Called in:
            -self.BilinearCutoff()"""

        ii = True
        z_idx = 0 # "depth" index
        while ii: # while ii == True:
                  # ii = True if the sum of dz[0] through  
                  # dz[z_idx] is less than (linear_ext*EkDepth).
            ii = abs(np.sum(self.dz[:z_idx],axis=0)) < abs(linear_ext*EkDepth)
            z_idx += 1 # iterate
        return z_idx

    def bilinear_cutoff(self,arg_tuple,d=.1,Sr=.01,Br=.01,swdens=1023):
        """This function creates a vertical bilinar cutoff eddy 
        viscosity profile as discussed [1]. NOTE The profile is linear with
        respect to topstress and botstress and extends this way through
        d percent of the turbulant boundary scale delta=(.4*u^*)/f where 
        .4 is Von Karms constant, u^*=sqrt(boundary stress/density) is the 
        shear velocity, and f is the coriolis parameter.
        INPUTS:
            -topstress = Surface Stress. Usually the windstress.
            -botstress = Sea-Floor stress. Currently we just assume 
                         (topstress = botstress)
            -(d=.1) = Percent of the Ekman Scale where the viscosity is changing
                      linearly. Defulted to .1 arbitrarily.
            -(Sr=.01) = "Surface Roughness discussed in [1]. Small
                         distance underneath surface boundary where 
                         the top gridpoint is evaluated.
            -(Br=.01) = "Bottom Roughness discussed in [1]. Same as 
                         surface roughness just for the sea-floor.
        REASSIGNS:
            -self. Av = A 1xNz array containing the bilinear cutoff form of Av
        Last Edited: 5/11/22
        -TH"""
        
        # Lambda: EkmanScale() calculates the turbulant boundary layer scale 
        topstress,botstress,lat = arg_tuple
        
        EkmanScale = lambda stress : (0.4*np.sqrt(stress/swdens))/coriolis(lat)
    
        Av = np.zeros((len(self.dz),len(topstress))) # initialize Av

        Dtop = EkmanScale(topstress) # Surface Ekman Scale
        Dbot = EkmanScale(botstress) # Bottom Ekman Scale
            
        i=0
        while i < len(Dtop):
            Dtop_idx = self.get_depth(Dtop[i]) # find max linear index
            Dbot_idx = self.get_depth(Dbot[i]) # same as above

            # z = array containing depth values, in the units of the
            # water column, for  each gridpoint where Av is linearly 
            # increasing/decreasing.
            z = np.linspace(Sr,.1*Dtop[i],Dtop_idx)

            # Assign linear values of Av in the surface layer.
            Av[:Dtop_idx,i] = (.4 * np.sqrt(topstress[i]/swdens) 
                    *np.array(self.dz[:Dtop_idx]/self.dz[Dtop_idx-1])*d*Dtop[i])
            # Assign linear values of Av in the bottom layer.

            Av[-Dbot_idx:,i] = (.4 * np.sqrt(botstress[i]/swdens)
                    *np.array(self.dz[:Dtop_idx]/self.dz[Dtop_idx-1])[::-1]*d*Dtop[i])

            # Assign intearior values (currently constant)
            Av[Dtop_idx:-Dbot_idx,i] = np.linspace(Av[Dtop_idx-1,i],
                    Av[-Dbot_idx,i],len(Av[Dtop_idx:-Dbot_idx,i]))
            
            i += 1
        return Av

    def constant(self,args):
        """            
        This function inputs a constant value of eddy
        viscosity in m**2/s and outputs a 1 x Nz array
        of that value."""
        value,lendata = args
        return np.full((len(self.dz),lendata),value)

    def eddylist(idx):
        """
        This function inputs an integer 
        idx corrosponding the the form of 
        eddy viscosity in viscdict ={} 
        defined in the makevisc() method 
        in class: Simple (ModelClass.py).
        OUTPUTS:
            -Av = form of eddy viscosity 
                  corospoinding to idx in
                  viscdict."""


        profiles = {'bilincut':0,'Constant:':1}
        visclist = [self.BilinearCutoff,
                    self.Constant]
        return visclist[idx]