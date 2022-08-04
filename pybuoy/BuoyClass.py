import numpy as np
from .WindClass    import Wind 
from .WaveClass    import Waves
from .CurrentClass import Currents 
from .BinClass     import Bins
from .OceanData    import NDBC,CORMP


class Buoy(Bins):
    """
    Write docstring."""

    def __init__(self,id=None):
        super().__init__()
        self.initialize(id)

    def initialize(self,id):

        self.bins = Bins()
        if id ==None:
            self.wind = Wind(None,None)
            self.waves = Waves(None,None)
            self.currents = Currents(None,None)
            self.depth = None
            self.lat = None
        else:
            from .BuoyHelp import readBuoy
            self.wind,self.waves,self.currents,self.depth,self.lat,rotation = readBuoy(id)
            self.waves.depth,self.currents.depth = [self.depth]*2
 #             self.wind.invert()
#             self.waves.invert()
            self.rotate_buoy(rotation)

    def rotate_buoy(self,y_displacement):
        """
        This function rotates the coordinate system 
        of the buoy y_displacement degrees clockwise 
        from the current y-axis. Wind and Wave angles 
        are adjusted to stay in [0,360], currents are 
        kept in Cartesian coordinates.
        INPUTS:
            -y_displacement = angle of rotation clockwise 
            -                 from the current positive y-axis
            -                 in degrees.
        REASSIGNS:
            -self.wind.j,self.waves.j = (self.wind.j,self.waves.j) - theta
            -self.currents.i,self.currents.j = u,v Cartesian coordinates in new 
            -                                  coordinate system."""

        self.wind.new_coordsys(y_displacement,cart=False)
        self.waves.new_coordsys(y_displacement,cart=False)
        self.currents.new_coordsys(y_displacement,cart=True)
        
    def savebuoy(self,buoy_id,path=''):
        """
        docstring."""
   
        if len(self.wind.i.shape) == 0 or len(self.wind.j.shape) == 0:
            print('no wind data')
        else:
            wind = np.zeros((len(self.wind.i),2))
            wind[range(2),:] = self.wind.i,self.wind.j
            np.save(path+buoy_id+'wind.npy',wind)
            
        if len(self.waves.T.shape) == 0 or len(self.waves.j.shape) == 0 or len(self.waves.swh.shape) == 0:
            print('no wave data')
        else:
            bulk_waves = np.zeros((3,len(self.waves.j)))
            bulk_waves[range(3),:] = self.waves.swh,self.waves.T,self.waves.j
            np.save(path+buoy_id+'waves.npy',bulk_waves)
        
        if len(self.waves.spec.shape) == 0:
            print('no wave spectrum data')
        else:
            np.save(path+buoy_id+'wvspec.npy',self.waves.spec)
            np.save(path+buoy_id+'fbins.npy',self.waves.fbins) 
        
        if len(self.currents.i.shape) == 0 or len(self.currents.j.shape) == 0:
            print('no current data')
        else:
            currents = np.zeros((self.currents.i.shape[0],self.currents.i.shape[1],2))
            currents[:,:,0] = self.currents.i
            currents[:,:,1] = self.currents.j
            np.save(path+buoy_id+'currents.npy',currents)
            
    
    def makebuoy(self):
        from .BuoyHelp import newBuoy
        newBuoy()

