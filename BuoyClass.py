import numpy as np
from WindClass import Wind
from CurrentClass import Currents
from WaveClass import Waves
from BinClass import Bins

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
            from BuoyHelp import readBuoy
            self.wind,self.waves,self.currents,self.depth,self.lat,rotation = readBuoy(id)
            self.waves.depth,self.currents.depth = [self.depth]*2
            self.wind.invert()
            self.waves.invert()
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

    def Model(self):
        """
        This funciton intializes a model class using wind stress
        and depth-averaged flow data from Buoy class."""
        from ModelClass import Simple

        flowda = self.currents.depth_average()
        wndstress = self.wind.windstress(cart=False)
        return Simple(self.depth,self.lat,wndstress,flowda)
    
    def makebuoy(self):
        from BuoyHelp import newBuoy
        newBuoy()

