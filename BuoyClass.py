import numpy as np
from WindClass import Wind
from CurrentClass import Currents
from WaveClass import Waves
from ModelClass import Simple
from Coriolis import coriolis
class Buoy:
    def __init__(self,h):
        #=================
        # These will be
        # assigned in
        # initialize(h)
        self.wind = []
        self.currents = []
        self.depth = 0
        #==================
        self.initialize(h)

    
    def Model(self):
        """
        This funciton intializes a model class using wind stress
        and depth-averaged flow data from Buoy class."""
        flowda = self.currents.depth_average()
        wndstress = self.wind.windstress()
        return Simple(self.depth,wndstress,flowda)
        
