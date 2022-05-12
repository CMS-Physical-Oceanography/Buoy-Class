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

    def initialize(self,h):

        wndpath = 'c:/users/twhes/onedrive/documents/ob27shelfflow/ndbc41037/WindDatalpf.txt'

        winddata = np.loadtxt(wndpath,dtype=object,delimiter=',')
        wind = np.array(winddata[:,1:],dtype=float)[list(winddata[:,0]).index('2013.1.1.0'):,:]

        flowpath = 'c:/users/twhes/onedrive/documents/ob27shelfflow/ob27currents/FilteredCurrentData_2013-2021.mat'
        
        from scipy.io import loadmat
        flowdata = loadmat(flowpath)
        csflow = flowdata['CrossShelf']
        asflow = flowdata['AlongShelf']

        self.wind = Wind(wind[:len(csflow[0]),0],wind[:len(csflow[0]),1])
        self.currents = Currents(csflow,asflow)
        self.depth = h
    
    def Model(self):
        """
        This funciton intializes a model class using wind stress
        and depth-averaged flow data from Buoy class."""
        flowda = self.currents.depth_average()
        wndstress = self.wind.windstress()
        return Simple(self.depth,wndstress,flowda)
        
