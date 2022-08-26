import numpy as np
from .WindClass    import Wind 
from .WaveClass    import Waves
from .CurrentClass import Currents
from .ClimateClass import Climate
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
            self.climate  = Climate()
            self.depth = None
            self.lat = None
            self.timestamps = None
        else:
            from .BuoyHelp import readBuoy
            self.wind,self.waves,self.currents,self.depth,self.lat,rotation = readBuoy(id)
            self.waves.depth,self.currents.depth = [self.depth]*2
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
            wind = np.zeros((2,len(self.wind.i)))
            wind[range(2),:] = self.wind.i,self.wind.j
            np.save(path+buoy_id+'wind.npy',wind)
            
        cli_data = self.climate.cat_data()
        if np.all(cli_data == 0):
            print('no climate data')
        else:
            np.save(path+buoy_id+'climate.npy',cli_data) 

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
        if type(self.timestamps)==None:
            print('no timestamps')
        else:
            np.save(path+buoy_id+'times.npy',self.timestamps)
            
    def readbuoy(self,buoy_id,path=''):
        """
        docstring"""
        def read_meta(buoy,data):
            buoy.timestamps = data
        def read_wind(buoy,data):
            buoy.wind.i = data[0,:]
            buoy.wind.j = data[1,:]
        def read_bulk_waves(buoy,data):
            data[data<0] = np.nan
            buoy.waves.swh = data[0,:]
            buoy.waves.T   = data[1,:]
            buoy.waves.j   = data[2,:]
            
        def read_spec(buoy,spec,fbins):
            buoy.waves.spec = spec
            buoy.waves.fbins= fbins
            
        def read_currents(buoy,data):
            buoy.currents.i = data[:,:,0]
            buoy.currents.j = data[:,:,1]
            
        def read_climate(buoy,data):          
            buoy.climate.sst         = data[:,0]
            buoy.climate.bottom_temp = data[:,1]
            buoy.climate.atm_pressure= data[:,2]
            buoy.climate.air_temp    = data[:,3]
            

        try:
            read_climate(self,np.load(path+buoy_id+'climate.npy',allow_pickle=True))
        except FileNotFoundError:
            print('no climate data')    
        try:
            read_wind(self,np.load(path + buoy_id + 'wind.npy',allow_pickle= True))
        except FileNotFoundError:
            print('no wind data')
        try:
            read_bulk_waves(self,np.load(path+buoy_id+'waves.npy',allow_pickle=True))
        except FileNotFoundError:
            print('no bulk waves data')
        try:
            read_spec(self,np.load(path+buoy_id+'wvspec.npy',allow_pickle=True),np.load(path + buoy_id+'fbins.npy',allow_pickle=True))
        except FileNotFoundError:
            print('no wave spectrum data')
        try:
            read_currents(self,np.load(path + buoy_id+'currents.npy',allow_pickle=True))
        except FileNotFoundError:
            print('no currents data')
        try:
            read_meta(self,np.load(path+buoy_id+'times.npy',allow_pickle=True))
        except FileNotFoundError:
            print('no time stamps') 
            
    def filter_currents(self,cut_freq,N):
        """
        docstring"""
        from pandas import DataFrame
        from scipy.signal import butter,filtfilt

        A,B = butter(N,cut_freq)
        out_i,out_j = [np.zeros_like(self.currents.i) for i in range(2)]
        nan_idxs = np.zeros_like(self.currents.i,dtype=bool)
        dims = self.currents.i.shape
        idx_map = np.array(range(max(dims)),dtype=int)

        for i in range(min(dims)):
            if np.all(np.isnan(self.currents.i[i,:])):
                pass
            else:
                start_idx = 0
                check = np.isnan(self.currents.i[i,:][start_idx])
                while check:
                    start_idx += 1
                    check = np.isnan(self.currents.i[i,:][start_idx])

                nan_idxs[i,:] = np.isnan(self.currents.i[i,:])
                interp_i = np.array(DataFrame(self.currents.i[i,start_idx:]).interpolate())[:,0]
                interp_j = np.array(DataFrame(self.currents.j[i,start_idx:]).interpolate())[:,0]
                filt_i = filtfilt(A,B,interp_i)
                filt_j = filtfilt(A,B,interp_j)

                self.currents.i[i,start_idx:] = filt_i
                self.currents.j[i,start_idx:] = filt_j

        self.currents.i[nan_idxs] = np.nan
        self.currents.j[nan_idxs] = np.nan
        
    def timeslice(self,start,end):

        idx_map = np.array(range(len(self.timestamps)))
        start_idx = idx_map[self.timestamps==start][0]
        end_idx = idx_map[self.timestamps==end][0] +1

        try:
            self.wind.i = self.wind.i[start_idx:end_idx]
            self.wind.j = self.wind.i[start_idx:end_idx]
        except IndexError:
            pass
        try:
            self.waves.j = self.waves.j[start_idx:end_idx]
            self.waves.swh = self.waves.swh[start_idx:end_idx]
            self.waves.T = self.waves.T[start_idx:end_idx]

            t_ax = list(self.waves.spec.shape).index(max(self.waves.spec.shape))
            if t_ax==0:
                self.waves.spec = self.waves.spec[start_idx:end_idx,:]
            else:
                self.waves.spec = self.waves.spec[:,start_idx:end_idx]
        except IndexError:
            pass
        self.timestamps = self.timestamps[start_idx:end_idx]

