from .datetimearr import datetime_array

class TimeSeries:
    """
    This class contains functions used to format time series
    data."""
    
    def __init__(self,data):
        self.timestamps = data[:,0]  
        self.data  = data[:,1:] 
        
    
    def sort_1Dtime_series(self,sep='-'):
        """
        This function inputs a time series array"""
                
        dims = self.data.shape
        
        out = datetime_array(years,cols=dims[1],sep=sep) 
        
        timemap = out[:,0]
        
        idx_map = np.linspace(0,len(timemap)-1,len(timemap),dtype=int)
        
        getidx = lambda timestamp : idx_map[timemap==timestamp]

        time_idxs = np.fromiter(map(getidx,self.timestamps),dtype=int)
#         burst = np.zeros((64,1))
        
#         for idx in range(dims[1]):
#             timeidx = getidx(idx_map,datatimes[idx]) 
#             burst[:,0] = data[:,idx]
#             out[:,timeidx] = burst
        out[time_idxs,1:] = self.data
        return out