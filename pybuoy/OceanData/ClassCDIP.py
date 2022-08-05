import numpy as np
import netCDF4 as nc
import time
from .datetimearr import datetime_array


class CDIP:
    
    def __init__(self,station):
        self.link = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/{}p1/{}p1_historic.nc'.format(station,station)
        self.data = nc.Dataset(self.link)
        self.time = self.epoch2stamp(np.array(self.data['waveTime']));
        
    def epoch2stamp(self,epoc,delim='-'):
        """
        docstring"""
        convert = lambda strtm, sp: '20'+strtm[6:8] + sp + strtm[0:2] + sp + strtm[3:5] + sp + strtm[-2:] 

        out = np.zeros_like(epoc,dtype=object)
        loc_times = map(time.gmtime,epoc)
        idx = 0
        for tm in loc_times:
            out[idx] = convert(time.strftime('%D,%H',tm),delim)
            idx += 1

        return out
    
    
    def get_spectrum(self):
        """
        docstring"""
        spec = np.array(self.data['waveEnergyDensity'])
        fbins = np.array(self.data['waveFrequency'])
        spec[spec==-999] = np.nan
        return spec,fbins
    
    def get_bulkstats(self,fields):
        """
        docstring"""
        data = np.zeros((len(self.time),len(fields)))
        
        idx = 0
        for field in fields:
            
            try:
                data[:,idx] = np.array(self.data[field])
            except IndexError:
                print('Invalid value at fields {}. Valid fields displayed below \n'.format(idx))
                for var in self.data.variables.values():
                    print(var)
                return None
            idx += 1
        
        data[data==-999] = np.nan
        
        return [data[:,i] for i in range(data.shape[1])]
        
    def make_timestamp(self,datehr,sep='-'):
        out = ''
        for i in range(len(datehr)):
            out += str(int(datehr[i])) if int(datehr[i]) >=10 else '0' + str(int(datehr[i]))
            out += sep
        return out[:-1]
        
    def make_hrly(self,data,cols,years):
#         years = np.linspace(data[0,0],data[-1,0],int(abs(data[0,0]-data[-1,0]))+1,dtype=int)
        out = datetime_array(years,len(cols)+1)
        idx_map = np.array(range(len(out)))
        
        data[data==99]= np.nan
        data[data==999]=np.nan
        
        idx = 0
        while idx < len(data)-1:
#             datehr = data[idx,:4]
            timestamp = data[idx,0]
            out_idx = idx_map[out[:,0]==timestamp][0]

            def check(step):
                try:
                    return np.all(timestamp == data[idx+step,0])
                except IndexError:
                    return False 
            step = 1
            same_hr = check(step)
            while same_hr == True:
                step += 1 
                same_hr = check(step)
            out[out_idx,1:] =  np.nanmean(data[idx:idx+step,cols],axis=0)
            idx += step        
        return out
                  
    def make_hrly2d(self,data,years,times):
#         years = np.linspace(data[0,0],data[-1,0],int(abs(data[0,0]-data[-1,0]))+1,dtype=int)
        time_map = datetime_array(years,1)
        idx_map = np.array(range(len(time_map)))
        
        data_shape = list(data.shape)
        out_shape = [0 for i in range(2)]
        out_shape[data_shape.index(max(data_shape))] = len(time_map)
        out_shape[data_shape.index(min(data_shape))] = min(data_shape)
        
        out = np.full(out_shape,np.nan)
        
        data[data==99]= np.nan
        data[data==999]=np.nan
        
        idx = 0
        while idx < max(data_shape)-1:
#             datehr = data[idx,:4]
            timestamp = times[idx]
            out_idx = idx_map[time_map[:,0]==timestamp][0]

            def check(step):
                try:
                    return np.all(timestamp == time_map[idx+step])
                except IndexError:
                    return False 
            step = 1
            same_hr = check(step)
            while same_hr == True:
                step += 1 
                same_hr = check(step)
            if data_shape.index(max(data_shape)) == 0:
                out[out_idx,:] =  np.nanmean(data[idx:idx+step,:],axis=0)
            else:
                out[:,out_idx] =  np.nanmean(data[:,idx:idx+step],axis=1)

            idx += step        
        return out
        
        
        
                