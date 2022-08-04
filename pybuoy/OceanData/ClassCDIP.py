import numpy as np
import netCDF4 as nc
import time


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
        
                    
        
        
        
                