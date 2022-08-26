import numpy as np


class Climate:
    
    def __init__(self):
        
        self.sst          = None
        self.bottom_temp  = None
        self.atm_pressure = None
        self.air_temp     = None

    def cat_data(self):
        """
        Assumes that all data are 1d and the same length"""
        data = [self.sst,self.bottom_temp,self.atm_pressure,self.air_temp]
        
        has_data = list(map(lambda arr : False if arr is None else True,data))
     
        out = np.zeros((len(self.sst),len(has_data)))
        
        for i in range(len(has_data)):
            if has_data[i]==False:
                pass
            else:
                out[:,i] = data[i]
                
        return out
                
