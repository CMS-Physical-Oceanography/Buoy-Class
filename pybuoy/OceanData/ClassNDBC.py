import requests
from bs4 import BeautifulSoup
from .ScrapeClass import DataScrape
from .TimeSeriesClass import TimeSeries
from .datetimearr import datetime_array
import numpy as np
import time

class NDBC(TimeSeries):
    
    def __init__(self,station):
        
        self.baselink   = 'https://www.ndbc.noaa.gov//'
        self.exts = {'history': 'station_history.php?station={}'.format(station)}
    
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
            datehr = data[idx,:4]
            timestamp = self.make_timestamp(datehr)
            out_idx = idx_map[out[:,0]==timestamp][0]

            def check(step):
                try:
                    return np.all(datehr == data[idx+step,:4])
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
    
    
            
    def get_histdata(self,datatypes,years,printlink=False):
        """
        This function returns a numpy array containing the raw? data
        call DataScrap here
        
        """
        lap0 = time.time()
        fields = {}
        for data in datatypes:
            fields[data] = list(map(str,years))

        scraper = DataScrape(self.baselink,fields=fields)
        
        goodhrefs = scraper.layer1_and_search(ext=self.exts['history'])
#         print('initial connection made',time.time()-lap0)

        output = np.zeros(len(goodhrefs),dtype=object)

        histfilt = lambda goodhrefs : filter(lambda href : 'histor' in href,
                                             goodhrefs)
        extensions = map(histfilt,
                         goodhrefs)
        
        field_idx = 0
        for field in extensions:
                i = 0
                has_data = False # avoids copies in output when theres no availible data.
                
                
                for ext in field:
                    has_data = True
                    
                    if int(years[i]) < 2007:
                        header = True
                    else:
                        header = False
                        
                    if i == 0:                        
                        href = scraper.layer1_search(fields=['view_text'],ext=ext)[0][0]
                        response = scraper.make_request(self.baselink,href)
                        
                        if header == False:
                            data = scraper.read_txtarry(response)
                        else:
                            data = scraper.read_txtarry(response)
                            data = np.array(data[1:,:],dtype=float)
                            
                        if printlink == True:    
                            print(self.baselink+href,time.time()-lap0)             
                    else:                        
                        href = scraper.layer1_search(fields=['view_text'],ext=ext)[0][0]
                        response = scraper.make_request(self.baselink,href) 
                        try:
                            if header == False:
                                
                                data = np.concatenate((data,scraper.read_txtarry(response)))
                            else:
                              
                                dat = scraper.read_txtarry(response)
                                dat = np.array(dat[1:,:],dtype=float)
                                data = np.concatenate((data,dat))
                                
                        except ValueError:
                            print('Inconsistent Array Size in the link below')
                            pass
                        if printlink == True:    
                            print(self.baselink+href,time.time()-lap0)                        
                    i += 1     
                    
                output[field_idx] = data if has_data == True else []
                field_idx += 1
        if np.all(years > output[0][0,0]):
            return output[0][1:]
        else:
            return output[0]
                

            
     
                                  
