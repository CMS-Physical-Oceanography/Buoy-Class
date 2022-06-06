import requests
from bs4 import BeautifulSoup
from .ScrapeClass import DataScrape
from .TimeSeriesClass import TimeSeries
import numpy as np
import time

class NDBC(TimeSeries):
    
    def __init__(self,station):
        
        self.baselink   = 'https://www.ndbc.noaa.gov//'
        self.exts = {'history': 'station_history.php?station={}'.format(station)}
    
            
    def get_histdata(self,datatypes,years):
        """
        This function returns a numpy array containing the raw? data
        call DataScrap here
        
        """
        lap0 = time.time()
        fields = {}
        for data in datatypes:
            fields[data] = years

        scraper = DataScrape(self.baselink,fields=fields)
        
        goodhrefs = scraper.layer1_and_search(ext=self.exts['history'])
        print('initial connection made',time.time()-lap0)

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
                    if i == 0:                        
                        href = scraper.layer1_search(fields=['view_text'],ext=ext)[0][0]
                        response = scraper.make_request(self.baselink,href)
                        data = scraper.read_txtarry(response)
                        print(self.baselink+href,time.time()-lap0)                        
                    else:                        
                        href = scraper.layer1_search(fields=['view_text'],ext=ext)[0][0]
                        response = scraper.make_request(self.baselink,href) 
                        try:
                            data = np.concatenate((data,scraper.read_txtarry(response)))
                        except ValueError:
                            print('Inconsistent Array Size in the link below')
                            pass
                            
                        print(self.baselink+href,time.time()-lap0)                        
                    i += 1     
                    
                output[field_idx] = data if has_data == True else []
                field_idx += 1
                
        return output
                

            
     
                                  
