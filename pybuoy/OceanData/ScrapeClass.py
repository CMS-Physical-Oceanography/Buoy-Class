import requests
from bs4 import BeautifulSoup
import numpy as np

class DataScrape:
    
    def __init__(self,baselink,fields={}):
        self.baselink = baselink
        self.fields = fields
        

    def make_request(self,link,ext=None):
        """
        This function uses requests.get() to 
        get html data from self.baselink\ext
        and returns a requests object.
        
        The way the base link and extension are 
        combined will probably cause problems."""
        
        if ext==None:
            return requests.get(link)
        else:
            return requests.get(link + ext)
    
    
    def get_hrefs(self,response):
        """
        This function returns a list of the 
        href links in a list of html <a> tags.
        from a requests response using 
        
        BeautifulSoup().find_all('a')

        We filter out values with NoneType 
        using the the list comprehension statement: 

        [atag.get('href') for atag in soup.find_all('a')]."""

        soup = BeautifulSoup(response.content,features="html.parser")

        return list(filter(
                           lambda href: href != None,
                           [atag.get('href') for atag in soup.find_all('a')]))
    
    def read_txtarry(self,response):
        """
        This function inputs a requests.get() object
        and returns a numpy array containing the data 
        array."""

        soup = BeautifulSoup(response.content,features='html.parser')
        
        with open('text.txt','w') as fil:
            fil.write(soup.contents[0])
            fil.close()
        return np.loadtxt('text.txt')
    
    
    def set_fields(self,fields=None):
        
        if fields == None:
            try:
                datafields = list(self.fields.keys())
            except AttributeError:
                datafields = self.fields
        else:
            try:
                datafields = list(fields.keys())
            except AttributeError:
                datafields = fields
                
        return datafields
    
    
    def layer1_search(self,fields=None,ext=None):
        """
        This function makes a request to the 
        base link + ext[extension] if 
        extension != None, and returns a list 
        of href strings for keys in the 
        datakeys dictionary."""

        response = self.make_request(self.baselink,ext)
        
        hrefs = self.get_hrefs(response)

        datafields = self.set_fields(fields=fields)  
    
        goodhrefs = [ [] for i in range(len(datafields))]

        for href in hrefs:

            checklist = list(map(lambda field : 1 if field in href else 0,
                                 datafields))

            try:
                fieldidx = checklist.index(1)
                goodhrefs[fieldidx].append(href)                
            except ValueError:
                pass
                
        return goodhrefs
    
    
    def layer1_and_search(self,ext=None):
        """
        This is a one level search that retruns hrefs 
        from a links <a> tags which satisfy 
        self.fields.keys() in href, and 
        self.fields['key'] in href 
        for all keys."""
        
        goodhrefs = self.layer1_search(ext=ext)
        
        output = [ [] for i in range(len(goodhrefs))]
        
        idx = 0
        for key in self.fields.keys():
            
            subfields = self.fields[key]
            hrefs = goodhrefs[idx]
            
            for href in hrefs:
                
                output[idx].append(href) if 1 in map(lambda subfield : 1 if subfield in href else 0,subfields) else 0
                
            idx += 1
            
        return output
        
        
    def layer2_search(self,goodhrefs):
        """
        This function is slow but pretty general.
        
        This is a two level search function which first 
        checks the <a> tags of a link and finds extensions
        containing self.fields.keys(), then queries those 
        extensions and find the results where self.fields[key] 
        is in the href and returns those links."""

        goodsubhrefs = [ [] for i in range(len(goodhrefs))]

        fieldkeys = list(self.fields.keys())

        for key_idx in range(len(fieldkeys)):

            key = fieldkeys[key_idx]
            
            fields_ = self.fields[key]
            
            good_subhrefs = []
    
            for href in goodhrefs[key_idx]:
            
                subhrefs_ = self.layer1_search(fields=fields_,ext=href)
                
                subhrefs =  list(filter(
                                   lambda subhrefs : subhrefs != [],
                                   subhrefs_))
                if subhrefs != []:
                    goodsubhrefs[key_idx].append(subhrefs)

            key_idx +=1

        return goodsubhrefs