import numpy as np
from scipy.io import loadmat
from ast import literal_eval
from WindClass import Wind 
from WaveClass import Waves
from CurrentClass import Currents

def readBuoy(id):    

    """
    This function reads information from a buoy_id.txt file
    and returns the attrubutes.
    Called In:
        -Simple.initialize(id != None)."""


    with open(id+'.txt','r') as buoyinfo:
        strinfo = buoyinfo.read()
        info = literal_eval(strinfo)


    winddata = loadmat(info[0])
    wind_ = winddata['wind']
    # times = winddata['time']

    wavedata = loadmat(info[1])
    waves_ = wavedata['waves']

    flowdata = loadmat(info[2])
    csflow_ = flowdata['CrossShelf']
    asflow_ = flowdata['AlongShelf']

    depth = info[3]
    lat = info[4]
    rotation = info[5]

    wind = Wind(wind_[:,0],wind_[:,1])
    waves = Waves(None,waves_[:,1],swh=waves_[:,0],T=waves_[:,2])
    currents = Currents(csflow_,asflow_)

    return wind,waves,currents,depth,lat,rotation

def newBuoy():
    """
    This script saves a dictionary containing data and parameters to make 
    a Buoy class object."""

    print("""This code inputs and stores Buoy parameters into a dictionary and saves them as name.txt.
    To use this information, intialize a buoy object with Buoy('name')""")

    print('Input Buoy name')
    name = input()
    print("Input wind data path or file name if its saved to the current directory.")
    wndpath = input()

    print("Input wave data path.")
    wvpath = input()

    print("Input current data path.")
    curpath = input()
    print('Input depth')
    depth = float(input())

    print("Input Lattitude in degrees.")
    lat = float(input())

    print("Would you like to define a coordinate system? (y/n)")
    check = input()

    yes = ['y','yes','yeah','ye','yea']
    if check in yes:
        print('Input the angle between the new positive y-axis and the existing positiv y-axis.')
        rotation = float(input())
    else:
        rotation = 0
    data = {0:wndpath,1:wvpath,2:curpath,3:depth,4:lat,5:rotation}

    with open(str(name) + '.txt','w') as file:
        file.write(str(data))
        file.close()

    print('Buoy Saved as ' + str(name) + '.txt')
    print('Call Buoy object as Buoy("'+str(name)+'") to use this info.')