import numpy as np
import pandas as pd
from scipy.signal import butter,filtfilt

def filter_1d(data,cut_freq,N,btype='lowpass'):
    A,B = butter(N,cut_freq,btype=btype)
    out = np.zeros_like(data)
    nan_idxs = np.isnan(data)
    start_idx = 0
    is_nan = True if np.isnan(data[start_idx])==True else False
    while is_nan:
        start_idx +=1
        is_nan = True if np.isnan(data[start_idx])==True else False
    interp = np.array(pd.DataFrame(data[start_idx:]).interpolate())[:,0]
    filt = filtfilt(A,B,interp,axis=0)
    out[start_idx:] = filt
    out[nan_idxs] = np.nan
    
    return out

def filter_comps(data,A,B):
    comps      = np.zeros_like(data)
    comps[:,0] = np.cos(np.deg2rad(data[:,0]))*data[:,1]
    comps[:,1] = np.sin(np.deg2rad(data[:,0]))*data[:,1]
    
    idx_map = np.array(range(len(data)),dtype=int)
    
    nan_idxs = idx_map[np.isnan(data[:,0])]
    
    df = pd.DataFrame(comps)
    filtcomps = filtfilt(A,B,np.array(df.interpolate()),axis=0)
    
    mag  = np.sqrt(filtcomps[:,0]**2 + filtcomps[:,1]**2)
    traj = np.rad2deg(np.arctan2(filtcomps[:,1],filtcomps[:,0]))

    traj = np.fromiter(map(lambda x : x if x >= 0 else x+360,traj),dtype=float)
    mag[nan_idxs] = np.nan
    traj[nan_idxs] = np.nan
    
    return mag,traj