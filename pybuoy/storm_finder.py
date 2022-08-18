import numpy as np
import pandas as pd

def find_continuity(data,condition,duration):
    out = []
    idx_map = np.array(range(len(data)),dtype=int)
        
    idxs = idx_map[condition(data)]
    i=0
    while i < len(idxs):
        try:
            if idxs[i+duration] - idxs[i] == duration:
                step = duration
                check = condition(data[idxs[i]+step])

                while check == True:
                    step += 1
                    check = condition(data[idxs[i]+step])

                out.append((idxs[i],idxs[i]+step))

                i += step
            else:
                i +=1

        except IndexError:
            i = len(idxs)

            return out
    return out

def check_continuity(data,condition,duration):
    idx_map = np.array(range(len(data)),dtype=int)
    idxs = idx_map[condition(data)]
    if len(data) < duration:
        return False 
    else:
        j = 0
        while j < len(idxs):
            try:
                if idxs[j+duration] - idxs[j] == duration:
                    return True
                else:
                    j += 1
            except IndexError:
                return False 
        return False 

def storm_search(data,thresh,dur,end_thresh,end_dur):
         

    events = find_continuity(data,lambda x : x > thresh,dur) # find continuous events that meet thresh
    
    i = 0
    storms=[]
    check = []
    while i < len(events)-1:
        sli = data[events[i][1]:events[i+1][0]]
        storm_passed = check_continuity(sli,lambda x : x < end_thresh,end_dur) # check if storm ende
        check.append(storm_passed)
        if storm_passed == True:
            storms.append(events[i])
            i += 1
        else:
            storm_events = []
            while storm_passed == False and i < len(events)-2:
                storm_events.append(events[i])
                i += 1
                sli = data[events[i][1]:events[i+1][0]]
                storm_passed = check_continuity(sli,lambda x : x < end_thresh,end_dur)
                check.append(storm_passed)
            storm_events.append(events[i])
            
            summs = []
            for evnt in storm_events:
                summs.append(np.nansum(data[evnt[0]:evnt[1]]))
                
            storms.append(storm_events[summs.index(max(summs))]) #append the period with the largest integrated stress to 
                                                                # storms
            i += 1
    storms.append(events[-1]) 
    print(len(storms),check.count(True))
    return storms

def integrate_storms(check,data,storms,timestamps,dt=1):
    sums = np.zeros(len(storms))
    i = 0
    for storm in storms:
        if storm == 0:
            pass
        else:
            sums[i]  = np.nansum(data[storm[0]:storm[1]])
            real_count = list(np.isnan(check[storm[0]:storm[1]])).count(False)
            if real_count/len(data[storm[0]:storm[1]]) < .1:
                sums[i] = 0
        i += 1
    storms = np.array(storms,dtype=object)

    sorted_sums = sums.copy()
    sorted_sums.sort()  
    sorted_storms = list(map(lambda ranked_sum : storms[sums==ranked_sum][0],sorted_sums))
    
    return sorted_sums[::-1],sorted_storms[::-1]

def wvstorm_table(buoy,thresh,dur,end_thresh,end_dur,wndbuoy,roughness=.000675,Nstorms=15,storm_names=None):
    """
    docstring"""
    
    header = ['Start Date','Storm Name','Duration (hrs)',r'$\int\tau_{wvs} dt$  (Pa h)',r'Max $\tau_{wvs}$ (Pa)',
              r'Avg. $u_{br}$ (m/s)','Avg. SWH (m)', r'Avg. $T_p$ (s)','Avg. MWD (deg)',r'Avg. $\tau_{wnd}$ (Pa)',
              'Avg. Wind Dir. (deg)']
    
    def avg_angles(angles,storm):
        traj = np.deg2rad(angles[storm[0]:storm[1]])

        x = np.nanmean(np.cos(traj))
        y = np.nanmean(np.sin(traj))
        ang = np.rad2deg(np.arctan2(y,x))
        return ang if ang >=0 else 360 +ang
        
    storm_mean = lambda storm,data : np.nanmean(data[storm[0]:storm[1]])
    fmt_stamp = lambda stamp : stamp[5:7]+ '-' + stamp[8:11] + stamp[:4]
    
    bot_strs,ubr,T_b = buoy.waves.bottom_stress(buoy.waves.depth,roughness)
    int_strs = np.array(pd.DataFrame(bot_strs).interpolate())[:,0]
    wv_storms = storm_search(int_strs,thresh,dur,end_thresh,end_dur)
    Idata,ordered_storms = integrate_storms(buoy.waves.swh,int_strs,wv_storms,buoy.timestamps)
    
    if storm_names is None:
        storm_names = ['']*len(Idata)
    
    peak_strs = list(map(lambda storm : np.nanmax(bot_strs[storm[0]:storm[1]]) if len(storm) > 0 else 0,ordered_storms))
    duration = list(map(lambda storm : storm[1]-storm[0],ordered_storms))
    wnd_strs = wndbuoy.wind.windstress().mag()
    stamps = list(map(fmt_stamp,buoy.timestamps))
    times = list(map(lambda storm : stamps[storm[0]],ordered_storms))
    
    mean_strs,mean_ubr,mean_swh,mean_tp,mean_mwd,mean_wndstrs,mean_wnddir = [ [] for i in range(7)]
    
    for storm in ordered_storms:
        mean_strs.append(storm_mean(storm,bot_strs))
        mean_swh.append(storm_mean(storm,buoy.waves.swh))
        mean_tp.append(storm_mean(storm,buoy.waves.T))
        mean_wndstrs.append(storm_mean(storm,wnd_strs))
        mean_ubr.append(storm_mean(storm,ubr))
        mean_mwd.append(avg_angles(buoy.waves.j,storm))
        mean_wnddir.append(avg_angles(wndbuoy.wind.j,storm))
        
        
    
    Ncols = len(header)
    table_arr = np.zeros((len(ordered_storms),Ncols),dtype=object)
    table_arr[:,0] = times
    table_arr[:len(storm_names),1] = storm_names
    table_arr[:,2] = duration
    table_arr[:,3] = Idata
    table_arr[:,4] = peak_strs
    table_arr[:,5] = mean_ubr
    table_arr[:,6] = mean_swh
    table_arr[:,7] = mean_tp
    table_arr[:,8] = mean_mwd
    table_arr[:,9] = mean_wndstrs
    table_arr[:,10] = mean_wnddir
    
    table = pd.DataFrame(table_arr[:Nstorms])
    table.columns = header
    table.index = range(1,Nstorms+1)
    return table,ordered_storms

def wndstorm_table(buoy,thresh,dur,end_thresh,end_dur,wvbuoy,Nstorms=15,storm_names=None):
    
    header = ['Start Date','Storm Name','Duration (hrs)',r'$\int\tau_{wnd} dt$  (Pa h)',r'Max $\tau_{wnd}$ (Pa)',
              r'Avg. $\tau_{wnd}$ (Pa)','Avg. Wind Dir. (deg)','Avg. SWH (m)',r'Avg. $T_p$ (s)','Avg. $D_p$ (deg)']
    
    def avg_angles(angles,storm):
        traj = np.deg2rad(angles[storm[0]:storm[1]])

        x = np.nanmean(np.cos(traj))
        y = np.nanmean(np.sin(traj))
        ang = np.rad2deg(np.arctan2(y,x))
        return ang if ang >=0 else 360 +ang
        
    storm_mean = lambda storm,data : np.nanmean(data[storm[0]:storm[1]])
    fmt_stamp = lambda stamp : stamp[5:7]+ '-' + stamp[8:11] + stamp[:4]

    wnd_strs = buoy.wind.windstress().mag()
    int_strs = np.array(pd.DataFrame(wnd_strs).interpolate())[:,0]
    wnd_storms = storm_search(int_strs,thresh,dur,end_thresh,end_dur)
    Idata,ordered_storms = integrate_storms(buoy.wind.i,int_strs,wnd_storms,buoy.timestamps)
    
    
    if storm_names is None:
        storm_names=['']*len(Idata)
    
    peak_strs = list(map(lambda storm : np.nanmax(wnd_strs[storm[0]:storm[1]]),ordered_storms))
    duration = list(map(lambda storm : storm[1]-storm[0],ordered_storms))
    stamps = list(map(fmt_stamp,buoy.timestamps))
    times = list(map(lambda storm : stamps[storm[0]],ordered_storms))
    
    mean_strs,mean_ubr,mean_swh,mean_wnddir,mean_wvdir,mean_Tp = [[] for i in range(6)]
    for storm in ordered_storms:
        mean_strs.append(storm_mean(storm,wnd_strs))
        mean_swh.append(storm_mean(storm,wvbuoy.waves.swh))
        mean_wnddir.append(avg_angles(buoy.wind.j,storm))
        mean_wvdir.append(avg_angles(wvbuoy.waves.j,storm))
        mean_Tp.append(storm_mean(storm,wvbuoy.waves.T))
    
    Ncols = len(header)
    table_arr = np.zeros((len(ordered_storms),Ncols),dtype=object)
    table_arr[:,0] = times
    table_arr[:len(storm_names),1] = storm_names
    table_arr[:,2] = duration
    table_arr[:,3] = Idata
    table_arr[:,4] = peak_strs
    table_arr[:,5] = mean_strs
    table_arr[:,6] = mean_wnddir
    table_arr[:,7] = mean_swh
    table_arr[:,8] = mean_Tp
    table_arr[:,9] = mean_wvdir
    
    
    table = pd.DataFrame(table_arr[:Nstorms])
    table.columns = header
    table.index = range(1,Nstorms+1)
    return table