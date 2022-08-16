import numpy as np

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