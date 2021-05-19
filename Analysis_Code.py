"""CODE AUTHOR: WILLIAM CROOK"""
#---------------------------------------------------------------------------------------------------------------#
"""IMPORTS"""
#---------------------------------------------------------------------------------------------------------------#

import lightkurve as lk
import numpy as np
from scipy.signal import find_peaks
import copy
import pandas as pd

#---------------------------------------------------------------------------------------------------------------#
"""MAIN PROGRAM"""
#---------------------------------------------------------------------------------------------------------------#

def MainProgram(Star,lc_list,DataType,mode):
    global starttime, noisemask, exclusion, df_list, strongest_period, refined_outliers
    df_list = []
    flattened_list = []
    need_new_options = False
    
    for lightcurve in lc_list:
        lc = copy.deepcopy(lightcurve[2])
        exclusion = [False]*len(lc.time.value)
        
        if DataType == 'FFI':
            BasicClean(lc.flux.value)
            RemoveExtremeOutliers(lc.time.value,lc.flux.value,0.95)
            Exclude(lc)
            xcluded = [lc.time.value[i] if exclusion[i] == True else None for i in range(len(exclusion))]
            ycluded = [lc.flux.value[i] if exclusion[i] == True else None for i in range(len(exclusion))]
            lc.flux.value[exclusion] = 1
            lc_flat = lc
            
        elif DataType == 'SPOC':
            xcluded = [None]
            ycluded = [None]
            noisemask = 1
            new_time, new_flux, _, _ = MovingAverage(lc.time.value, 
                                               lc.flux.value, 
                                               5)
            new_errors = np.array([])
            for ind in range(0,len(lc.flux_err.value)-5,6):
                new_errors = np.append(new_errors,np.mean(lc.flux_err.value[ind:ind+5]))
            lc = lk.LightCurve(time=new_time,flux=new_flux,flux_err=new_errors)
            lc_flat = lc.flatten()
            
        else:
            xcluded = [None]
            ycluded = [None]
            lc_flat = lc.flatten()
        
        BasicClean(lc_flat.flux.value)
        flattened_list.append([lightcurve[0],int(lightcurve[1]),noisemask,lc_flat,xcluded,ycluded])
       
    con_sectors =  FindConsecutive(flattened_list)
    
    grouped_lightcurves = []
    start = 0
    for ind in range(len(con_sectors)):
        grouped_lightcurves.append(flattened_list[start:start+len(con_sectors[ind])])
        start = start + len(con_sectors[ind])
        
    for group in range(len(con_sectors)):   
        strongest_period = 0
        if len(con_sectors[group]) > 1:
            time_range = ([[x[3].time.value[0],x[3].time.value[-1]] for x in grouped_lightcurves[group]])
            noisemask = np.mean([lightcurve[2] for lightcurve in grouped_lightcurves[group]])
            lc_collection = lk.collections.LightCurveCollection([lightcurve[3] for lightcurve in con_sectors[group]])
            lc_flat = lc_collection.stitch()
        else:
            time_range = [[grouped_lightcurves[group][0][3].time.value[0],grouped_lightcurves[group][0][3].time.value[-1]]]
            noisemask = grouped_lightcurves[group][0][2]
            lc_flat = con_sectors[group][0][3]
        
        timespan = lc_flat.time.value[-1]-lc_flat.time.value[0]
        try:
            period_options = PeriodOptions(lc_flat,timespan)
        except:
            period_options = []
        if len(period_options) == 0:
            need_new_options = True
        
        xnew, ynew, peaks, _ = MovingAverage(lc_flat.time.value, 
                                             lc_flat.flux.value,  
                                             2)
        all_peaks = CheckPeaks(xnew,ynew,peaks)
        Transits = GroupTransitPoints(ynew,xnew,all_peaks)
     
        matches_per_batch = []
        for batch in Transits:
            if need_new_options == True and len(batch[1]) > 1:
                refined_outliers = False
                not_outlier, _, _ = RemoveOutliers(np.diff(batch[1]))
                strongest_period = np.mean(not_outlier)
                period_options = [strongest_period]
                
            match_list = BuildTree(batch[0],batch[1],batch[2],period_options,timespan)
            matches_per_batch.append(match_list)

        period_list, peak_mask_list, raw_peaks = CheckValidity(matches_per_batch,all_peaks,xnew,ynew,timespan,time_range)
            
        for ind in range(len(period_list)):
            peak_mask = peak_mask_list[ind]
            peaks = raw_peaks[ind][peak_mask]
            
            lc_forfolding = copy.deepcopy(lc_flat)
            
            time_peaks = [np.where(lc_flat.time.value==xnew[peak])[0][0] for peak in peaks]
            exclusion = [False]*len(lc_flat.time.value)
            for peak in time_peaks:
                StepBack(lc_flat.time.value,lc_flat.flux.value,[peak],1)
                
            isolated_transits = np.where(np.array(exclusion) == True)[0]
            below_noise = [ind for ind in range(len(lc_flat.flux.value)) if lc_flat.flux.value[ind] < noisemask and ind not in isolated_transits]
            lc_forfolding.flux.value[below_noise] = 1
                    
            start_time = xnew[peaks][0]
            
            if need_new_options == False:
                target_period = period_list[ind]
                period_range = np.linspace(target_period*0.95,target_period*1.05,int(500*target_period*0.1))
                pg = lc_forfolding.to_periodogram(method='bls', period=period_range, frequency_factor=500)
                
                period = pg.period_at_max_power.value
            else:
                pg = None
                period = period_list[ind]
            drop_in_flux, drop_in_flux_err, drop_in_flux_graphs = DropInFlux(lc_forfolding,
                                                                             start_time,
                                                                             period,
                                                                             xnew[peaks],
                                                                             ynew[peaks])                   
            
            if mode == 1 and drop_in_flux != False:
                SaveToDataFrame(mode,
                                period,
                                drop_in_flux,
                                drop_in_flux_err,
                                grouped_lightcurves[group],
                                peaks,
                                xnew,
                                ynew,
                                drop_in_flux_graphs,
                                pg)
            elif drop_in_flux != False:
                SaveToDataFrame(mode,
                                period,
                                drop_in_flux,
                                drop_in_flux_err,
                                grouped_lightcurves[group],
                                None,None,None,None,None)
           
    Possible_Exoplanets = GroupExoplanets(df_list)
    results_list = AverageValues(Possible_Exoplanets, Star, mode)
       
    return Possible_Exoplanets, results_list

def SaveToDataFrame(mode,period,drop_in_flux,drop_in_flux_err,lightcurves,peaks,xnew,ynew,drop_in_flux_graphs,pg):
    global df_list
    if len(lightcurves) == 1:
        sector_range = str(lightcurves[0][0])
    else:
        sector_range = str(lightcurves[0][0])+' - '+str(lightcurves[-1][0])
    for data in lightcurves:
        start = data[3].time.value[0]
        end = data[3].time.value[-1]
        period_err = CalculateErrors(data[3].time.value,data[3].flux.value,drop_in_flux,end-start)
        df =  pd.DataFrame(columns = ['Sector', 
                                      'Period',
                                      'Period Error',
                                      'Drop in Flux',
                                      'Drop in Flux Error'])
        if mode == 1:
            x_fit = [x for x in xnew if start<x<end]
            y_fit = [ynew[ind] for ind in range(len(ynew)) if start<xnew[ind]<end]
            x_peaks = [xnew[peak] for peak in peaks if start<xnew[peak]<end]
            y_peaks = [ynew[peak] for peak in peaks if start<xnew[peak]<end]
            
            if len(x_peaks) == 0:
                continue
            
            if pg == None:
                pg_data = [['No data avaiable'],None] 
            else:
                pg_data = [
                           [pg.period.value,pg.power.value,'line','darkblue','Period (days)','Power','Periodogram for sectors '+sector_range],
                           [pg.period_at_max_power.value,pg.max_power.value,'x','red']
                          ]
            
            df = df.append({'Sector': data[0],
                            'Period': period,
                            'Period Error': period_err,
                            'Drop in Flux': drop_in_flux,
                            'Drop in Flux Error': drop_in_flux_err,
                            'lc graph':[
                                        [data[3].time.value,data[3].flux.value,'line','darkblue','Time (days)','Normalised Flux','Lightcurve for sector '+str(data[0])],
                                        [x_peaks,y_peaks,'x','red'],
                                        [x_fit,y_fit,'line','lime'],
                                        [[start,end],[data[2],data[2]],'line','blue'],
                                        [data[4],data[5],'line','red']
                                       ],
                            'lc fold':[
                                       [drop_in_flux_graphs[0],drop_in_flux_graphs[1],'line','darkblue','Phase','Normalised Flux','Folded lightcurve for sectors '+sector_range],
                                       [drop_in_flux_graphs[2],drop_in_flux_graphs[3],'x','red'],
                                       [drop_in_flux_graphs[4],drop_in_flux_graphs[5],'line','lime']
                                      ],
                            'periodogram': pg_data
                            },ignore_index=True)
        else:
            df = df.append({'Sector': data[0],
                            'Period': period,
                            'Period Error': period_err,
                            'Drop in Flux': drop_in_flux,
                            'Drop in Flux Error': drop_in_flux_err
                            },ignore_index=True) 
        
        df_list.append(df)
    return

#---------------------------------------------------------------------------------------------------------------#
"""GROUPING AND AVERAGING"""
#---------------------------------------------------------------------------------------------------------------#

def GroupExoplanets(df_list):
    Possible_Exoplanets = [] 
    if len(df_list) == 1:
        return [df_list[0]]
    
    periods = [df['Period'].iloc[0] for df in df_list]
    for period in periods:
        if period == 0:
            continue
        grouped_df = pd.DataFrame(columns=['Sector', 
                                            'Period',
                                            'Period Error',
                                            'Drop in Flux',
                                            'Drop in Flux Error'
                                            ])
        Upperlim = period*1.05
        Lowerlim = period*0.95
        matches = [ind for ind in range(len(periods)) if Lowerlim<periods[ind]<Upperlim and periods[ind] != 0]
        
        for match in matches:
            grouped_df = grouped_df.append(df_list[match].iloc[0])
            periods[match] = 0
            
        grouped_df = grouped_df.reset_index(drop= True)
        Possible_Exoplanets.append(grouped_df)                
    
    return Possible_Exoplanets

def AverageValues(Exo,Star,mode):
    global refined_outliers
    results_list = []
    if len(Exo) == 0:
        return [pd.DataFrame([{'Star': Star,'Period': 0,'Period Error':0,'Drop in Flux':0,'Drop in Flux Error':0}])]
    
    for i in range(len(Exo)):  
        results =  pd.DataFrame(columns = ['Star', 
                                            'Period',
                                            'Period Error'
                                            'Drop in Flux'
                                            'Drop in Flux Error'])
        P = np.mean(Exo[i]['Period'])
        P_err = np.sqrt(sum(np.square(Exo[i]['Period Error'])))
        
        Drops = np.array(list(dict.fromkeys(np.array(Exo[i]['Drop in Flux']))))
        Drop_errs = np.array(list(dict.fromkeys(np.array(Exo[i]['Drop in Flux Error']))))
        
        if len(Drops) > 2:
            refined_outliers = False
            not_outlier, drop_std, drop_mean = RemoveOutliers(Drops)
            indices = [ind for ind in range(len(Drops)) if Drops[ind] in not_outlier]
            D = np.mean(not_outlier)
            if len(indices) > 1:
                D_err = np.std(Drop_errs[indices])/np.sqrt(len(indices))
            else:
                D_err = Drop_errs[indices][0]
        elif len(Drops) > 1:
            D = np.mean(Drops)
            D_err =  np.std(Drop_errs)/np.sqrt(len(Drop_errs))
        else:
            D = Drops[0]
            D_err = Drop_errs[0]
        
        print(P,P_err,D,D_err)
        results = results.append({'Star': Star,'Period': P,'Period Error': P_err,'Drop in Flux': D,
                                          'Drop in Flux Error': D_err}, ignore_index=True)
        results_list.append(results)
    
    return results_list     

#---------------------------------------------------------------------------------------------------------------#
"""ERROR CALCULATIONS"""
#---------------------------------------------------------------------------------------------------------------#

def CalculateErrors(time,flux,drop_in_flux,T):
    N=len(flux)
    sq_residual_flux = [abs(y-1)**2 for y in flux]
    RMS = np.sqrt(sum(sq_residual_flux)/N)
    
    if drop_in_flux != 0:
        period_err = np.sqrt(6/N)*(1/(np.pi*T))*(RMS/drop_in_flux)
    else:
        period_err = None    
        
    return period_err         

#---------------------------------------------------------------------------------------------------------------#
"""DETERMINE POSSIBLE PERIODS"""
#---------------------------------------------------------------------------------------------------------------#

def PeriodOptions(lc_flat,timespan):
    global strongest_period
    period_range = np.linspace(0.5,timespan/2,int(500*(timespan/2-0.5)))
    bls = lc_flat.to_periodogram(method='bls', period=period_range, frequency_factor=500)
    strongest_period = bls.period_at_max_power.value
    period = bls.period.value
    power = bls.power.value
    max_range = len(bls.power.value)
    
    x_cutoff, y_cutoff, _, cutoff_std = MovingAverage(period, 
                                                      power,
                                                      int(max_range/40))
    
    y_cutoff = [y_cutoff[y]+2*cutoff_std[y] for y in range(len(y_cutoff))]
    cutoff_mean = np.mean(y_cutoff)
    x_cutoff = np.insert(x_cutoff,[0,len(x_cutoff)],[period[0],period[-1]])
    y_cutoff = np.insert(y_cutoff,[0,len(y_cutoff)],[cutoff_mean,cutoff_mean])
    
    bls_peaks, _ = find_peaks(power)
    
    refined_peaks = []
    for ind in range(len(x_cutoff)-1):
        gradient = (y_cutoff[ind+1]-y_cutoff[ind])/(x_cutoff[ind+1]-x_cutoff[ind])
        intercept = y_cutoff[ind]-gradient*x_cutoff[ind]
        section_indices = np.asarray(np.where((x_cutoff[ind] < period)&(period <= x_cutoff[ind+1])))[0]
        
        peaks_in_bin = np.array([val for val in bls_peaks if x_cutoff[ind]<period[val]<x_cutoff[ind+1]])
        
        pause = False
        for point in section_indices:
            if (gradient*period[point]+intercept > power[point] or point == section_indices[-1]) and pause == False:
                current_peaks = peaks_in_bin[peaks_in_bin<point]
                if len(current_peaks) != 0:
                    refined_peaks.append(current_peaks[np.argmax([power[x] for x in current_peaks])])
                pause = True
                peaks_in_bin = peaks_in_bin[peaks_in_bin>point]
            elif gradient*period[point]+intercept < power[point] and pause == True:
                pause = False
            elif point in peaks_in_bin and gradient*period[point]+intercept > power[point]:
                peaks_in_bin = peaks_in_bin[peaks_in_bin>point]
    
    period_options = period[refined_peaks]
    power_options = [x for x in power[refined_peaks]]
    power_power_order = sorted(power_options,reverse=True)
    period_power_order = [period_options[list(power[refined_peaks]).index(x)] for x in power_power_order]
    
    for target in period_power_order:
        if target == 0:
            continue
        possible_multiples = [list(period_options).index(x) for x in period_options if x!= 0 and CheckIfSimilar(x,target) == True]
        powers = np.array(power_options)[possible_multiples]
        for i in range(1,len(powers)):
            if powers[i] < powers[i-1] and target != period_options[possible_multiples[i]]:
                period_options[possible_multiples[i]] = 0
    
    period_options = np.array(period_options)
    period_options = period_options[period_options!=0]            
    
    bls = None
    period = None
    power = None
    
    return list(sorted(period_options))

#---------------------------------------------------------------------------------------------------------------#
"""DETERMINE POSSIBLE TRANSITS"""
#---------------------------------------------------------------------------------------------------------------#

def GroupTransitPoints(flux,time,peaks):
    Transits = []
    if len(peaks) != 0:
        min_flux = min(flux[peaks]) - np.std(flux)
        bins_list = np.arange(min_flux,1,(1-min_flux)/10)
        
        counts, bin_edges = np.histogram(flux[peaks], bins=bins_list)
        
        counts = np.array(counts)
        
        centre_points = np.diff(bin_edges)/2
        centre_points = np.array([bin_edges[i]+centre_points[i] for i in range(len(bin_edges)-1)])
        step = bin_edges[1]-bin_edges[0]
        centre_points = np.insert(centre_points,[0,len(centre_points)],[centre_points[0]-step,centre_points[-1]+step])
        counts = np.insert(counts,[0,len(counts)],[0,0])
        
        histpeaks, _ = find_peaks((counts)) 
        
        start = []
        end = []
        for peak in histpeaks:
            limits = []
            direction = [-1,1]
            for opt in direction:
                counter = 0
                prev = peak
                while True:
                    counter +=1
                    current = peak+(counter*opt)
                    if current == 0 or current == len(centre_points)-1 or counts[current] == 0 or (counts[current] == counts[prev] and prev == peak):
                        limits.append(centre_points[current])
                        break
                    elif counts[current] > counts[prev] or (counts[current] == counts[prev] and prev != peak):
                        limits.append(centre_points[prev])
                        break
                    else:
                        prev = current
            start.append(limits[0])
            end.append(limits[1])
       
        for i in range(len(histpeaks)):
            peaks_ind = peaks[np.where(np.logical_and(flux[peaks]>=start[i], flux[peaks]<=end[i]))[0]]
            Transits.append([peaks_ind,time[peaks_ind],flux[peaks_ind],start[i],end[i]])
        
    return Transits    

def CheckPeaks(time,flux,peaks):
    corrected_peaks = []
    for i in range(len(peaks)):
        if peaks[i] == None:
                continue
        peaks_list = [peaks[i]]
        counter = 0
        while True:
            counter += 1
            if flux[peaks[i]+counter] > 1 or peaks[i]+counter == len(flux)-1:
                break
            if peaks[i]+counter in peaks:
                peaks_list.append(peaks[i]+counter)
                peaks = np.where(peaks==peaks[i]+counter, None, peaks) 
        
        corrected_peaks.append(peaks_list[np.array([flux[x] for x in peaks_list]).argmin()])

    return np.array(corrected_peaks)

#---------------------------------------------------------------------------------------------------------------#
"""MAKE POSSIBLE TRANSITS TREE"""
#---------------------------------------------------------------------------------------------------------------#
           
def BuildTree(peaks, transits, fluxes, periods, timespan):
    global repeat, max_period, branch, root, for_later
    None_buffer = [None]
    differences = []
    
    for i in transits:
        branch_list = None_buffer + [x-i for x in transits if x-i==abs(x-i) and x-i != 0]
        differences.append(branch_list)
        None_buffer.append(None)
   
    match_list = []
    for period in periods:
        root = []
        for_later = []
        for k in range(len(differences)):
            repeat = 1
            max_period = timespan
            branch = [k]
            FindOthers(differences,period,period,k)
            if len(branch) > 1:
                root.append(branch)
        
        counter = 0
        if len(for_later) != 0:
            while(True):
                repeat = 1
                branch = for_later[counter]
                FindOthers(differences,period,period,branch[-1])
                if len(branch) > 1:
                    root.append(branch)
                counter +=1
                if counter == len(for_later):
                    break
        
        root = sorted(root, key=len, reverse=True)
        unique = []
        for branch in range(len(root)):
            if root[branch] == [None]:
                continue
            unique.append(root[branch])
            root = CheckForOverlap(root[branch],root)
            root[branch] = [None]
            
        compare_points = []
        compare_gaps = []
        compare_est = []
        compare_flux = []
        for leaf in unique:
            correct_transits, gaps, Period_diff, mean_flux = MakeComparables(period,transits[leaf],fluxes[leaf])
            compare_points.append(correct_transits)
            compare_gaps.append(gaps)
            compare_est.append(Period_diff)
            compare_flux.append(mean_flux)
        
        correct = FindBestMatch(compare_points,compare_gaps,compare_est,compare_flux,None)
        if correct != None:
            match_list.append([period,
                               unique[correct],
                               compare_points[correct],
                               compare_gaps[correct],
                               compare_est[correct],
                               compare_flux[correct],
                               peaks]) 
      
    return match_list

def MakeComparables(period,x_peaks,y_peaks):
    correct_transits = [x for x in np.diff(x_peaks) if x<period*1.05]
    gaps = sum([round(x/period)-1 for x in np.diff(x_peaks) if x>period*1.05])
    Period_diff = abs(period - np.mean(correct_transits))
    correct_transits = len(correct_transits)
    if correct_transits == 1:
        correct_transits = 0
    mean_flux = np.mean(y_peaks)
    
    return correct_transits, gaps, Period_diff, mean_flux

def FindBestMatch(compare_points,compare_gaps,compare_est, compare_flux, flags):    
    if flags != None:
        true_flags = [ind for ind in range(len(flags)) if flags[ind] == True]
        if len(true_flags) != 0:
            compare_points = [compare_points[ind] if ind in true_flags else 0 for ind in range(len(compare_points))]
            
    max_points = [ind for ind in range(len(compare_points)) if compare_points[ind] == max(compare_points) and compare_points[ind] != 0]
     
    if len(max_points) == 0:
            correct = None
    elif len(max_points) > 1:
        masked_gaps = np.array(compare_gaps)[max_points]
        min_gaps = [ind for ind in range(len(compare_gaps)) if compare_gaps[ind] == min(masked_gaps) and ind in max_points]
        if len(min_gaps) > 1:
            masked_flux = np.array(compare_flux)[min_gaps]
            min_flux = min(masked_flux)
            flux_check = [ind for ind in range(len(compare_flux)) if min_flux*0.99<compare_flux[ind]<min_flux*1.01 and ind in min_gaps]
            if len(flux_check) > 1:
                compare_est = np.array(compare_est)[flux_check]
                correct = np.array(compare_est).argmin()
            else:
                correct = flux_check[0]
        else:
            correct = min_gaps[0]
    else:
        correct = max_points[0]
        
    return correct
      

def FindOthers(differences,og_period,period,row):
    global repeat, max_period, for_later, branch
    if period > max_period:
        return
    else:
        matching = np.array([None if x==None else x if WithinRange(x, period)==True else None for x in differences[row]])
        matching_transit = np.where(matching!=None)[0]
        if len(matching_transit) != 0:
            if len(matching_transit) > 1: 
                new_branches = [branch+[x] for x in matching_transit[1:]]
                for_later.extend(new_branches)
            if repeat != 1: 
                period = period/repeat
                repeat = 1
            branch.append(matching_transit[0])
            FindOthers(differences, og_period, period, matching_transit[0])
        else:
            repeat += 1
            FindOthers(differences, og_period, og_period*repeat, row)
            
def CheckForOverlap(branch,root):
    for ind in range(len(root)):
        check = root[ind]
        if check == [None] or branch == [None]:
            continue
        if check == branch:
            root[ind] = [None]
            continue
        if len(check) > len(branch):
            continue
        for i in range(len(branch)-len(check)+1):
            newlist = branch[:i]+check+branch[i+len(check):]
            if newlist == branch:
                root[ind] = [None]
                break
    return root
 
def AboveExpected(batch,timespan,xnew,ynew,time_range):
    compare_periods = []
    compare_peaks = []
    all_data = []
    for lot in batch:
        expected = round(timespan/lot[0])+1
        if len(lot[1]) > expected*0.5:
            compare_periods.append(lot[0])
            compare_peaks.append(lot[1])
            lot.append(True)
            all_data.append(lot)
            continue
        
        for i in range(len(time_range)):
            sector = time_range[i]
            segmented_list = []
            expected = round((sector[1]-sector[0])/lot[0])+1
            peaks_in_sector = [peak for peak in lot[1] if sector[0]<xnew[lot[6][peak]]<sector[1]]
            if len(peaks_in_sector) > expected*0.5:
                segmented_list.append([i,peaks_in_sector])
        
        if len(segmented_list) > 1:
            if len(segmented_list) >= len(time_range)/2:
                flag = True
            else:
                flag = False
            con_sectors = FindConsecutive([[x[0],None,None,x[1]] for x in segmented_list])
        
            total_correct = 0
            total_gaps = len(con_sectors) - 1
            total_diff = []
            total_mean = []
            for sec in con_sectors:
                peaks_in_sector = np.concatenate(sec[3])
                correct_transits, gaps, Period_diff, mean_flux = MakeComparables(lot[0],xnew[lot[6][peaks_in_sector]],ynew[lot[6][peaks_in_sector]])
                
                total_correct = total_correct + correct_transits
                total_gaps = total_gaps + gaps
                total_diff.append(Period_diff)
                total_mean.append(mean_flux)
            
            compare_periods.append(lot[0])
            total_peaks = np.concatenate([x[1] for x in segmented_list])
            compare_peaks.append(total_peaks)
            all_data.append([lot[0],total_peaks,total_correct,total_gaps,np.mean(total_diff),np.mean(total_mean),lot[6],flag])
            
        elif len(segmented_list) != 0:
            peaks_in_sector = segmented_list[0][1]
            correct_transits, gaps, Period_diff, mean_flux = MakeComparables(lot[0],xnew[lot[6][peaks_in_sector]],ynew[lot[6][peaks_in_sector]])
            compare_periods.append(lot[0])
            compare_peaks.append(peaks_in_sector)
            all_data.append([lot[0],peaks_in_sector,correct_transits,gaps,Period_diff,mean_flux,lot[6],False])
          
    for ind in range(len(compare_periods)-1):
        target = compare_periods[ind]
        if target == 0:
            continue            
        for jnd in range(len(compare_periods)):
            if ind == jnd or compare_periods[jnd] == 0:
                continue
            comparables = [compare_peaks[ind],compare_peaks[jnd]]
            best = FindMostLikely([all_data[ind],all_data[jnd]])
            if best == None:
                lengths = np.array([len(compare_peaks[ind]),len(compare_peaks[jnd])])
                best = lengths.argmax()
            
            root = comparables[best]
            leaf = comparables[1-best]
            matching_peaks = [x for x in leaf if x in root]
            
            if len(matching_peaks) != 0 and best == 0:
                compare_periods[jnd] = 0
            if len(matching_peaks) != 0 and best == 1:
                compare_periods[ind] = 0
    
    return compare_periods, compare_peaks, all_data

def CheckValidity(matches_per_batch,all_peaks,xnew,ynew,timespan,time_range):
    global strongest_period
    overall = []
    for batch in matches_per_batch:
        compare_periods, compare_peaks, all_data = AboveExpected(batch,timespan,xnew,ynew,time_range)   
        overall.extend([all_data[lot] for lot in range(len(all_data)) if compare_periods[lot] != 0])
        
    target_periods = [x[0] for x in overall]
    if strongest_period not in target_periods and len(all_peaks) != 0:
        strongest_list = BuildTree(all_peaks,xnew[all_peaks],ynew[all_peaks],[strongest_period],timespan)
        compare_periods, compare_peaks, all_data = AboveExpected(strongest_list,timespan,xnew,ynew,time_range)
        overall.extend([all_data[lot] for lot in range(len(all_data)) if compare_periods[lot] != 0])
        target_periods = [x[0] for x in overall]
    
    for target in target_periods:
        if target == 0:
            continue
        counts = [ind for ind in range(len(target_periods)) if target_periods[ind] == target]
        
        if len(counts) > 1:
            overall = np.array(overall,dtype=list)
            correct = FindMostLikely(overall[counts])
            if correct != None:
                del counts[correct]
                target_periods = [0 if ind in counts else target_periods[ind] for ind in range(len(target_periods))]
         
    overall = [overall[ind] for ind in range(len(overall)) if target_periods[ind] != 0]
    period_list = [x[0] for x in overall]
    peak_mask_list = [x[1] for x in overall]
    peaks = [x[6] for x in overall]
    
    return period_list, peak_mask_list, peaks

def FindMostLikely(comparables):
    compare_points = [x[2] for x in comparables]
    compare_gaps = [x[3] for x in comparables]
    compare_est = [x[4] for x in comparables]
    compare_flux = [x[5] for x in comparables]
    flags = [x[7] for x in comparables]
    best = FindBestMatch(compare_points, compare_gaps, compare_est, compare_flux, flags)
    
    return best

def CheckIfSimilar(x,target):
    multiple = round(x/target)
    upper = (x/target)*1.02
    lower = (x/target)*0.98
    if lower<multiple<upper:
        return True
    elif lower<multiple-0.5<upper:
        return True
    elif lower<multiple+0.5<upper:
        return True
    else:
        return False
    
def WithinRange(x,period):
    upper = x*1.05
    lower = x*0.95
    if period - lower > 0.5: lower = period - 0.5
    if upper - period > 0.5: upper = period + 0.5
    
    if lower<period<upper:
        return True
    else:
        return False
          
#---------------------------------------------------------------------------------------------------------------#
"""EXCLUSION FUNCTIONS"""
#---------------------------------------------------------------------------------------------------------------#

def BasicClean(flux):
    global noisemask, fluxstd
    flux[np.isnan(flux)] = 1
    
    fluxmean = np.mean(flux[(flux < 1) & (flux > 0.9)])
    fluxstd = np.std(flux[(flux < 1) & (flux > 0.9)])
    noisemask = fluxmean - fluxstd
    
    flux[flux > fluxmean + 2*fluxstd] = 1
    
    return

def RemoveExtremeOutliers(time,flux,threshold):   
    global fluxstd
    outliers = [[x,None,None] for x in range(len(flux)) if flux[x] < threshold]
    con_outliers = FindConsecutive(outliers)
    for i in range(len(con_outliers)):
        current_outliers = [x[0] for x in con_outliers[i]]
        StepBack(time,flux,current_outliers,1)
    if flux[0] < 1-5*fluxstd:
        StepBack(time,flux,[0,0],1)
    if flux[-1] < 1-5*fluxstd:
        StepBack(time,flux,[len(flux)-1,len(flux)-1],1)    
    return

def Exclude(lc):
    global noisemask 
    time = lc.time.value
    flux = lc.flux.value
    xnew, ynew, peaks, _ = MovingAverage(time, flux, 3)

    if len(peaks) == 0:
        flux[exclusion] = 1
        BasicClean(flux)
        xnew, ynew, peaks, _ = MovingAverage(time, flux, 3)
    
    Transits = GroupTransitPoints(ynew, xnew, peaks)
    ends = [batch[4] for batch in Transits]
    At_noise_floor = np.where(ends>noisemask*0.995)[0]
    for i in list(At_noise_floor):
        Transits[i] = None
            
    for batch in Transits:
        if batch == None:
            continue
        peaks = batch[0]
        start = round(time[0])
        end = round(time[-1])
        step = (end-start)/len(peaks)
        bins_list = np.arange(start,end,step)
        bins_list = np.append(bins_list,end)
        counts, bin_edges = np.histogram(xnew[peaks], bins=bins_list)
        
        mask_list = np.where(counts > 2)[0]
        
        for i in range(len(mask_list)):
            mask = (time > bins_list[mask_list[i]]) & (time < bins_list[mask_list[i]+1]) & (flux < noisemask)
            mask_index = np.where(mask==True)[0]
            StepBack(time, flux, [mask_index[0],mask_index[-1]], 1)
        
    x_diff = [x for x in np.diff(time)]
    x_diff = np.append(x_diff,0)
    x_space = 10*np.mean(x_diff)

    prev_diff = 0
    for j in range(len(time)):
        if (x_diff[j] > x_space or prev_diff > x_space) and flux[j] < noisemask:
           StepBack(time, flux, [j,j], 1)
        else:
            prev_diff = x_diff[j]
   
    return lc

def FindConsecutive(thelist):
    consec = []
    constep = []
    noisemasks = []
    last = False
    if len(thelist) == 1: 
        consec.append(thelist)
        return consec
    for i in range(len(thelist)-1):
        if thelist[i][2] == None:
            new_masks = []
            range_check = True
        else:
            new_masks = [mask for mask in [thelist[i][2],thelist[i+1][2]] if mask not in noisemasks]
            lowest = min(noisemasks+new_masks)
            highest = max(new_masks)
            if lowest not in new_masks:
                range_check = [True if lowest*0.99<new_masks[0]<lowest*1.01 else False][0]
            else:
                range_check = [True if lowest*0.99<highest<lowest*1.01 else False][0]
            
        if thelist[i][0] + 1 == thelist[i+1][0] and thelist[i][1] == thelist[i+1][1] and range_check == True:
            constep.append(thelist[i])
            noisemasks.extend(new_masks)
            last = True
        elif last == True:
            constep.append(thelist[i])
            consec.append(constep)
            constep = []
            noisemasks = []
            last = False
        else:
            consec.append([thelist[i]])
    if last == True:
        constep.append(thelist[-1])
        consec.append(constep)
        
    elif len(thelist) != 0: 
        consec.append([thelist[-1]])
    
    return consec

#---------------------------------------------------------------------------------------------------------------#
"""FITTING FUNCTIONS"""
#---------------------------------------------------------------------------------------------------------------#

def RemoveOutliers(values):
    global refined_outliers
    std = np.std(values)
    mean = np.mean(values)
    distance_from_mean = abs(values - mean)
    no_outliers = distance_from_mean <= std
    not_outlier = values[no_outliers]
    
    if len(not_outlier) <= 2:
        return not_outlier, std, mean
    if std > 0.5 and refined_outliers == False:
        refined_outliers = True
        not_outlier, std, mean = RemoveOutliers(not_outlier)
    
    return not_outlier, std, mean
    
def DropInFlux(lc_forfolding,start_time,Period,xfit_peaks,yfit_peaks):
    global refined_outliers
    num_peaks = len(xfit_peaks)-1
    
    refined_outliers = False
    not_outlier, peak_flux_std, peak_flux_mean = RemoveOutliers(yfit_peaks)
    Drop_in_Flux = np.mean(not_outlier)
    Drop_in_Flux_err = peak_flux_std/np.sqrt(len(not_outlier))
    
    lc_fold = lc_forfolding.fold(Period, epoch_time=start_time)
    lc_fold = lc_fold[(-1<lc_fold.phase.value)&(lc_fold.phase.value<1)]
   
    xnew = np.array([])
    ynew = np.array([])
    for i in range(0,len(lc_fold.flux.value)-num_peaks,1+num_peaks):
        ynew = np.append(ynew,np.mean(lc_fold.flux.value[i:i+num_peaks]))
        xnew = np.append(xnew,np.mean(lc_fold.phase.value[i:i+num_peaks]))
        
    min_point = np.where(ynew==min(ynew))[0][0]
    difference = ynew[min_point] - Drop_in_Flux
    fold_mask = ynew[min_point]+2*np.std(ynew)
    ynew[ynew<fold_mask] = ynew[ynew<fold_mask] - difference
    
    bins_list = np.arange(-1,1.1,0.1)
    counts, bin_edges = np.histogram(xnew[ynew<fold_mask], bins=bins_list)
    consec = []
    total_drops = []
    new_group = True
    for ind in range(len(counts)):
        if counts[ind] != 0:
            consec.append([ind,counts[ind]])
            new_group = False
        if (counts[ind] == 0 or ind == 18) and new_group == False:
            total_drops.append(consec)
            consec = []
            new_group = True          
    
    if len(total_drops) > 1 or len(total_drops[0]) > 10:
        return False, None, None
    else:
        max_counts = max([x[1] for x in total_drops[0]])
        index = [x[0] for x in total_drops[0] if x[1] == max_counts][0]
        if index < 2 or index > 17:
            return False, None, None
          
    return 1 - Drop_in_Flux, Drop_in_Flux_err, [lc_fold.phase.value,lc_fold.flux.value,xnew[min_point],Drop_in_Flux,xnew,ynew]

def StepBack(time,flux,outlier,target):
    global exclusion
    lower = outlier[0]
    upper = outlier[-1]
    counter = 0
    while True:
        counter +=1
        lowerlim = lower-counter
        if lowerlim <= 0:
            lowerlim = 0
            counter = 0
            break
        elif flux[lowerlim] > target:
            counter = 0 
            break
    while True:
        counter +=1
        upperlim = upper+counter
        if upperlim >= len(time)-1:
            upperlim = len(time)-1
            break
        if flux[upperlim] > target:
            break
    
    if lowerlim == 0: 
        mask = (time >= time[lowerlim]) & (time < time[upperlim])
    elif upperlim == len(time)-1: 
        mask = (time > time[lowerlim]) & (time <= time[upperlim])
    else:
        mask = (time > time[lowerlim]) & (time < time[upperlim])
        
    exclusion = [mask[i] if mask[i]==True else exclusion[i] for i in range(len(mask))] 
    return

def MovingAverage(time,flux,weight):
    global noisemask, exclusion
    
    xnew = np.array([])
    ynew = np.array([])
    std = np.array([])
    for i in range(0,len(flux)-weight,1+weight):
        std = np.append(std,np.std(flux[i:i+weight]))
        ynew = np.append(ynew,np.mean(flux[i:i+weight]))
        xnew = np.append(xnew,np.mean(time[i:i+weight]))
      
    peaks, _ = find_peaks((ynew*-1)+1)
    peaks = peaks[ynew[peaks]<noisemask]   
    
    refined_peaks = []
    for j in range(len(peaks)):
        peak_in_time = np.abs(time - xnew[peaks[j]]).argmin()
        counter = 0
        lowest_peak = [peak_in_time,peak_in_time]
        while True:
            counter +=1
            if peak_in_time+counter == len(time) -1:
                counter = 0
                break
            elif flux[peak_in_time+counter] > noisemask :
                counter = 0
                break
            elif flux[peak_in_time+counter] < flux[lowest_peak[0]]:
                lowest_peak[0] = peak_in_time+counter
        while True:
            counter +=1
            if peak_in_time-counter == 0:
                break
            if flux[peak_in_time-counter] > noisemask: 
                break
            elif flux[peak_in_time-counter] < flux[lowest_peak[1]]: 
                lowest_peak[1] = peak_in_time-counter 
        index = lowest_peak[(flux[lowest_peak]).argmin()]
        refined_peaks.append(index)
        
    
    xnew[peaks] = [x for x in time[refined_peaks]]
    ynew[peaks] = [y for y in flux[refined_peaks]]
    
    return xnew, ynew, peaks, std
