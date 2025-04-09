# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 10:44:47 2024

@author: pmxmb14
"""

import get_cleared_data
#pip install changefinder
import changefinder

############ Find the local maximas in the scores:
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def find_significant_local_maxima(time_series, threshold, min_distance=1, prominence=2):
    # Add prominence filtering
    peaks, properties = find_peaks(time_series, distance=min_distance, prominence=prominence)
    
    significant_peaks = [peak for peak in peaks if time_series[peak] > threshold]
    peak_values = [time_series[peak] for peak in significant_peaks]
    
    return significant_peaks, peak_values

#For 2008 financial crisis:
counted_2008=get_cleared_data.counted_2008
points=list(counted_2008.values())

cf = changefinder.ChangeFinder()
scores = [cf.update(p) for p in points]

time_series=scores
threshold=10
significant_peaks, peak_values = find_significant_local_maxima(time_series, threshold)
largest_scores = significant_peaks
#Find the changepoint here:
f, (ax1, ax2) = plt.subplots(2, 1)
f.subplots_adjust(hspace=0.4)
ax1.plot(counted_2008.keys(),counted_2008.values(),color='black')
for i in largest_scores:
    ax1.axvline(x=list(counted_2008.keys())[i],color='red')
#add the important events:
ax1.axvline(x=datetime(2007,8,1,0,0),color='orange',linestyle='-',label='Interbank market \nfreeze')
ax1.axvline(x=datetime(2008,3,15,0,0),color='blue',linestyle='--',label='Bear Stearns \ncollapse')
ax1.axvline(x=datetime(2008,9,6,0,0),color='red',linestyle=':',label='Lehman Brothers \nbankruptcy')
ax1.axvline(x=datetime(2008,11,3,0,0),color='purple',dashes=[3, 5, 2, 5],label='TARP')
ax1.axvline(x=datetime(2009,1,20,0,0),color='green',linestyle='-.',label='ARRA')
ax1.legend(bbox_to_anchor=(1.05, 0.5))
ax1.set_ylabel("Infected companies")
N = 2  # 1 tick every 2
xticks_pos = ax1.get_xticks()
xticks_labels = ax1.get_xticklabels()
myticks = [j for i,j in enumerate(xticks_pos) if not i%N]  # index of selected ticks
newlabels = [label for i,label in enumerate(xticks_labels) if not i%N]
ax1.set_xticks(myticks,newlabels)
#Initiate changefinder function
cf = changefinder.ChangeFinder()
scores = [cf.update(p) for p in points]
ax2.plot(scores,color='black')
#ax2.set_title("anomaly score")
ax2.set_xlabel('Day')
ax2.set_ylabel('Anomaly score')
for i in largest_scores:
    ax2.axvline(x=i,color='red')
#ax2.set_xlim(0,500)
ax2.axvline(x=15,color='orange',linestyle='-',label='Interbank market \nfreeze')
ax2.axvline(x=250,color='blue',linestyle='--',label='Bear Stearns \ncollapse')
ax2.axvline(x=420,color='red',linestyle=':',label='Lehman Brothers \nbankruptcy')
ax2.axvline(x=475,color='purple',dashes=[3, 5, 2, 5],label='TARP')
ax2.axvline(x=550,color='green',linestyle='-.',label='ARRA')
plt.savefig('Change_point_2008_with_events.pdf',bbox_inches='tight')
plt.show()


#For 2020 financial crisis:
counted_2020=get_cleared_data.counted_2020
points=list(counted_2020.values())

cf = changefinder.ChangeFinder()
scores = [cf.update(p) for p in points]

time_series=scores
threshold=7.5
significant_peaks, peak_values = find_significant_local_maxima(time_series, threshold)

largest_scores = significant_peaks
#largest_scores=np.array([13,32,52,138,161,219,227,271,369,384,417,439,446,485]) #2020
from datetime import datetime
f, (ax1, ax2) = plt.subplots(2, 1)
f.subplots_adjust(hspace=0.4)
ax1.plot(counted_2020.keys(),counted_2020.values(),color='black')
for i in largest_scores:
    ax1.axvline(x=list(counted_2020.keys())[i],color='red')
#add the important events:
ax1.axvline(x=datetime(2020,2,25,0,0),color='orange',linestyle='-',label='Stock market crash')
ax1.axvline(x=datetime(2020,3,11,0,0),color='blue',linestyle='--',label='COVID-19 declared \na global pandemic')
ax1.axvline(x=datetime(2020,3,25,0,0),color='red',linestyle=':',label='Governments offer \nstimulus packages')
ax1.axvline(x=datetime(2020,6,20,0,0),color='purple',dashes=[3, 5, 2, 5],label='Governments ease \nlockdown restrictions')
ax1.axvline(x=datetime(2020,7,15,0,0),color='green',linestyle='-.',label='Significant increase \nin COVID-19 cases')
#plt.axvline(x=datetime(2008,9,15,0,0),color='green',label='Stock Market crashes')
ax1.legend(bbox_to_anchor=(1.05, 0.5))
ax1.set_ylabel("Infected companies")
N = 2  # 1 tick every 2
xticks_pos = ax1.get_xticks()
xticks_labels = ax1.get_xticklabels()
myticks = [j for i,j in enumerate(xticks_pos) if not i%N]  # index of selected ticks
newlabels = [label for i,label in enumerate(xticks_labels) if not i%N]
ax1.set_xticks(myticks,newlabels)
#Initiate changefinder function
cf = changefinder.ChangeFinder()
scores = [cf.update(p) for p in points]
ax2.plot(scores,color='black')
ax2.set_xlabel('Day')
ax2.set_ylabel('Anomaly score')
for i in largest_scores:
    ax2.axvline(x=i,color='red')
ax2.axvline(x=10,color='orange',linestyle='-',label='Stock market crash')
ax2.axvline(x=26,color='blue',linestyle='--',label='COVID-19 declared \na global pandemic')
ax2.axvline(x=42,color='red',linestyle=':',label='Governments offer \nstimulus packages')
ax2.axvline(x=128,color='purple',dashes=[3, 5, 2, 5],label='Governments ease \nlockdown restrictions')
ax2.axvline(x=152,color='green',linestyle='-.',label='Significant increase \nin COVID-19 cases')
plt.savefig('Change_point_2020_with_events.pdf',bbox_inches='tight')
plt.show()
