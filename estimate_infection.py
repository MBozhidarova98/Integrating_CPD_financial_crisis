# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:00:53 2023

@author: mbozh
"""

import sys
import csv
import glob
import pandas as pd
import os 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#Read the data:
df=pd.read_csv('Stock_price_data.csv')
df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m-%d')
dates=list(df['Date'])
df=df.drop('Date',axis=1)
#df=df.drop('Unnamed: 0',axis=1)

#Get log returns:
def get_log_returns(df):
    new_df=pd.DataFrame(columns=df.columns)
    for i in range(len(df)-1):
        my_array=np.array(df.iloc[i]/df.iloc[i+1]).astype(float)
        log_ret=np.log(my_array)
        new_df.loc[i] = list(log_ret)
    return(new_df)

log_returns=get_log_returns(df)
log_returns['Date']=dates[1:]

#get volatility:
def get_volatility_curve(company_name,data_frame_of_log_returns,T):
    sds=list()
    dates=list(data_frame_of_log_returns[[company_name,'Date']].dropna()['Date'])[::-1]
    date=list()
    vector_of_values=list(data_frame_of_log_returns[company_name].dropna())[::-1]
    for i in range(len(vector_of_values)-T-1):
        data=vector_of_values[i:i+T]
        sds.append(np.std(data))
        date.append(dates[i+T])
    return(sds,date)

#get mean
def get_mean_curve(company_name,data_frame_of_log_returns,T):
    means=list()
    dates=list(data_frame_of_log_returns[[company_name,'Date']].dropna()['Date'])[::-1]
    date=list()
    vector_of_values=list(data_frame_of_log_returns[company_name].dropna())[::-1]
    for i in range(len(vector_of_values)-T-1):
        data=vector_of_values[i:i+T]
        means.append(np.mean(data))
        date.append(dates[i+T])
    return(means,date)

#Get the infections: 
T=21 #set the rolling time period
threshold=0.9 #set the volatility threshold
dates_list=list()
c=(list(log_returns.columns))[:-1]
for i in c:
    volatilities=dict()
    volatility,date=get_volatility_curve(i,log_returns,T)
    for j in range(len(volatility)):
        volatilities[date[j]]=volatility[j]
    t=np.quantile(volatility,threshold)
    means=dict()
    mean,date_means=get_mean_curve(i,log_returns,T)
    for j in range(len(mean)):
        means[date_means[j]]=mean[j]
    for j in date:
        if volatilities[j]>=t  and means[j]<=0:
            dates_list.append(j) #a company is infected whenever high volatility+negative mean

#Count the number of infected companies at each day:
from collections import Counter
infected=dict(Counter(sorted(dates_list)))

#start_date = pd.to_datetime('2002-01-18', format='%Y-%m-%d')
#end_date = pd.to_datetime('2022-01-17', format='%Y-%m-%d')
#infected = {date: value for date, value in infected.items() if start_date <= date <= end_date}


plt.figure()
plt.plot(list(infected.keys()),list(infected.values()),color='black')
plt.plot(list(infected.keys())[0:450],list(infected.values())[0:450],label='Early 2000s recession')
plt.plot(list(infected.keys())[1400:2000],list(infected.values())[1400:2000],label='The Great Recession ')
plt.plot(list(infected.keys())[2300:2600],list(infected.values())[2300:2600],label='Cypriot financial crisis')
plt.plot(list(infected.keys())[3300:3800],list(infected.values())[3300:3800],label='2015-16 stock market \nselloff')
plt.plot(list(infected.keys())[4350:4800],list(infected.values())[4350:4800],label='COVID-19 recession')
plt.xlabel('Date',fontsize=14)
plt.ylabel('Infected companies',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.ylim(0,350)
plt.legend(loc='upper left',fontsize=10)
#ax = plt.gca()
#ax.legend(loc='center right',bbox_to_anchor=(1.25, 1.1),ncol=3)
plt.show()

################ Count company as infected if it is infected in consequitive days:
#####2008 crisis:              
data_with_companies=dict()
companies_with_dates=dict()
for i in c:
    volatility,date=get_volatility_curve(i,log_returns,T)
    mean,date=get_mean_curve(i,log_returns,T)
    t=np.quantile(volatility,threshold)
    for j in range(len(date)):
        if volatility[j]>=t and mean[j]<=0:
            if date[j] in data_with_companies.keys():
                data_with_companies[date[j]].append(i)
            else:
                data_with_companies[date[j]]=[i]
            if i in companies_with_dates.keys():
                companies_with_dates[i].append(date[j])
            else:
                companies_with_dates[i]=[date[j]]

#Get the start and end date for 2008 financial crisis:  
start_date=list(infected.keys())[1300]
end_date=list(infected.keys())[1930]


from datetime import date, timedelta
def get_between_dates(sdate,edate):
    delta = edate - sdate
    dates=list()       # as timedelta
    for i in range(1,delta.days):
        day = sdate + timedelta(days=i)
        dates.append(day)
    return(dates)

dates_2008=get_between_dates(start_date, end_date)
dates_2008_companies=dict()
for i in companies_with_dates.keys():
    for j in companies_with_dates[i]:
        if j in dates_2008:
            if i in dates_2008_companies.keys():
                dates_2008_companies[i].append(j)
            else:
                dates_2008_companies[i]=[j]
                

for i in dates_2008_companies.keys():
    dates_2008_companies[i]=get_between_dates(min(dates_2008_companies[i]),max(dates_2008_companies[i]))    


###### Now plot the infections:
dates_companies_2008=dict()
for i in dates_2008_companies.keys():
    for j in dates_2008_companies[i]:
        if j in dates_companies_2008.keys():
            dates_companies_2008[j].append(i)
        else:
            dates_companies_2008[j]=[i]


####Count them:
counted_2008=dict()
for i in dates_companies_2008.keys():
    counted_2008[i]=len(dates_companies_2008[i])

import collections
counted_2008 = collections.OrderedDict(sorted(counted_2008.items()))

#Plot the crisis+important events
from datetime import datetime
plt.figure()
plt.plot(counted_2008.keys(),counted_2008.values(),color='black')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Infected companies')
plt.axvline(x=datetime(2007,8,1,0,0),color='orange',linestyle='-',label='Interbank market \nfreeze')
plt.axvline(x=datetime(2008,3,15,0,0),color='blue',linestyle='--',label='Bear Stearns \ncollapse')
plt.axvline(x=datetime(2008,9,6,0,0),color='red',linestyle=':',label='Lehman Brothers \nbankruptcy')
plt.axvline(x=datetime(2008,11,3,0,0),color='purple',dashes=[3, 5, 2, 5],label='TARP')
plt.axvline(x=datetime(2009,2,1,0,0),color='green',linestyle='-.',label='ARRA')
#plt.axvline(x=datetime(2008,9,15,0,0),color='green',label='Stock Market crashes')
#plt.legend(loc='upper left',fontsize=7.5)
ax = plt.gca()
ax.legend(loc='upper right',fontsize=6.7)
#ax.legend(loc='center right',bbox_to_anchor=(1.6, 0.5),ncol=1)
#plt.savefig(2008crisis.pdf',bbox_inches='tight')


#####Do for 2020 financial crisis:
end_date=list(infected.keys())[4770]
start_date=list(infected.keys())[4410]



from datetime import date, timedelta
def get_between_dates(sdate,edate):
    delta = edate - sdate
    dates=list()       # as timedelta
    for i in range(1,delta.days):
        day = sdate + timedelta(days=i)
        dates.append(day)
    return(dates)

dates_2020=get_between_dates(start_date, end_date)
dates_2020_companies=dict()
for i in companies_with_dates.keys():
    for j in companies_with_dates[i]:
        if j in dates_2020:
            if i in dates_2020_companies.keys():
                dates_2020_companies[i].append(j)
            else:
                dates_2020_companies[i]=[j]
                

for i in dates_2020_companies.keys():
    dates_2020_companies[i]=get_between_dates(min(dates_2020_companies[i]),max(dates_2020_companies[i]))    


dates_companies_2020=dict()
for i in dates_2020_companies.keys():
    for j in dates_2020_companies[i]:
        if j in dates_companies_2020.keys():
            dates_companies_2020[j].append(i)
        else:
            dates_companies_2020[j]=[i]


####Count them:
counted_2020=dict()
for i in dates_companies_2020.keys():
    counted_2020[i]=len(dates_companies_2020[i])

import collections
counted_2020 = collections.OrderedDict(sorted(counted_2020.items()))

#Make nice plot with events:
r'''
https://www.mckinsey.com/capabilities/strategy-and-corporate-finance/our-insights/the-impact-of-covid-19-on-capital-markets-one-year-in
https://www.investopedia.com/government-stimulus-efforts-to-fight-the-covid-19-crisis-4799723
r'''
plt.figure()
plt.plot(counted_2020.keys(),counted_2020.values(),color='black')
plt.xlabel('Date')
plt.xticks(rotation=45)
plt.ylabel('Infected companies')
plt.axvline(x=datetime(2020,2,25,0,0),color='orange',linestyle='-',label='Stock market crash')
plt.axvline(x=datetime(2020,3,11,0,0),color='blue',linestyle='--',label='COVID-19 declared \na global pandemic')
plt.axvline(x=datetime(2020,3,25,0,0),color='red',linestyle=':',label='Governments offer \nstimulus packages')
plt.axvline(x=datetime(2020,6,20,0,0),color='purple',dashes=[3, 5, 2, 5],label='Governments ease \nlockdown restrictions')
plt.axvline(x=datetime(2020,7,15,0,0),color='green',linestyle='-.',label='Significant increase in \nCOVID-19 cases worldwide')
#plt.axvline(x=datetime(2008,9,15,0,0),color='green',label='Stock Market crashes')
#plt.legend(loc='upper right',fontsize=8)
ax = plt.gca()
ax.legend(loc='upper right',fontsize=7.5)
#ax.legend(loc='center right',bbox_to_anchor=(1.5, 0.5),ncol=1)
#plt.savefig('2020crisis.pdf',bbox_inches='tight')

#Get a list of infected, susceptible, recovered at each time step:
od=collections.OrderedDict(sorted(dates_companies_2008.items()))
dates_2008=list(od.keys())
days_2008=[i for i in range(1,len(dates_2008)+1)]
infected_companies_2008=list(od.values())

recovered_companies_2008 = [[]]
for i in range(len(infected_companies_2008) - 1):
    recovered_companies = set(infected_companies_2008[i]) - set(infected_companies_2008[i + 1])
    recovered_companies_2008.append(list(recovered_companies))

od=collections.OrderedDict(sorted(dates_companies_2020.items()))
dates_2020=list(od.keys())
days_2020=[i for i in range(1,len(dates_2020)+1)]
infected_companies_2020=list(od.values())

recovered_companies_2020 = [[]]
for i in range(len(infected_companies_2020) - 1):
    recovered_companies = set(infected_companies_2020[i]) - set(infected_companies_2020[i + 1])
    recovered_companies_2020.append(list(recovered_companies))
    
