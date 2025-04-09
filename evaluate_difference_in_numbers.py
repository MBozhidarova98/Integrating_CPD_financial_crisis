# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 15:09:46 2023

@author: mbozh
"""

import random
import heapq
import numpy as np
import EoN
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import pandas as pd
import csv
def weights(pars,PMFG_layer, continents_layer, sectors_layer, global_layer):
    loc=pars[0]
    cont=pars[1]
    sec=pars[2]
    glob=pars[3]
    for edge in global_layer.edges():
        global_layer[edge[0]][edge[1]]['weight']=glob
        if edge in PMFG_layer.edges():
            global_layer[edge[0]][edge[1]]['weight']+=loc*PMFG_layer[edge[0]][edge[1]]['weight']
        if edge in continents_layer.edges():
            global_layer[edge[0]][edge[1]]['weight']+=cont
        if edge in sectors_layer.edges():
            global_layer[edge[0]][edge[1]]['weight']+=sec
    return(global_layer)

def discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,recovery_probabilities,tmax,initial_infecteds,return_full_data=False):
    #make the initial things:
    G1=weights(pars,PMFG_layer, continents_layer, sectors_layer, global_layer)
    G=nx.DiGraph()
    for e in G1.edges():
        G.add_edge(e[0],e[1],weight=G1[e[0]][e[1]]['weight'])
        G.add_edge(e[1],e[0],weight=G1[e[0]][e[1]]['weight'])
    per_edge_probabilities=dict()
    for e in G.edges():
        per_edge_probabilities[e]=G[e[0]][e[1]]['weight']
    t,S,I,R=[0],[len(G.nodes())-len(initial_infecteds)],[len(initial_infecteds)],[0]
    infecteds,susceptibles,recovereds=[initial_infecteds],[list(set(G.nodes)-set(initial_infecteds))],[]
    statuses=dict() #dictionary of statuses which will change on each time step
    for u in G.nodes():
        statuses[u]='susceptible'
    Q=list() #queue of events
    for u in G.nodes():
        if u in initial_infecteds:
            statuses[u]='infected'
            Event={'node':u,'action':'gets infected','time':0}
            Q.append(Event)
            for k in range(1,len(recovery_probabilities)):
                if np.random.binomial(1,recovery_probabilities[k])==1:
                    recovery_time=k
                    Event={'node':u,'action':'recovers','time':recovery_time}
                    Q.append(Event)
                    break
    total_infecteds=(Counter(statuses.values()))['infected']
    time=1
    while total_infecteds!=0 and time<tmax:
        for u in G.nodes():
            if statuses[u]=='infected':
                for v in G.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<per_edge_probabilities[(u,v)]:
                            Event={'node':v,'action':'gets infected','time':time}
                            Q.append(Event)
                            for k in range(time+1,len(recovery_probabilities)):
                                if np.random.binomial(1,recovery_probabilities[k])==1:
                                    recovery_time=k
                                    Event={'node':v,'action':'recovers','time':recovery_time}
                                    Q.append(Event)
                                    break
        #print(Q)
        for event in Q:
            if event['time']==time:
                if event['action']=='gets infected':
                    statuses[event['node']]='infected'
                if event['action']=='recovers':
                    statuses[event['node']]='recovered'
        t.append(time)
        S.append((Counter(statuses.values()))['susceptible'])
        I.append((Counter(statuses.values()))['infected'])
        R.append((Counter(statuses.values()))['recovered'])
        if return_full_data==True:
            susceptible,infected,recovered=[],[],[]
            for i in statuses.keys():
                if statuses[i]=='susceptible':
                    susceptible.append(i)
                if statuses[i]=='recovered':
                    recovered.append(i)
                if statuses[i]=='infected':
                    infected.append(i)
            infecteds.append(infected)
            susceptibles.append(susceptible)
            recovereds.append(recovered) 
        time=time+1
        total_infecteds=(Counter(statuses.values()))['infected']  
    if return_full_data==True:
        return(t,susceptibles,infecteds,recovereds)
    return(t,S,I,R)

#For 2020 financial crisis:
#Construct the network
import estimate_infection
l=estimate_infection.log_returns
nodes=list(l.columns)[:-1]


G=nx.Graph()
G.add_nodes_from(nodes)
edges=pd.read_csv('PMFG_chis_95_before_2020.csv')
for i in range(len(edges)):
    G.add_edge(edges.iloc[i]['Source'],edges.iloc[i]['End'],weight=edges.iloc[i]['weight'])
PMFG_layer=G

global_layer=nx.Graph()
global_layer.add_nodes_from(nodes)
global_layer.add_edges_from(itertools.permutations(nodes, 2))

sectors_layer=nx.Graph()
sectors_layer.add_nodes_from(nodes)
Sectors_edges=pd.read_csv('Sectors_network.csv')
source=list(Sectors_edges['Company1'])
target=list(Sectors_edges['Company2'])
for i in range(len(Sectors_edges)):
    sectors_layer.add_edge(source[i],target[i])

continents_layer=nx.Graph()
continents_layer.add_nodes_from(nodes)
Continents_edges=pd.read_csv('Continents_network.csv')
source=list(Continents_edges['Company1'])
target=list(Continents_edges['Company2'])
for i in range(len(Continents_edges)):
    continents_layer.add_edge(source[i],target[i])

infected_companies_2020=estimate_infection.infected_companies_2020
days_2020=estimate_infection.days_2020
counted_2020=estimate_infection.counted_2020

#N is number of future days we simulate for 
def evaluate_recoveries(N,number_previous_days):
    absolute_difference_in_numbers=list()
    df=pd.read_csv('2020_{}.xls'.format(number_previous_days))
    b1=list(df['beta1'])
    b2=list(df['beta2'])
    b3=list(df['beta3'])
    b4=list(df['beta4'])
    p=list(df['p'])
    for i in range(number_previous_days,len(days_2020)-N-1):
        newly_recovered=len(set(infected_companies_2020[i])-set(infected_companies_2020[i+N]))
        tmax=N
        iterations = 10000
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-number_previous_days],b2[i-number_previous_days],b3[i-number_previous_days],b4[i-number_previous_days]]
        a=list()
        for counter in range(iterations):
            t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-number_previous_days:i+N],tmax,initial_infecteds,return_full_data=True)
            predicted_newly_recovered=len(set(initial_infecteds)-set(I[-1]))
            a.append(np.abs(predicted_newly_recovered-newly_recovered))
        absolute_difference_in_numbers.append(np.mean(a))
    return(absolute_difference_in_numbers) 

for n in [1,2,3,4,5,6,7,8,9,10,20,30]:
    num=list()
    for k in [1,10,20,30]:
        num.append(evaluate_recoveries(k,n))
    with open("data_2020_{}_absolute_difference_newly_recovered.csv".format(n), "w") as f:
        wr = csv.writer(f)
        wr.writerows(num)

def evaluate_total_numbers(N,number_previous_days):
    absolute_difference_in_numbers=list()
    df=pd.read_csv('2020_{}.xls'.format(number_previous_days))
    b1=list(df['beta1'])
    b2=list(df['beta2'])
    b3=list(df['beta3'])
    b4=list(df['beta4'])
    p=list(df['p'])
    for i in range(number_previous_days,len(days_2020)-N-1):
        tmax=N
        iterations = 10000
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-number_previous_days],b2[i-number_previous_days],b3[i-number_previous_days],b4[i-number_previous_days]]
        common_continents=list()
        for counter in range(iterations):
            t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-number_previous_days:i+N],tmax,initial_infecteds,return_full_data=False)
            tshift = EoN.get_time_shift(t, I+R, 50)
            min_time_after_shift=min(np.array(t)-tshift)
            min_time_after_shifting=min(min_time_after_shifting,min_time_after_shift)
            max_time_after_shift=max(np.array(t)-tshift)
            max_time_after_shifting=max(max_time_after_shift,max_time_after_shifting)
            #Find the mean of the shifted curves:
            xs.append(list(np.array(t)-tshift))
            ys.append(list(I))
        x_values=np.linspace(min_time_after_shifting,max_time_after_shifting,1000)
        ys_interp = [np.interp(x_values, xs[i], ys[i]) for i in range(len(xs))]
        mean_y_axis = np.mean(ys_interp, axis=0)
        absolute_difference_in_numbers.append(np.abs(mean_y_axis[-1]-list(counted_2020.values())[i+N]))
    return(absolute_difference_in_numbers)   

for n in [1,2,3,4,5,6,7,8,9,10,20,30]:
    num=list()
    for N in [1,10,20,30]:
        num.append(evaluate_total_numbers(k,n))
    with open("data_2020_{}_absolute_difference_total_numbers.csv".format(n), "w") as f:
        wr = csv.writer(f)
        wr.writerows(num)


def evaluate_newly_infected(N,number_previous_days):
    absolute_difference_in_numbers=list()
    df=pd.read_csv('2020_{}.xls'.format(number_previous_days))
    b1=list(df['beta1'])
    b2=list(df['beta2'])
    b3=list(df['beta3'])
    b4=list(df['beta4'])
    p=list(df['p'])
    for i in range(number_previous_days,len(days_2020)-N-1):
        newly_infected=len(set(infected_companies_2020[i+N])-set(infected_companies_2020[i]))
        tmax=N
        iterations = 10000
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-number_previous_days],b2[i-number_previous_days],b3[i-number_previous_days],b4[i-number_previous_days]]
        a=list()
        for counter in range(iterations):
            t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-number_previous_days:i+N],tmax,initial_infecteds,return_full_data=True)
            predicted_newly_infected=len(set(I[-1])-set(initial_infecteds))
            a.append(np.abs(predicted_newly_infected-newly_infected))
        absolute_difference_in_numbers.append(np.mean(a))
    return(absolute_difference_in_numbers)  

for i in [1,2,3,4,5,6,7,8,9,10,20,30]:
    num=list()
    for N in [1,10,20,30]:
        num.append(evaluate_newly_infected(k,n))
    with open("data_2020_{}_absolute_difference_newly_infected.csv".format(n), "w") as f:
        wr = csv.writer(f)
        wr.writerows(num)
        
