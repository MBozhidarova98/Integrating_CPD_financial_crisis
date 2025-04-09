# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 13:37:09 2023

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
import pandas as pd
import itertools

#Simulate SIR model on the multilevel network:
def multiplex_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,recovery_probabilities,tmax,initial_infecteds,return_full_data=False):
    #make the initial things:
    for edge in PMFG_layer.edges():
        PMFG_layer[edge[0]][edge[1]]['weight']=PMFG_layer[edge[0]][edge[1]]['weight']*pars[0]
    for edge in continents_layer.edges():
        continents_layer[edge[0]][edge[1]]['weight']=pars[1]
    for edge in sectors_layer.edges():
        sectors_layer[edge[0]][edge[1]]['weight']=pars[2]
    for edge in global_layer.edges():
        global_layer[edge[0]][edge[1]]['weight']=pars[3]
    t,S,I,R=[0],[len(G.nodes())-len(initial_infecteds)],[len(initial_infecteds)],[0]
    infecteds,susceptibles,recovereds=[initial_infecteds],[list(set(G.nodes)-set(initial_infecteds))],[]
    statuses=dict() #dictionary of statuses which will change on each time step
    for u in global_layer.nodes():
        statuses[u]='susceptible'
    Q=list() #queue of events
    for u in global_layer.nodes(): 
        if u in initial_infecteds:
            statuses[u]='infected'
            Event={'node':u,'action':'gets infected','time':0}
            Q.append(Event)
            for k in range(1,len(recovery_probabilities)): #check when a company recovers
                if np.random.binomial(1,recovery_probabilities[k])==1:
                    recovery_time=k
                    Event={'node':u,'action':'recovers','time':recovery_time}
                    Q.append(Event)
                    break
    total_infecteds=(Counter(statuses.values()))['infected']
    time=1
    while total_infecteds!=0 and time<tmax:
        newly_infected=list()
        for u in global_layer.nodes():
            if statuses[u]=='infected':
                for v in global_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<global_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in PMFG_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<PMFG_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in continents_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<continents_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in sectors_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<sectors_layer[u][v]['weight']:
                            newly_infected.append(v)
        for v in set(newly_infected):
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

#Simulate SIR model on the multilevel network (weights happen outside function, so that we don't do it per each simulation):
def multiplex_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,recovery_probabilities,tmax,initial_infecteds,return_full_data=False):
    #make the initial things:
    t,S,I,R=[0],[len(G.nodes())-len(initial_infecteds)],[len(initial_infecteds)],[0]
    infecteds,susceptibles,recovereds=[initial_infecteds],[list(set(G.nodes)-set(initial_infecteds))],[]
    statuses=dict() #dictionary of statuses which will change on each time step
    for u in global_layer.nodes():
        statuses[u]='susceptible'
    Q=list() #queue of events
    for u in global_layer.nodes(): 
        if u in initial_infecteds:
            statuses[u]='infected'
            Event={'node':u,'action':'gets infected','time':0}
            Q.append(Event)
            for k in range(1,len(recovery_probabilities)): #check when a company recovers
                if np.random.binomial(1,recovery_probabilities[k])==1:
                    recovery_time=k
                    Event={'node':u,'action':'recovers','time':recovery_time}
                    Q.append(Event)
                    break
    total_infecteds=(Counter(statuses.values()))['infected']
    time=1
    while total_infecteds!=0 and time<tmax:
        newly_infected=list()
        for u in global_layer.nodes():
            if statuses[u]=='infected':
                for v in global_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<global_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in PMFG_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<PMFG_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in continents_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<continents_layer[u][v]['weight']:
                            newly_infected.append(v)
                for v in sectors_layer.successors(u):
                    if statuses[v]=='susceptible':
                        if random.random()<sectors_layer[u][v]['weight']:
                            newly_infected.append(v)
        for v in set(newly_infected):
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


#Faster (but equivallent) version with aggregated network:
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

def discrete_fast_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,recovery_probabilities,tmax,initial_infecteds,return_full_data=False):
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

##################Time window n: ################
#Get the estimated parameters:
df=pd.read_csv('2020_10.xls')
b1=list(df['beta1'])
b2=list(df['beta2'])
b3=list(df['beta3'])
b4=list(df['beta4'])
p=list(df['p'])

infected_companies_2020=estimate_infection.infected_companies_2020
days_2020=estimate_infection.days_2020
counted_2020=estimate_infection.counted_2020

#Predictions: #n is the window size, k is the number of future days
def plot_infections(n,k):
    plt.figure()
    for i in range(n,len(days_2020)-k-1):
        tmax=k
        iterations = 10000
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-n],b2[i-n],b3[i-n],b4[i-n]]
        for counter in range(iterations):
            #t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n:i+k],tmax,initial_infecteds,return_full_data=False)
            t,S,I,R=discrete_fast_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n]*np.ones(k),tmax,initial_infecteds,return_full_data=False)
            tshift = EoN.get_time_shift(t, I+R, 50)
            min_time_after_shift=min(np.array(t)-tshift)
            min_time_after_shifting=min(min_time_after_shifting,min_time_after_shift)
            max_time_after_shift=max(np.array(t)-tshift)
            max_time_after_shifting=max(max_time_after_shift,max_time_after_shifting)
            xs.append(list(np.array(t)-tshift))
            ys.append(list(I))
        x_values=np.linspace(min_time_after_shifting,max_time_after_shifting,1000)
        ys_interp = [np.interp(x_values, xs[i], ys[i]) for i in range(len(xs))]
        mean_y_axis = np.mean(ys_interp, axis=0)
        x=np.linspace(i,i+tmax,1000)
        plt.plot(x,mean_y_axis,linewidth=2)
    plt.plot(counted_2020.values(),color='black')
    plt.xlabel('Day')
    plt.ylabel('Infected companies')
    plt.savefig('n{}_k{}_2020_means_predictions_cleared.pdf'.format(n,k))

###### Also get the variance #####
def plot_infections(n,k):
    plt.figure()
    upper_CIs=list()
    lower_CIs=list()
    x_axis=list()
    for i in range(n,len(days_2020)-k-1,30):
        tmax=k
        iterations = 10
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-n],b2[i-n],b3[i-n],b4[i-n]]
        for counter in range(iterations):
            #t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n:i+k],tmax,initial_infecteds,return_full_data=False)
            t,S,I,R=discrete_fast_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n]*np.ones(k),tmax,initial_infecteds,return_full_data=False)
            tshift = EoN.get_time_shift(t, I+R, 50)
            min_time_after_shift=min(np.array(t)-tshift)
            min_time_after_shifting=min(min_time_after_shifting,min_time_after_shift)
            max_time_after_shift=max(np.array(t)-tshift)
            max_time_after_shifting=max(max_time_after_shift,max_time_after_shifting)
            xs.append(list(np.array(t)-tshift))
            ys.append(list(I)) 
        x_values=np.linspace(min_time_after_shifting,max_time_after_shifting,1000)
        ys_interp = [np.interp(x_values, xs[i], ys[i]) for i in range(len(xs))]
        mean_y_axis = np.mean(ys_interp, axis=0)
        upper_CI=np.quantile([y[-1] for y in ys],0.95)
        lower_CI=np.quantile([y[-1] for y in ys],0.05)
        upper_CIs.append(upper_CI)
        lower_CIs.append(lower_CI)
        x=np.linspace(i,i+tmax,1000)
        x_axis.append(x[-1])
        plt.plot(x,mean_y_axis,linewidth=2)
    plt.fill_between(x_axis, upper_CIs, lower_CIs, color='grey', alpha=0.2)
    plt.plot(counted_2020.values(),color='black')
    plt.xlabel('Day')
    plt.ylabel('Infected companies')
    plt.show()
    plt.savefig('n{}_k{}_2020_means_predictions_with_CI.pdf'.format(n,k))

n=30
k=20
df=pd.read_csv('2020_{}_cleared_NEW.xls'.format(n))
b1=list(df['beta1'])
b2=list(df['beta2'])
b3=list(df['beta3'])
b4=list(df['beta4'])
p=list(df['p'])

plot_infections(n,k)

#Get full CIs of the curves:
def plot_infections(n,k):
    plt.figure()
    for i in range(n,len(days_2020)-k-1,30):
        tmax=k
        iterations = 10
        min_time_after_shifting=0
        max_time_after_shifting=0
        xs=list()
        ys=list()
        initial_infecteds=infected_companies_2020[i]
        pars=[b1[i-n],b2[i-n],b3[i-n],b4[i-n]]
        for counter in range(iterations):
            #t,S,I,R=discrete_fast_SIR_2(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n:i+k],tmax,initial_infecteds,return_full_data=False)
            t,S,I,R=discrete_fast_SIR(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,p[i-n]*np.ones(k),tmax,initial_infecteds,return_full_data=False)
            tshift = EoN.get_time_shift(t, I+R, 50)
            min_time_after_shift=min(np.array(t)-tshift)
            min_time_after_shifting=min(min_time_after_shifting,min_time_after_shift)
            max_time_after_shift=max(np.array(t)-tshift)
            max_time_after_shifting=max(max_time_after_shift,max_time_after_shifting)
            xs.append(list(np.array(t)-tshift))
            ys.append(list(I)) 
        x_values=np.linspace(min_time_after_shifting,max_time_after_shifting,1000)
        ys_interp = [np.interp(x_values, xs[i], ys[i]) for i in range(len(xs))]
        mean_y_axis = np.mean(ys_interp, axis=0)
        upper_CI=np.quantile(ys_interp,0.95,axis=0)
        lower_CI=np.quantile(ys_interp,0.05,axis=0)
        x=np.linspace(i,i+tmax,1000)
        plt.plot(x,mean_y_axis,linewidth=2)
        plt.fill_between(x, upper_CI, lower_CI, color='grey', alpha=0.2)
    plt.plot(counted_2020.values(),color='black')
    plt.xlabel('Day')
    plt.ylabel('Infected companies')
    plt.show()