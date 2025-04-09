# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 14:27:37 2023

@author: mbozh
"""

import numpy as np
import networkx as nx
import estimate_infection
import pandas as pd
import itertools

#Define the negative log-likelihood function to be optimized:
def negative_log_likelihood(pars,PMFG_layer, continents_layer, sectors_layer, global_layer,infected_companies,days):
    l=0
    loc=pars[0]
    cont=pars[1]
    sec=pars[2]
    glob=pars[3]
    recovereds=[[]]
    for i in range(len(days)-1):
        infected_nodes=set(infected_companies[i+1])-set(infected_companies[i])
        recovered_nodes=list(set(infected_companies[i])-set(infected_companies[i+1]))
        for j in recovereds[i]:
            recovered_nodes.append(j)
        recovereds.append(recovered_nodes)
        #print(recovereds)
        for node in PMFG_layer.nodes():
            prob_non_infection=1
            if node not in recovereds[i] and node not in infected_companies[i]:
                for neighbour in PMFG_layer.neighbors(node):
                    if neighbour in infected_companies[i]:
                        prob_non_infection=prob_non_infection*(1-PMFG_layer[neighbour][node]['weight']*loc)
                for neighbour in continents_layer.neighbors(node):
                    if neighbour in infected_companies[i]:
                        prob_non_infection=prob_non_infection*(1-cont)
                for neighbour in sectors_layer.neighbors(node):
                    if neighbour in infected_companies[i]:
                        prob_non_infection=prob_non_infection*(1-sec)
                for neighbour in global_layer.neighbors(node):
                    if neighbour in infected_companies[i]:
                        prob_non_infection=prob_non_infection*(1-glob)
            if node in infected_nodes:
                l=l+np.log(1-prob_non_infection)
            else:
                l=l+np.log(prob_non_infection)
    return -l

# Construct the network
l=estimate_infection.log_returns
nodes=list(l.columns)[:-1]
#Get PMFG layer:
G=nx.Graph()
G.add_nodes_from(nodes)
edges=pd.read_csv('PMFG_chis_95_before_2020.csv')
for i in range(len(edges)):
    G.add_edge(edges.iloc[i]['Source'],edges.iloc[i]['End'],weight=edges.iloc[i]['weight'])
PMFG_layer=G

#GLobal layer:
global_layer=nx.Graph()
global_layer.add_nodes_from(nodes)
global_layer.add_edges_from(itertools.permutations(nodes, 2))

#Sectors layer:
sectors_layer=nx.Graph()
sectors_layer.add_nodes_from(nodes)
Sectors_edges=pd.read_csv('Sectors_network.csv')
source=list(Sectors_edges['Company1'])
target=list(Sectors_edges['Company2'])
for i in range(len(Sectors_edges)):
    sectors_layer.add_edge(source[i],target[i])


#Continents layer:
continents_layer=nx.Graph()
continents_layer.add_nodes_from(nodes)
Continents_edges=pd.read_csv('Continents_network.csv')
source=list(Continents_edges['Company1'])
target=list(Continents_edges['Company2'])
for i in range(len(Continents_edges)):
    continents_layer.add_edge(source[i],target[i])

#Get infections data:
infected_companies_2020=estimate_infection.infected_companies_2020
days_2020=estimate_infection.days_2020

#Estimate parameters for different windows size+rolling window:
from scipy.optimize import minimize
pars=[0.01,0.01,0.01,0.01] #set initial values
bnds=((1e-10,0.99),(1e-10,0.99),(1e-10,0.99),(1e-10,0.99)) #boundary cond
for n in [1,2,3,4,5,6,7,8,9,10,20,30]: #window width
    w_10_1=list()
    for i in range(len(infected_companies_2020)-n-1):
        res = minimize(negative_log_likelihood, pars, args=(PMFG_layer, continents_layer, sectors_layer, global_layer,infected_companies_2020[i:i+n+1],days_2020[i:i+n+1]),bounds=bnds)
        w_10_1.append(res.x)
    #Save the results into a dataframe to use later:
    d=pd.DataFrame()
    b1=list()
    b2=list()
    b3=list()
    b4=list()
    for i in w_10_1:
        b1.append(i[0])
        b2.append(i[1])
        b3.append(i[2])
        b4.append(i[3])
    d['beta1']=b1
    d['beta2']=b2
    d['beta3']=b3
    d['beta4']=b4
    d.to_csv('2020_{}_cleared.xls'.format(n))

#Evaluate the recovery probability:
def negative_log_likelihood_recoveries(pars,recoveries,infected):
    p1=pars[0]
    l=1
    for i in range(len(recoveries)):
        l=l+recoveries[i]*np.log(p1)+(infected[i]-recoveries[i])*np.log(1-p1)
    return(-l)

recoveries_2020=list()
for i in range(len(infected_companies_2020)-1):
    recovered=len(set(infected_companies_2020[i])-set(infected_companies_2020[i+1]))
    recoveries_2020.append(recovered)

counted_2020=estimate_infection.counted_2020

#Recovries for window width n
for n in [1,2,3,4,5,6,7,8,9,10,20,30]:
    #here we read the previous dataframe to save a new column to it
    df=pd.read_csv('2020_{}_cleared.xls'.format(n))
    pars=[0.001] #initial value
    bnds=[(1e-15,0.99)] #boundary conditions
    recs=list()
    for i in range(len(infected_companies_2020)-n-1):
        res = minimize(negative_log_likelihood_recoveries, pars, args=(recoveries_2020[i:i+n+1],list(counted_2020.values())[i:i+n+1]),bounds=bnds,method='Nelder-Mead')
        res = res.x[0]
        recs.append(res)
    df['p']=recs
    df.to_csv('2020_{}_cleared.xls'.format(n))
    
########### Estimate parameters for CPD_model:
#Reestimate parameters after changepoints found by changefinder:

largest_scores=[25, 119, 176, 283, 356, 372, 416, 434, 559, 849] #2008
#largest_scores=[32, 52, 138, 161, 227, 245, 271, 369, 394, 417, 446] #2020

infected_companies_2008=estimate_infection.infected_companies_2008
days_2008=estimate_infection.days_2008

from scipy.optimize import minimize
pars=[0.01,0.01,0.01,0.01]
bnds=((1e-10,0.1),(1e-10,0.1),(1e-10,0.1),(1e-10,0.1))
for n in [1,2,3,4,5,6,7,8,9,10,20,30]:
    df=pd.read_csv('2020_{}_cleared.xls'.format(n))
    b1=list(df['beta1'])
    b2=list(df['beta2'])
    b3=list(df['beta3'])
    b4=list(df['beta4'])
    w_10_1=list()
    for i in range(len(b1)):
        res=[b1[i],b2[i],b3[i],b4[i]]
        for change_point in largest_scores:
            #if there is a change point in the previous n days:
            if change_point in range(i,i+n):
                res = minimize(negative_log_likelihood, pars, args=(PMFG_layer, continents_layer, sectors_layer, global_layer,infected_companies_2008[change_point:i+n+1],days_2008[change_point:i+n+1]),bounds=bnds)
                res=res.x
        w_10_1.append(res)
    d=pd.DataFrame()
    b1=list()
    b2=list()
    b3=list()
    b4=list()
    for i in w_10_1:
        b1.append(i[0])
        b2.append(i[1])
        b3.append(i[2])
        b4.append(i[3])
    d['beta1']=b1
    d['beta2']=b2
    d['beta3']=b3
    d['beta4']=b4
    d.to_csv('2020_{}_cleared_CPD.xls'.format(n))


#Collect the recoveries:
def log_likelihood_recoveries_2(pars,recoveries,infected):
    p1=pars[0]
    l=1
    for i in range(len(recoveries)):
        l=l+recoveries[i]*np.log(p1)+(infected[i]-recoveries[i])*np.log(1-p1)
    return(-l)

recoveries_2008=estimate_infection.recoveries_2008
counted_2008=estimate_infection.counted_2008

recoveries_2008=list()
for i in range(len(infected_companies_2008)-1):
    recovered=len(set(infected_companies_2008[i])-set(infected_companies_2008[i+1]))
    recoveries_2008.append(recovered)

#Recovries for window width n
for n in [1,2,3,4,5,6,7,8,9,10,20,30]:
    df=pd.read_csv('2020_{}_cleared_CPD.xls'.format(n))
    pars=[0.001]
    bnds=[(1e-15,0.99)]
    recs=list()
    df2=pd.read_csv('2020_{}_cleared.xls'.format(n))
    p=list(df2['p'])
    for i in range(len(df)):
        res=p[i]
        for change_point in largest_scores:
            #Reestimate if there is a change point in the previous n days:
            if change_point in range(i,i+n):
                res = minimize(log_likelihood_recoveries_2, pars, args=(recoveries_2008[change_point:i+n+1],(list(counted_2008.values()))[change_point:i+n+1]),bounds=bnds,method='Nelder-Mead')
                res=res.x[0]
        recs.append(res)
    df['p']=recs
    df.to_csv('2020_{}_cleared_CPD.xls'.format(n))
   