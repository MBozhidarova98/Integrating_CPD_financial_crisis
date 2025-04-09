# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 14:38:16 2023

@author: mbozh
"""

import pandas as pd
import estimate_infection_github_code
import itertools
import networkx as nx
import planarity

#Get the upper tail dependence coefficients:
df_tail_correlation=pd.read_csv('chis_95_before_2020.csv')
x=list(df_tail_correlation['chis'])
weights=list()
for i in x:
    weights.append(1-i**2)
df_tail_correlation['Weights']=weights

l=estimate_infection_github_code.log_returns
nodes=list(l.columns)[:-1]
G=nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(itertools.permutations(nodes, 2))

for edge in G.edges():
    if ((df_tail_correlation['Company1']== edge[0]) & (df_tail_correlation['Company2']== edge[1])).any():
        G[edge[0]][edge[1]]['weight']=df_tail_correlation[(df_tail_correlation['Company1']== edge[0]) & (df_tail_correlation['Company2']== edge[1])]['Weights'].values[0]
    elif ((df_tail_correlation['Company1']== edge[1]) & (df_tail_correlation['Company2']== edge[0])).any():
        G[edge[0]][edge[1]]['weight']=df_tail_correlation[(df_tail_correlation['Company1']== edge[1]) & (df_tail_correlation['Company2']== edge[0])]['Weights'].values[0]
    else:
        G[edge[0]][edge[1]]['weight']=1

def sort_graph_edges(G):
    sorted_edges = []
    for source, dest, data in sorted(G.edges(data=True),
                                     key=lambda x: x[2]['weight']):
        sorted_edges.append({'source': source,
                             'dest': dest,
                             'weight': data['weight']})
        
    return sorted_edges

def compute_PMFG(sorted_edges, nb_nodes):
    PMFG = nx.Graph()
    for edge in sorted_edges:
        PMFG.add_edge(edge['source'], edge['dest'])
        if not planarity.is_planar(PMFG):
            PMFG.remove_edge(edge['source'], edge['dest'])
            
        if len(PMFG.edges()) == 3*(nb_nodes-2):
            break
    
    return PMFG   
    
#MST=nx.minimum_spanning_tree(G)
#nx.local_efficiency(MST) #0.0
#nx.global_efficiency(MST) #0.1397096650927733

sorted_edges = sort_graph_edges(G)
PMFG_layer = compute_PMFG(sorted_edges, len(G.nodes))
#nx.draw(PMFG_layer,node_size=10)

#Put the edges to be the tail dependence coefficients:
for edge in PMFG_layer.edges():
    PMFG_layer[edge[0]][edge[1]]['weight']=0
    if ((df_tail_correlation['Company1']== edge[0]) & (df_tail_correlation['Company2']== edge[1])).any():
        PMFG_layer[edge[0]][edge[1]]['weight']=abs(df_tail_correlation[(df_tail_correlation['Company1']== edge[0]) & (df_tail_correlation['Company2']== edge[1])]['chis'].values[0])
    elif ((df_tail_correlation['Company1']== edge[1]) & (df_tail_correlation['Company2']== edge[0])).any():
        PMFG_layer[edge[0]][edge[1]]['weight']=abs(df_tail_correlation[(df_tail_correlation['Company1']== edge[1]) & (df_tail_correlation['Company2']== edge[0])]['chis'].values[0])

edges_1=list()
edges_2=list()
weights=list()
for edge in PMFG_layer.edges():
    edges_1.append(edge[0])
    edges_2.append(edge[1])
    weights.append(PMFG_layer[edge[0]][edge[1]]['weight'])

df_PMFG=pd.DataFrame()
df_PMFG['Source']=edges_1
df_PMFG['End']=edges_2
df_PMFG['weight']=weights    

#Save the PMFG into a dataframe:
df_PMFG.to_csv('PMFG_chis_95_before_2020.csv')