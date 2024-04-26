#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 12:17:49 2022

@author: diego
"""

import pandas as pd
import numpy as np
import pylab as plt
import seaborn as sns
import os
import networkx as nx

import shelve as sl
import gzip
import json

from scipy.optimize import curve_fit
from scipy.stats import kendalltau



#%%   ------ functions ------------

def get_gen(G,u,g):
    gen = {}
    if len(list(G.neighbors(u)))>0:
        for v in G.neighbors(u):
            gen[v] = g +1
            gen.update( get_gen(G,v,g+1) ) 
    return gen
def get_dictgen(r):

    u0 = r[ 'external_infections'][0]    
    
    Gedges = r['Gedges']
    Gnodes = r['Gnodes']
    
    G =nx.DiGraph()
    #G.add_nodes_from( [ (int(u),{'time': t})  for u, t in Gnodes if t==t]   )
    G.add_edges_from( [ (int(u),int(v),{'time': t})  for u, v, t in Gedges]   )
    
    dgen = get_gen(G,u0,0)
    return dgen,dict(G.out_degree())


def get_Dfg(filename,index0=0 ):
    with gzip.open(filename, 'r') as infile:
        jread = infile.read()
        jread= jread.decode('utf-8')
    js = json.loads( jread  )
    Df = pd.DataFrame()
    Rs = js['runs']
    Sizes =[]
    #print(js['Param']['beta'])
    for ind, r in enumerate(Rs):
        x  = np.array( r['Gedges'] )
        if x.size:
            df = pd.DataFrame(x, columns = ['u','v','time',])
            df['run'] = ind + index0 + 1
            
            dgen,dR = get_dictgen(r)
            df['generation'] = df['u'].map(dgen)
            df['R'] = df['u'].map(dR)
            df['size'] = r['final_size']

            Df = pd.concat([Df,df])
        #Sizes.append((ind + index0 + 1, r['final_size'], r['R0'], get_timesind( np.array( r['size(t)'] ) )  )   )
        Sizes.append((ind + index0 + 1, r['final_size'] )   )
    return Df, Sizes, ind + index0 + 1  


def get_data(files):
    data = None
    index0 = 0
    Df = pd.DataFrame()
    Sizes = []
    R0 = []

    for file in files:
        df, s,  index0 = get_Dfg(file,index0=index0)
        Df = pd.concat([Df,df])
        Sizes += s
    Df = Df.reset_index()
    return Df,Sizes

def get_Cij(df,Npop):
    A = np.zeros((Npop,Npop))
    for ind in df.index:
        u,v = df.loc[ind][['u','v']]
        if v==v:
            A[int(v),int(u)] += 1
    return A      

def get_Cijt(df0,times,ntrial,Npop):
    dCij_t = []
    for t0,tf in zip(times[:-1],times[1:]):
        #print(t0,tf)
        df =df0[ (df0['time']>=t0) & (df0['time']<tf)]
        A = get_Cij( df , Npop)/ntrial
        dCij_t.append(A)
    return dCij_t 

def get_accCij(dCij_t):
    Cij_t = []
    Asum = np.zeros((Npop,Npop))
    for A in dCij_t:
        Asum  += A
        Cij_t.append(Asum.copy() )
    return Cij_t  

def get_hist(c):
    bins = np.linspace(2, 250, 100)
    hist0 = np.histogram(c, bins=bins)[0]
    hist0 = hist0**2
    return hist0/hist0.sum()

def fun_size(c):
    return np.mean(c[c>1])

def to_list(X):
    return [[ float(y) for y in x] for x in X]
#%%   --------------------- TO READ ONE set of conditions and obtain the infection pattern.
# this can be used to generate figure 3.a, and also needed for figure 4

folder = './SimResults/test/'
folder_out  = './SimResults/test/'
file0 = 'test_dis=SIR_prot=NoTest_R=3.5_n=100_imm=0.0_t=60.0_wden=0_ctest=auto.gzip'
key = 'SIR_R0=3.5'

# you can read all the files you want, but they are meant to have the same parameters
list_of_files = [ folder + file0 ] 
data,S = get_data( list_of_files  )

# introduce the population of the dataset
Npop = 75
cij = get_Cij(data,  Npop ) / data['run'].max()
#
np.savetxt(folder_out + f'Cij_{key}.csv' , cij, delimiter=";")

#df, s,  index0 = get_Dfg(folder + file0,index0=0)


#%% -------------- TO make the plot of figure
model0 = 'SIR_R0=4.0'
model1 = 'SIR_R0={}'
file0 = 'test_dis=SIR_prot=NoTest_R=3.5_n=100_imm=0.0_t=60.0_wden=0_ctest=auto.gzip'



times = np.arange(1,15,0.1)

dCij_t  = get_Cijt(data, times, data['run'].max(), Npop )
Cij_t = get_accCij(dCij_t) 


size_ref_t = np.array( [  np.array(
          data.query(f'time < {t}').groupby('run').count()['size'] 
          )+1 for t in times ] )


hist0_t = np.array([ get_hist(c) for c in size_ref_t ])
size0 = np.array([ fun_size(c) for c in size_ref_t ])[:-1]/Npop



plots_data = {'time': list(times) , 
              'reference_' + model0: 
                 {'size(t)': list(size0), 
                  'size_dist(t)': to_list(hist0_t),
                  'dist_size(t)': to_list(size_ref_t ),
                  'Cij(t)': [to_list(c) for c in Cij_t]}}


for indx, (R0,color) in enumerate( [(1.5,'C0'),(2.0,'C1'),(2.5,'C2'),(3.0,'C3'),(3.5,'C4')]):
    key = model1.format(R0)
    label = '$R_0^1={}$'.format(R0)

    file_path = folder_out + f'Cij_{key}.csv'
    Cij1 = np.loadtxt(file_path, delimiter=";")
    y= np.array([ cossim( c,Cij1 ) for c in Cij_t ])
    

    #itmin = np.absolute(size0 - size1).argmin()
    itmin = y[1:].argmax() +1    
    
    dis_size1 = np.array([v[1] for v in dict_size[key] ] )
    size1 = fun_size(dis_size1)/Npop
    
    plots_data[key]={ 'Cij': to_list(Cij1),'sim': list(y), 'size_tmin': float(size1),
                     'tmin': float(time[itmin]), 'itmin': int(itmin), 
                     'dist_size': [int(s) for s in dis_size1] ,
                    }
   
import json
with open('./fig2data.json', 'w') as json_file:
    json.dump(plots_data, json_file)


#%%
