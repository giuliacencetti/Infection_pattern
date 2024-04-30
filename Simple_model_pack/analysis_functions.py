#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
import networkx as nx

import shelve as sl
import gzip
import json


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
    # this function should be modified for a dataset different that the office
    bins = np.linspace(2, 250, 100)
    hist0 = np.histogram(c, bins=bins)[0]
    hist0 = hist0**2
    return hist0/hist0.sum()

def fun_size(c):
    return np.mean(c[c>1])

def to_list(X):
    return [[ float(y) for y in x] for x in X]