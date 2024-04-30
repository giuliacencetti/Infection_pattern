#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import os
import networkx as nx

import shelve as sl
import gzip
import json


from analysis_functions import *



#%% -------------- TO make the plot of figure 4

# reference model
folder = './SimResults/office/' # folder with the output of the simulation
model0 = 'SIR_R0=4.0'
# file with the runs done with the reference parameter, to recover the temporal detail
file0 = 'InfPatt_dis=SIR_prot=NoTest_R=4.0_n=1000_imm=0.0_t=60.0_wden=0_ctest=auto.gzip'

#models to compare
folder_matrices = './SimResults/office/' # folder with the matrices to use to compare
model1 = 'SIR_R0={}'

times = np.arange(1,15,0.1)

data,S = get_data( [ folder + file0 ]  )
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

    file_path = folder_matrices + f'Cij_{key}.csv'
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
    
# Saving the data to do the plot   
with open('./fig2data.json', 'w') as json_file:
    json.dump(plots_data, json_file)


