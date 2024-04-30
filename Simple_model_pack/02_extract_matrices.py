#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 12:17:49 2022

@author: diego
"""

import pandas as pd
import numpy as np
import os
import networkx as nx

import shelve as sl
import gzip
import json

import argparse

from analysis_functions import *

parser = argparse.ArgumentParser(description='Analyse simulations.')

# --- simulation parameters
parser.add_argument('--analysis', type=int, default=1,
                    help='which analysis to run')
parser.add_argument('--Npop', type=int,
                    help='size of the output matrix')

parser.add_argument('--folder', type=str, default='./SimResults/office/',
                    help='folder with results of simulations')
parser.add_argument('--folder_out', type=str, default='./SimResults/office/',
                    help='folder to save analysis')
parser.add_argument('--file', type=str,
                    help='file to analyse')
parser.add_argument('--label', type=str, 
                    help='''label that acts as decriptor for the analysis file. 
                            The output with have the format Cij_{key}.csv''')


args_read = parser.parse_args()

#%%   ------ functions ------------

#%%   --------------------- TO READ ONE set of conditions and obtain the infection pattern.
# this can be used to generate figure 3.a, and also needed for figure 4
folder = args_read.folder
folder_out  = args_read.folder_out

if args_read.analysis == 1:

    file_name = args_read.file
    label = args_read.label
    Npop = args_read.Npop
    
    # you can read all the files you want, but they are meant to have the same parameters
    list_of_files = [ folder + file_name ] 
    data,S = get_data( list_of_files  )
    
    cij = get_Cij(data,  Npop ) / data['run'].max() 
    np.savetxt(folder_out + f'Cij_{label}.csv' , cij, delimiter=";")



#%% -------------- TO make the plot of figure 4
elif args_read.analysis == 2:
    folder = './SimResults/office/'
    model0 = 'SIR_R0=4.0'
    # file with the runs done with the reference parameter, to recover the temporal detail
    file0 = 'InfPatt_dis=SIR_prot=NoTest_R=4.0_n=1000_imm=0.0_t=60.0_wden=0_ctest=auto.gzip'
    
    #older to save analysis
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
        
    # Saving the data to do the plot   
    
    with open('./fig2data.json', 'w') as json_file:
        json.dump(plots_data, json_file)


#%%
