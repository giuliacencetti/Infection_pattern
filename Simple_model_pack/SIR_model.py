#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 06:48:21 2022

@author: diego
"""
import numpy as np
import EpiKit as epi
import dictfit
S         =            ['S',  'I', 'R']
dict_of_rel_inf   = {
    'adult': np.array([0,  1.00, 0, ]),
    'children': np.array([0,  1.00, 0, ]),
    }


def beta_from_R0(R0, data, net, r_period, weekdays, generation=1):
    if  data == 'office' and weekdays==7:
        DictFit = dictfit.get_dict['office']
    elif data == 'office' and r_period!=1:
        Dkey = {0.25:'office_r025',0.5:'office_r05',1.5:'office_r15',}
        DictFit = dictfit.get_dict[Dkey[args_read.r_period]]
    else:
        DictFit = dictfit.get_dict[data]
    
    if net in DictFit[generation].keys():
        A,f = DictFit[generation][net]
    else:
        print('---> net {} not in the list'.format(args_read.net))
        A,f = DictFit[1]['DYN']
    return -np.log( 1 - R0/A)/ (f  /(60*24) )

def get_disease_model(beta,ParamSim, data):
    # if ParamSim['tauI']==0:
    #     dict_of_rel_inf   = {
    #         'adult': np.array([0,  1.00, 1.0, ]),
    #         'children': np.array([0,  1.00, 1.0, ]),
    #         }
    
    
    if data in ['office','hospital','conf']:
        relative_infectiousness = np.array(
              [dict_of_rel_inf['adult'],
               dict_of_rel_inf['adult'] ])
    elif data in ['cprepa', 'school']:
        relative_infectiousness = np.array(
              [dict_of_rel_inf['adult'],
               dict_of_rel_inf['adult'] ])   
    S_inf = relative_infectiousness * beta 
    r_period = ParamSim['r_period'] 
    Param = { 'tauI' :  ParamSim['tauI'],
              'pgamma': ( 1/ParamSim['pvar'] )**2,
              'sigma' : np.array([1.0, 1.0  ]),
              'weekden': ParamSim['weekden'] ,
              
              'beta': beta,
              'r_period':r_period ,
              }
    #Param.update(ParamVax)
    if Param['tauI']>0:
        React = [ epi.Reaction('S',  'I', beta ), #0
                  epi.Reaction('I',  'R',  1/Param['tauI'] , gamma=Param['pgamma'] ) ] #2
    else:
        React = [ epi.Reaction('S',  'I', beta ), #0
                  epi.Reaction('I',  'R', 0  ) ] #2


    Dict_states = {'healthy':'S',
                   'after_infection':'I',
                   'recovered':'R'}


    return Param, React, S_inf,Dict_states


    # S         =            ['S',  'E',  'Ip',  'Ic',  'Isc', 'Rp', 'R']
    # S_inf0    = np.array([ [ 0,    0,0.55*0.63, 0.63,0.55*0.63,  0,   0, ],
    #                        [ 0,    0,   0.55,   1.00,   0.55,    0,   0, ],
    #                        ])
    # SIGMA_STUD = 0.5
    # PROD_DET = 0.3
# Param_out, React, Sinf = get_covid_model(beta,Param)



