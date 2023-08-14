import numpy as np
import networkx as nx
import random

import json
import csv
import os
#import time

#random.seed(1)

class Simplagion_model:


    def __init__(self, edgelist, nodes_list, final_result):
        '''
        Constructor.

        The method defines the setup for the simulations.

        Parameters
        ----------
        edgelist: dict
            all interactions
            keys: edges, values: strength
        edgelist2: dict
            purely pairwise interactions
            keys: edges, values: strength
        nodes_list: list
            nodes
        final_results: str
            'inf_tree', 'superspreader','sing_traj', 'sing_traj_link_tri'
        '''

        self.nodes_list = nodes_list
        self.contacts =  [ [] for _ in range(len(self.nodes_list)) ] # pairwise contacts
        self.edgelist = edgelist
        self.weights =  [ [] for _ in range(len(self.nodes_list)) ] # pairwise weights

        self.final_result = final_result

        self.store_contacts()




        #self.initialize_time0()

        #self.simulate()





    def store_contacts(self):
        '''
        Store list of contacts for each node from edgelist.

        '''
        # PAIRWISE CONTACTS
        for edge in self.edgelist:
            self.contacts[edge[0]].append(edge[1])
            self.contacts[edge[1]].append(edge[0])
            self.weights[edge[0]].append(self.edgelist[edge])
            self.weights[edge[1]].append(self.edgelist[edge])







    def init_contagions(self):

        self.empty_C_L = dict()
        for edge in self.edgelist:
            self.empty_C_L[tuple(edge)] = 0
            self.empty_C_L[tuple(edge[::-1])] = 0








    def initialize_time0(self,seed_size):
        '''
        Initialize the status of the initial infected people and set initial edges weigths at 0.

        '''
        # Initialize patient 0:

        if self.final_result == 'inf_tree':
            self.init_contagions()

        elif self.final_result == 'R0':
            self.R0 = 0

        self.Istate = np.full(len(self.nodes_list),0)
        self.Rstate = np.full(len(self.nodes_list),0)

        if self.final_result == 'inf_tree':
            for n in random.sample(list(self.nodes_list), seed_size):
                self.Istate[n] = 1
            self.C_L = self.empty_C_L.copy()

        elif self.final_result == 'R0':
            self.pat0 = random.randint(0,len(self.nodes_list)-1)
            self.Istate[self.pat0] = 1



    def simulate(self, beta, betaT, mu, save_temp_ev, t_i):
        '''
        Run the simulation.

        Returns
        ----------
        describe...
        '''
        #filename = 'fileprova.txt'
        #file = open(filename, 'a+')
#        inf_thresh = 2
#        flag_nb_inf = False
        flag_inf = False
        if save_temp_ev:
            Inf = [self.Istate.sum()]
            Rec = [self.Rstate.sum()]

            def save_on_csv(filename, variable_list, writing_operation):
                with open(filename, writing_operation) as csvfile:
                    writer = csv.writer(csvfile)
                    try:
                        [writer.writerow(s) for s in variable_list]
                    except:
                        writer.writerow(variable_list)

        #print('elementi div da 0',len([el for el in list(self.C_L.values()) if el!=0]))

        nb_pat0r = self.Istate.sum()
        tmax = 6000
        nn = 0
        for t in range(tmax):

            Sstate = 1 - self.Istate - self.Rstate
            if any(c == -1 for c in Sstate):
                print('error')
                break
            old_susc = list(np.where(Sstate==1)[0])
            set_old_susc = set(old_susc)
            old_inf = list(np.where(self.Istate==1)[0])

            for i in old_inf:
                Sneighs_ix = np.array([(ind,x) for (ind, x) in enumerate(self.contacts[i]) if x in set_old_susc],dtype=object)
                if Sneighs_ix.size == 0:
                    indexes = []
                    Sneighs = []
                else:
                    indexes = list(Sneighs_ix[:,0])
                    Sneighs = list(Sneighs_ix[:,1])
                for n in range(len(indexes)):
                    j = Sneighs[n]
                    w = self.weights[i][indexes[n]]
                    r = np.random.uniform(0, 1)
                    if beta*w < 0 or beta*w > 1:
                        print('errore')
                    if r < beta*w:
                        self.Istate[j] = 1
                        if self.final_result == 'inf_tree':
                            norm = beta*sum(Ineighs_strength)
                            for n in range(len(Ineighs)):
                                self.C_L[tuple((Ineighs[n],i))] += beta*Ineighs_strength[n]/norm
                        elif self.final_result == 'R0':
                            if i == self.pat0:
                                self.R0 += 1

                    #print("3--- %.2e seconds ---" % (time.time() - start_time))
                #print('-------------------------------')
            if mu != 0:

                for i in old_inf:
                    r = np.random.uniform(0, 1)
                    if r < mu:
                        self.Istate[i] = 0
                        self.Rstate[i] = 1




            if save_temp_ev:
                Inf.append(self.Istate.sum())
                Rec.append(self.Rstate.sum())

#            NotSusc_nb = self.Istate.sum() + self.Rstate.sum()
#            if mu > 0 and flag_nb_inf == False and NotSusc_nb >= inf_thresh:
#                flag_nb_inf = True # ho superato la soglia: qualcuno si è infettato

            # Motivi per interrompere:
            # 1- Non ci sono più suscettibili
            # 2- Non ci sono più infetti
            # 3- Stallo, la soglia che ho scelto non mi permette di andare avanti

            # 1
            Sstate = 1 - self.Istate - self.Rstate
            if Sstate.sum() == 0: # non ci sono più suscettibili
                #print('No more susc')
                break
            # 2
            if self.Istate.sum() == 0: # non ci sono più infetti
                #print('No more inf')
                break
            # 3
            NN = 300
            if list(np.where(Sstate==1)[0]) == old_susc:
                nn += 1
            else:
                nn = 0
            if nn == NN:
                #print('stuck at t=%d with %d susceptible'%(t,Sstate.sum()))
                break


        if Sstate.sum() != 0 and self.Istate.sum() != 0 and nn != NN:
            print('not finished: ci sono ancora suscettibili e infetti')
        #elif self.Istate.sum() != 0 and list(np.where(Sstate==1)[0]) != old_susc:
        #    print('not finished: ci sono ancora infetti')


        #if t == tmax-1:
            #print('t_max infetti=%d nn=%d'%(self.Istate.sum(),nn))




        if self.final_result == 'inf_tree':
            C_L_list = np.array(list(self.C_L.values()))
            return C_L_list, flag_inf, t_i

        elif self.final_result == 'R0':
            a = self.Istate.sum() + self.Rstate.sum()
            return self.R0,a
