import numpy as np
import networkx as nx
#from collections import deque
import random

import json
import csv
import os


class Threeshold_strength_resp_model:
    '''
    Describe...

    Attributes
    ----------
    I: dict
        details of infected people
    infected: list
        infected people


    Methods
    -------
    describe scheme
    '''


    def __init__(self, edgelist):
        '''
        Constructor.

        The method defines the setup for the simulations.

        Parameters
        ----------
        edgelist: list
            interacting nodes
        threshold: float
            threshold parameter (the same for each node)
        '''

        self.nodes_list = list(np.unique(list(edgelist.keys())))
        self.contacts = dict()
        self.edgelist = edgelist # è un dizionario
        self.C = dict()

        self.store_contacts()

        #self.initialize_time0()






    def store_contacts(self):
        '''
        Store list of contacts for each node from edgelist.

        '''

        for w_edge in self.edgelist:

            if w_edge[0] in self.contacts:
                self.contacts[w_edge[0]].append((w_edge[1],self.edgelist[w_edge])) # append(neighbour, strength)
            else:
                self.contacts[w_edge[0]] = [(w_edge[1],self.edgelist[w_edge])]
            if w_edge[1] in self.contacts:
                self.contacts[w_edge[1]].append((w_edge[0],self.edgelist[w_edge]))
            else:
                self.contacts[w_edge[1]] = [(w_edge[0],self.edgelist[w_edge])]








    def initialize_time0(self, pat0):
        '''
        Initialize the status of the initial infected people and set initial edges weigths at 0.

        '''

        self.pat0 = pat0 # initial infected
        self.Istate = np.full(len(self.nodes_list),0)
        self.Istate[self.pat0] = 1
        self.Rstate = np.full(len(self.nodes_list),0)
        self.nb_Ineighs = []
        for edge in self.edgelist:
            self.C[tuple(edge)] = 0
            self.C[tuple(edge[::-1])] = 0



    def simulate(self,threshold,mu,t_i):
        '''
        Run the simulation.

        Returns
        ----------
        describe...
        '''
        tmax = 300
        mean_deg_t = np.full(tmax,0)
        flag_inf = False


        for t in range(tmax):

            Sstate = 1 - self.Istate - self.Rstate
            if any(c == -1 for c in Sstate):
                print('error')
                break
            old_susc = list(np.where(Sstate==1)[0])
            old_inf = list(np.where(self.Istate==1)[0])

            for i in old_susc:
                strength = 0
                Istrength = 0
                Ineighs = []
                for neigh in self.contacts[i]:
                    strength += neigh[1]
                    if neigh[0] in old_inf:
                        Istrength += neigh[1] # strength of i-neigh[0]
                        Ineighs.append(neigh[0])
                        #print('neigh',i,Ineighs,Istrength)
                        #print('cond',Istrength, threshold*strength,threshold,strength)
                if Istrength > threshold*strength:
#                    self.infected.append(i)
#                    self.susceptible.remove(i)
                    self.Istate[i] = 1
                    flag_inf = True
                    t_i[i].append(t+1)
                    self.nb_Ineighs.append(len(Ineighs))

                    for j in Ineighs:
                        self.C[tuple((j,i))] = self.edgelist[tuple(np.sort((i,j)))]/Istrength
                        #print('weight change',self.C[tuple(np.sort((i,j)))])

            if mu != 0:
                for i in old_inf:
                    r = np.random.uniform(0, 1)
                    if r < mu:
                        self.Istate[i] = 0
                        self.Rstate[i] = 1


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
            if list(np.where(Sstate==1)[0]) == old_susc:
                #print('stacked at t %d'%t) # stallo
                break

        if Sstate.sum() != 0 and list(np.where(Sstate==1)[0]) != old_susc:
            print('not finished: ci sono ancora suscettibili')
        #elif self.Istate.sum() != 0 and list(np.where(Sstate==1)[0]) != old_susc:
        #    print('not finished: ci sono ancora infetti')

        return self.C, flag_inf, t_i, self.nb_Ineighs
