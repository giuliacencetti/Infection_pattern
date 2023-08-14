import numpy as np
import networkx as nx
import random

import json
import csv
import os
#import time

#random.seed(1)

class Simplagion_model:


    def __init__(self, edgelist, edgelist3, nodes_list, weighted, final_result):
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
        edgelist3: dict
            3-body interactions
            keys: triplets, values: strength
        nodes_list: list
            nodes
        weighted: bool
            True if weighted network
        final_results: str
            'inf_tree', 'superspreader','sing_traj', 'sing_traj_link_tri'
        '''

        self.nodes_list = nodes_list
        self.contacts =  [ [] for _ in range(len(self.nodes_list)) ] # pairwise contacts
        self.contactsT = [ [] for _ in range(len(self.nodes_list)) ] # triangular contacts
        self.edgelist = edgelist
        self.edgelist3 = edgelist3
        self.weighted = weighted
        if self.weighted:
            self.weights =  [ [] for _ in range(len(self.nodes_list)) ] # pairwise weights
            self.weightsT =  [ [] for _ in range(len(self.nodes_list)) ] # triangular weights

        self.final_result = final_result

        self.store_contacts()

        if final_result == 'inf_tree':
            self.init_contagions()


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
            if self.weighted:
                self.weights[edge[0]].append(self.edgelist[edge])
                self.weights[edge[1]].append(self.edgelist[edge])

        '''
        if self.final_result == 'sing_traj':
            if self.weighted:
                self.degrees = [sum(w) for w in self.weights]
            else:
                self.degrees = [len(c) for c in self.contacts]
        '''

        # TRIANGULAR CONTACTS
        for triangle in self.edgelist3:
            self.contactsT[triangle[0]].append((triangle[1],triangle[2]))
            self.contactsT[triangle[1]].append((triangle[2],triangle[0]))
            self.contactsT[triangle[2]].append((triangle[0],triangle[1]))
            if self.weighted:
                self.weightsT[triangle[0]].append(self.edgelist3[triangle])
                self.weightsT[triangle[1]].append(self.edgelist3[triangle])
                self.weightsT[triangle[2]].append(self.edgelist3[triangle])






    def init_contagions(self):

        self.empty_C_L = dict()
        for edge in self.edgelist:
            self.empty_C_L[tuple(edge)] = 0
            self.empty_C_L[tuple(edge[::-1])] = 0

        self.empty_C_T = self.empty_C_L.copy() # such that couples appear in the same order








    def initialize_time0(self,seed_size):
        '''
        Initialize the status of the initial infected people and set initial edges weigths at 0.

        '''
        # Initialize patient 0:



        self.Istate = np.full(len(self.nodes_list),0)
        self.Rstate = np.full(len(self.nodes_list),0)

        if self.final_result == 'inf_tree':
            for n in random.sample(list(self.nodes_list), seed_size):
                self.Istate[n] = 1
            self.C_L = self.empty_C_L.copy()
            self.C_T = self.empty_C_T.copy()



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
        #flag_inf = False
        '''
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
        '''
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
            old_inf = list(np.where(self.Istate==1)[0])

            for i in old_susc:

                set_old_inf = set(old_inf)
                if self.weighted:
                    # pairwise infected contacts
                    Ineighs_ix = np.array([(ind,x) for (ind, x) in enumerate(self.contacts[i]) if x in set_old_inf],dtype=object)
                    if Ineighs_ix.size == 0:
                        indexes = []
                        Ineighs = []
                    else:
                        indexes = list(Ineighs_ix[:,0])
                        Ineighs = list(Ineighs_ix[:,1])
                    Ineighs_strength = np.array([self.weights[i][n] for n in indexes])
                    # triangular infected contacts
                    IneighsT_ix = np.array([(ind,x) for (ind, x) in enumerate(self.contactsT[i]) if x[0] in set_old_inf and x[1] in set_old_inf],dtype=object)
                    #IneighsT_ix = np.array([(ind,x) for (ind, x) in enumerate(self.contactsT[i]) if self.Istate[x[0]]==1 and self.Istate[x[1]]==1])
                    if IneighsT_ix.size == 0:
                        IneighsT = []
                        IneighsT_strength = np.array([])
                        p = 1 - np.prod(1-beta*Ineighs_strength)
                    else:
                        indexesT = list(IneighsT_ix[:,0])
                        IneighsT = list(IneighsT_ix[:,1])
                        IneighsT_strength = np.array([self.weightsT[i][n] for n in indexesT])
                        p = 1 - np.prod(1-beta*Ineighs_strength)*np.prod(1-betaT*IneighsT_strength) # infection probability
                    if p > 1 or p<0:
                        print('problema! p =',p)

                else:
                    Ineighs = [neigh for neigh in self.contacts[i] if neigh in set_old_inf]
                    IneighsT = [couple for couple in self.contactsT[i] if couple[0] in set_old_inf and couple[1] in set_old_inf]

                    p = 1 - (1-beta)**len(Ineighs)*(1-betaT)**len(IneighsT) # infection probability
                    if p < 0 or p > 1:
                        print('error: p = ',p)
                r = np.random.uniform(0, 1)
                if r < p:
                    self.Istate[i] = 1
                    #flag_inf = True
                    t_i[i].append(t+1)

                    if self.final_result == 'inf_tree':
                        if self.weighted:
                            norm = beta*sum(Ineighs_strength) + betaT*sum(IneighsT_strength)
                            for n in range(len(Ineighs)):
                                self.C_L[tuple((Ineighs[n],i))] += beta*Ineighs_strength[n]/norm
                            for n in range(len(IneighsT)):
                                self.C_T[tuple((IneighsT[n][0],i))] += betaT*IneighsT_strength[n]/2/norm
                                self.C_T[tuple((IneighsT[n][1],i))] += betaT*IneighsT_strength[n]/2/norm

                        else:
                            for j in Ineighs:
                                #self.C_L[j,i] += beta/(len(Ineighs)*beta + len(IneighsT)*betaT)
                                self.C_L[tuple((j,i))] += beta/(len(Ineighs)*beta + len(IneighsT)*betaT)
                            for couple in IneighsT:
                                #self.C_T[couple[0],i] += betaT/(len(Ineighs)*beta + len(IneighsT)*betaT)/2
                                #self.C_T[couple[1],i] += betaT/(len(Ineighs)*beta + len(IneighsT)*betaT)/2
                                self.C_T[tuple((couple[0],i))] += betaT/(len(Ineighs)*beta + len(IneighsT)*betaT)/2
                                self.C_T[tuple((couple[1],i))] += betaT/(len(Ineighs)*beta + len(IneighsT)*betaT)/2



                    #print("3--- %.2e seconds ---" % (time.time() - start_time))
                #print('-------------------------------')
            if mu != 0:

                for i in old_inf:
                    r = np.random.uniform(0, 1)
                    if r < mu:
                        self.Istate[i] = 0
                        self.Rstate[i] = 1



            '''
            if save_temp_ev:
                Inf.append(self.Istate.sum())
                Rec.append(self.Rstate.sum())
            '''
            '''
            if mu > 0 and self.Istate.sum() == 0: # everybody has recovered
                #print('all recovered')
                break

            if mu == 0 and self.Istate.sum() == len(self.nodes_list):
                #print('all infected at t=%d'%t)
                break
            '''


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
        #print('2',t_i)

        if Sstate.sum() != 0 and self.Istate.sum() != 0 and nn != NN:
            print('not finished: ci sono ancora suscettibili e infetti')
        #elif self.Istate.sum() != 0 and list(np.where(Sstate==1)[0]) != old_susc:
        #    print('not finished: ci sono ancora infetti')


        #if t == tmax-1:
            #print('t_max infetti=%d nn=%d'%(self.Istate.sum(),nn))


        '''
        if save_temp_ev:
            if mu > 0:
                filenameI = 'results/temp_evolution/SIR_temp_ev_Inf_hospital_beta_%.4f_betaT_%.4f_mu_%.2f.csv'%(beta,betaT,mu)
                filenameR = 'results/temp_evolution/SIR_temp_ev_Rec_hospital_beta_%.4f_betaT_%.4f_mu_%.2f.csv'%(beta,betaT,mu)
            elif mu==0:
                filenameI = 'results/temp_evolution/SI_temp_ev_Inf_hospital_beta_%.4f_betaT_%.4f_mu_%.2f.csv'%(beta,betaT,mu)
                filenameR = 'results/temp_evolution/SI_temp_ev_Rec_hospital_beta_%.4f_betaT_%.4f_mu_%.2f.csv'%(beta,betaT,mu)
            save_on_csv(filenameI,Inf,'a')
            save_on_csv(filenameR,Rec,'a')
        '''


        C_L_list = np.array(list(self.C_L.values()))
        C_T_list = np.array(list(self.C_T.values()))
        a = self.Istate.sum() + self.Rstate.sum()

        return C_L_list, C_T_list, t_i, a
        #return C_L_list, C_T_list, t_i, a, t
