# -*- coding: utf-8 -*-

"""
Created on June 14, 2016

Written in Python 3


"""


import numpy as np
from numpy import random
import itertools

class NK_model(object):
    
    def __init__(self, n, k):
        self.n = n
        self.k = k
        self.genotype_list = [''.join(i) for i in itertools.product('01', repeat = self.n)] # generate all the possible 
                                                                   #genotypes with length n (two possible alleles in each site)
        self.start_type = 0  
        self.get_allelic_fitness() # get the individual fitness of each position
        self.get_epistic_site()  # get the epistic network among the sites
        self.cal_fitness()  # calculate the fitness of every genotype
        self.landscape = dict(zip(self.genotype_list, self.w_fitness)) # generate the landscape
        
        
    def get_allelic_fitness(self):  
        '''Generate the individual fitness of each alleles'''
        
        self.fitness_list = []
        
        epistic_geno = [''.join(i) for i in itertools.product('01', repeat = self.k+1)]  # get all possible combinations of 
                                # the site i with its epistic sites (because including site i itself, repeat = self.k+1)
        for j in range(self.n): # for each site in the genotype, generate the corresponding fitness
            epistic_fit = random.random(len(epistic_geno)) 
            epistic_map = dict(zip(epistic_geno, epistic_fit)) # generate the fitness of each site-epistic site pairs
            self.fitness_list.append(epistic_map)  # generate the fitness list
            
    
    def get_epistic_site(self):
        '''Generate the epistic networks among the sites'''

        self.epistic_site = []
        for i in range(self.n):
            all_site = [x for x in range(self.n) if x != i]  # # get the indexes of all sites excluding i
            
            indices = set()   # will be used to store the indexes of epistic sites
            while len(indices) < self.k:
                indices.add(random.randint(0, len(all_site))) # randomly pick the indexes of k elements from all_site

            indices = list(indices)  # convert the set into list
            indices.sort()
            
            one_site = []
            for j in indices:  # get the epistic sites of site i
                one_site.append(all_site[j])
            
            self.epistic_site.append(one_site)
            
    
    
    def cal_fitness(self):
        '''Calculate the fitness of every genotype'''
        
        self.w_fitness = []  # will be used to store the fitness of each genotype in the landscape
        
        for i in self.genotype_list:
            one_fitness = []  # will be used to store the fitness of each site in the genotype i
            for j in range(self.n):
                epis_pair = i[j]  # the calculated site itself
                for s in self.epistic_site[j]: # the epistic site
                    epis_pair +=i[s]  # get all epistic site into epis_pair (will be used to get the fitness of this pair)
                    
                one_fitness.append(self.fitness_list[j][epis_pair])  # get the fitness of the epis_pair (fitness of site j in
                                                                    # genotype i)
            
            i_fitness = sum(one_fitness)/self.n  # calcualte the fitness of genotype i and store into self.w_fitness
            self.w_fitness.append(i_fitness)
            
    
    @staticmethod
    def creat_neighbor(genotype, mutant_num):
        '''creat the mutant_num-mutant neighbors of the given genotype
        
        parameters: 
        ----------
            genotype: str that repesents genotype, eg: "01", "00"
            mutant_num: int, the number of mutant site to creat neighbors
            
        output: 
        -------
            list that store the mutant_num-mutant neighbors of the given genotype
        '''
        
        total_neigh_list= []  # will be used to store the created neighbors of the given genotype
    
        geno_list=list(genotype) # convert the given genotype (str) into a list
    
        indices = list(range(len(geno_list))) # get all the indices for the given genotype
    
        mutant_indices = list(itertools.combinations(indices, mutant_num)) # get all possible combinations of picking mutant_num
                                                                    # intergers from indices (i.e. mutant sites)
    
        for i in mutant_indices:  # get all its mutant_num -mutant neighbors of the given genotype
            neigh_list = geno_list.copy()
            for j in i:
                if geno_list[j] == '0':
                    neigh_list[j] = '1'
                else:
                    neigh_list[j] = '0'
            
            total_neigh_list.append(''.join(neigh_list))
    
        return total_neigh_list
    
    
    def random_walk(self):
        '''Simulate one step  of random walk along the landscape'''
        
        if self.start_type==0:  # self.start_type==0 (the initial given value), means initial step, randomly selects a genotype
                                # as start_type
            self.start_type = self.genotype_list[random.randint(0, len(self.genotype_list))]
        else: # start!='', not initial step, just pass
            pass
        
        fitness = self.landscape  # landscape
        
        neigh_list = [self.start_type]  # this list will be used to store the self.start_type and all its neighbors
                                # (next step will select from the neighbors of self.start_type)
        
        neigh_list += self.creat_neighbor(self.start_type, 1) # get all 1-mutant neighbors of start_type
                
        neigh_fitness = {i : fitness[i] for i in neigh_list[1:]}  # creat a dictionary to store the genotypes and fitness of 
                                                            # the 1-mutant neighbors of self.start_type (start_type is not 
                                                            # included as it won't be selected for self.next_step) 
                                                            # key: genotype; value: fitness

        self.start_type_fitness = fitness[self.start_type] # get the fitness of self.start_type
        
        possible_neigh = {k : v for k, v in neigh_fitness.items() if v >= self.start_type_fitness} # get all possible   
                                                                    # neighbors for self.next_step (only if the fitness >= 
                                                                    # self.start_type fitness can be possible neighbors 
                                                                    # for self.next_step)
                                                                                        
        if possible_neigh:   # if possible_neigh is not empty, it means that possible genotypes for self.next_step exist;
                             # select the genotype for self.next_step according to their fitness
            pop_gens = list(possible_neigh.keys())   # get all the possible genotypes of self.next_step
            self.next_type = pop_gens[random.randint(0, len(pop_gens))] # randomly select one genotype as self.next_step
            self.next_type_fitness = fitness[self.next_type] # get the fitness of self.next_step
                                                                                
        else: # if possible_neigh is empty, it means that there if no option for self.next_step-->cannot evolve
              #-->reach local optima
            self.next_type = ''
            self.next_type_fitness = 0
        
    
    def repeat_random_walk(self):
        '''Simulate the random walk in the landscape'''
        
        self.walk_list = []  # this list will be used to store the steps of random walk
        
        fit = self.random_walk()  # run the random_walk method
        
        self.walk_list.append([self.start_type, self.start_type_fitness])  # the initial step
        
        
        while self.next_type and self.next_type_fitness:  # this means that while next_step!='' and fitness[next_step]!=0.0,  
                                                         #not reach local optima
            self.walk_list.append([self.next_type, self.next_type_fitness]) # the next_step
            self.start_type = self.next_type  # update the self.start_type and self.start_type_fitness
            self.start_type_fitness = self.next_type_fitness
            fit = self.random_walk() # run the random_walk method with the updated self.start_type
        
    
    def multiple_random_walk(self, times):
        
        self.multiple_walk_list = []
        
        for i in range(times):
            self.repeat_random_walk()
            self.multiple_walk_list.append(self.walk_list)
            self.start_type = self.genotype_list[random.randint(0, len(self.genotype_list))]
    
    
    def search_optimal(self):
        '''Search the genotype of local optimal'''
        
        self.optimal_list = []   # will be used to store the local optima
        for i in self.genotype_list:
            neighbor_fitness = [self.landscape[j] for j in self.creat_neighbor(i, 1)]  # the fitness of neigibors of i
            if self.landscape[i] > max(neighbor_fitness):  # if fitness of i is larger than all its neighbors, i is local optima
                self.optimal_list.append(i)   # append i into self.optimal_list
                
        
