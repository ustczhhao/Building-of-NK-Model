# -*- coding: utf-8 -*-

"""

Created on Aug 15, 2016

The K epistic sites are randomly selected from the total N sites. The model is additive.

Written in Python 3


"""




import numpy as np
from numpy import random
import itertools
from scipy import stats


class NK_model(object):
    
    def __init__(self, n, k):
        self.n = n
        self.k = k
        self.start_type = 0    # the initial given start genotype
        self.get_epistic_site()  # get the epistic network among the sites (including the site itself) 
        
        self.allelic_each_site = [[] for _ in range(self.n)]   # will be used to store the episitc combinations at each site
        self.fitness_list = [{} for _ in range(self.n)]   # will be used to store the fitness of epistic combinations at each
                                                          # site
        self.genotype_list = []   # will be used to store the genotypes generated (subset of all the genotypes)
        self.geno_fitness = []   # will be used to store the calculated fitness of genotypes generated
        self.landscape = {}  # will be used to store the generated landscape (subset of the whole landscape)
        
        self.repeat_adaptive_walk()
        self.search_different_allele()
        
        
    @staticmethod
    def generate_start_type(seq_length):
        '''generate a genotype as the start_type
        
        parameters: 
        ----------
            seq_length: int, the length of the wanted genotype sequence
        output: 
        -------
            str that represents the start_type
        '''
        
        a_list = []
        for i in range(seq_length):  # generate a list randomly composed of 0 and 1
            a_list.append(random.randint(0, 2))
            
        start_type_list = []   
        for i in a_list:   # convert the integer elements (0 and 1) in a_list into string elements ('0' and '1') and append 
            start_type_list.append(str(i))   # into start_type_list
            
        start_type = ''.join(start_type_list)  # join the elements in start_type_list into a string
        
        return start_type
    
    
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
    
        mutant_indices = list(itertools.combinations(indices, mutant_num)) # get all possible combinations of picking 
                                                                    # mutant_num intergers from indices (i.e. mutant sites)
    
        for i in mutant_indices:  # get all its mutant_num -mutant neighbors of the given genotype
            neigh_list = geno_list.copy()
            for j in i:
                if geno_list[j] == '0':
                    neigh_list[j] = '1'
                else:
                    neigh_list[j] = '0'
            
            total_neigh_list.append(''.join(neigh_list))
    
        return total_neigh_list
    
    
    def get_epistic_site(self):
        '''Generate the epistic networks among the sites'''

        self.epistic_site = []
        for i in range(self.n):
            all_site = [x for x in range(self.n) if x != i]  # get the indexes of all sites excluding i
            
            indices = set()   # will be used to store the indexes of epistic sites
            while len(indices) < self.k:
                indices.add(random.randint(0, len(all_site))) # randomly pick the indexes of k elements from all_site

            indices = list(indices)  # convert the set into list
            
            one_site = [i]   # first append the site under study into one_site (a little bit different from 0615 code)
            for j in indices:  # get the epistic sites of site i
                one_site.append(all_site[j])
            
            self.epistic_site.append(one_site) 
           
        
    @staticmethod     
    def get_epistic_combination(genotype, network):
        '''generate the epistic pairs in a given genotype according to the network
        
        parameters: 
        ----------
            genotype: str that repesents a genotype, 
            network: nested list that stores the epistic networks among the sites (self.epistic_site)
            
        output: 
        -------
            list that stores all the epistic combinations in the given genotype
        '''
        
        epistic_combination = []   # will be used to store the epistic combinations in the given genotype
        
        for j in range(len(genotype)):
            one_site = []   
            for i in network[j]:   # get all epistic alleles in one site (store in list one_site)
                one_site.append(genotype[i])
                
            epistic_combination.append(''.join(one_site))  # combine the epistic alleles and store into epistic_combination
                
        return epistic_combination
        
        
    @staticmethod      
    def cal_fitness(epistic_pair, fitness_list):
        '''calculate the fitness of a given genotype according to the fitness_list (here the epistic_pair is generate by
        the genotype, i.e. the epistic_combination in get_epistic_combination(genotype, network) method)
        
        parameters: 
        ----------
            epistic_pair: list that stores the epistic combinations in a genotype  
            fitness_list: list that stores the dictionary which indicates the fitness of each epistic combinations 
                          in each site
            
        output: 
        -------
            float, the calculated fitness of a genotype (the genotype which corresponded to epistic_pair)
        '''
        
        one_fitness = []    
        
        for i in range(len(epistic_pair)):  # fitness_list[i] is the fitness dictionary that corresponds to site i.  
            one_fitness.append(fitness_list[i][epistic_pair[i]])  # epistic_pair[i] is the epistic combination in site i (the
                                                                # key in the dictionary)
        i_fitness = sum(one_fitness)/len(one_fitness)  # calculate the fitness
            
        return i_fitness
        
     
    def adaptive_walk(self):
        '''Simulate one step  of adaptive walk along the landscape'''
        
        if self.start_type==0:  # self.start_type==0 (the initial given value), means initial step, randomly generates a 
                                # genotype as start_type (length is self.n)
            self.start_type = self.generate_start_type(self.n)
        else:  # start!='', not initial step, just pass
            pass
        
        neighbor_list = [self.start_type] # this list will be used to store the self.start_type and all its neighbors
                                        # (next step will select from the neighbors of self.start_type)
        
        neighbor_list += self.creat_neighbor(self.start_type, 1) # get all 1-mutant neighbors of start_type
                
        all_epistic_pair = [] # this list will be used to store all epistic combinations for genotypes in neighbor_list
        
        for i in neighbor_list:  # i is genotype in neighbor_list, self.epistic_site is network among the site
            all_epistic_pair.append(self.get_epistic_combination(i, self.epistic_site))
        
        
        allelic_set = []  # list that will be used to store all the epistic combinations in site 0, 1, 2, ..., self.n-1 
                        # (use set to exclude the repeat ones) for genotypes in neighbor_list
        
        for j in range(self.n): 
            pos_allelic_pair = [] # will be used to store all epistic combinations in site j for genotypes in neighbor_list
            for s in all_epistic_pair:
                pos_allelic_pair.append(s[j])   
            
            pos_allelic_set = list(set(pos_allelic_pair))  # convert the pos_allelic_pair to set then convert back to list to
                                                    # get the non-repeated epistic combinations in site j
                
            allelic_set.append(pos_allelic_set)  # append the non-repeated epistic combinations in site j into allelic_set to
                                        # get epistic combinations in all sites
            
        for s in range(self.n): # the length of allelic_set is self.n
            for i in allelic_set[s]: # i is the allelic combinations in site s, self.allelic_each_site is used to store the 
                                    #  episitc combinations at each site, self.allelic_each_site[s] is episitc combinations 
                                    # at site s that has already been generated
                if i not in self.allelic_each_site[s]: # if i is not in self.allelic_each_site[s], means epistic combination i
                                 # has not been generated, just append i into self.allelic_each_site[s] and generate a fitness  
                                 # for i and store into self.fitness_list[s](i will be the key)
                    self.allelic_each_site[s].append(i)
                    self.fitness_list[s][i] = random.random(1)[0]
            
        w_fitness = []   # will be used to store the calculated fitness of each genotype in neighbor_list
        
        for i in neighbor_list:
            if i not in self.genotype_list: # i not in self.genotype_list, means the fitness of genotype i has not been
                                            # calculated--> need to be calculated this time
                index_i = neighbor_list.index(i) # first get the index of i (this index is the same as the epistic 
                                                # combinations for genotype i in all_epistic_pair)
                i_fitness = self.cal_fitness(all_epistic_pair[index_i], self.fitness_list)
                w_fitness.append(i_fitness)
            else:  # otherwise it means that the fitness of genotype i has been calculated, just find it from the landscape
                w_fitness.append(self.landscape[i])
            
        neigh_fitness = dict(zip(neighbor_list[1:], w_fitness[1:]))  # creat a dictionary to store the genotypes and fitness of 
                                                            # the 1-mutant neighbors of self.start_type (start_type is not 
                                                            # included as it won't be selected for self.next_step) 
                                                            # key: genotype; value: fitness
    
        self.start_type_fitness = w_fitness[0]  # get the fitness of self.start_type
        
        possible_neigh = {k : v for k, v in neigh_fitness.items() if v >= self.start_type_fitness}  # get all possible   
                                                                    # neighbors for self.next_step (only if the fitness >= 
                                                                    # self.start_type fitness can be possible neighbors 
                                                                    # for self.next_step) 
        
        if possible_neigh:   # if possible_neigh is not empty, it means that possible genotypes for self.next_step exist;
                             # select the genotype for self.next_step according to their fitness
                

            pop_gens = []
            value = []

            for k, v in possible_neigh.items():
                pop_gens.append(k)  # get all the possible genotypes of self.next_step
                value.append(v)


            s_list = []

            for i in value:
                s = (i- self.start_type_fitness)/self.start_type_fitness
                s_list.append(s)

            sum_s_list = sum(s_list)

            p_list = []    

            for j in s_list:
                p_list.append(j/sum_s_list)

            # self.next_type = random.choice(pop_gens, 1, p_list)[0]
            next_type_index = list(np.random.multinomial(1, p_list)).index(1)

            self.next_type = pop_gens[next_type_index]

            self.next_type_fitness = possible_neigh[self.next_type]
                                                                                
        else: # if possible_neigh is empty, it means that there if no option for self.next_step-->cannot evolve
              #-->reach local optima
            self.next_type = ''
            self.next_type_fitness = 0   
           
    
        self.genotype_list.extend(neighbor_list)  # updating the self.genotype_list by storing all the genotype in neighbor_list
                                                # into self.genotype_list
            
        self.genotype_list = list(set(self.genotype_list))  # exclude the repeating genotypes in self.genotype_list
        
        landscape_subset = dict(zip(neighbor_list, w_fitness)) # the subset of landscape generated in this random_walk step
        
        self.landscape.update(landscape_subset)  # update self.landscape by the subset of landscape generated in this random
                                                 # _walk step
                
        
    def repeat_adaptive_walk(self):
        '''Simulate the random walk in the landscape'''
        
        self.walk_list = []  # this list will be used to store the steps of random walk
        
        fit = self.adaptive_walk()  # run the random_walk method
        
        self.walk_list.append([self.start_type, self.start_type_fitness])  # the initial step
        
        while self.next_type and self.next_type_fitness:  # this means that while next_step!='' and fitness[next_step]!=0.0,  
                                                         #not reach local optima
            self.walk_list.append([self.next_type, self.next_type_fitness]) # the next_step
            self.start_type = self.next_type  # update the self.start_type and self.start_type_fitness
            self.start_type_fitness = self.next_type_fitness
            fit = self.adaptive_walk() # run the random_walk method with the updated self.start_type 
        
        self.walk_genotype = list(set([i[0] for i in self.walk_list])) 
        self.new_walk_list = [[j, self.landscape[j]] for j in self.walk_genotype]
        self.new_walk_list.sort(key=lambda x: x[1])
        self.walk_genotype = [i[0] for i in self.new_walk_list]
 
        if len(self.new_walk_list) > 6:
            self.walk_genotype = self.walk_genotype[:6]
            self.new_walk_list = self.new_walk_list[:6]
    
    
    
    def search_different_allele(self):
        '''check all the alleles in each site for all genotypes involved in an adaptive walk(e.g. site 1: have '0' and '1',
        (means some genotypes have '0' in site 1 and some other have '1' in site 1),site 2: only have '0' (means all genotypes 
        involved in an adaptive walk have '0' in site 2))'''
        
        walk_geno_list = []  
        for i in self.walk_genotype:   # first convert the genotypes into lists (each genotype will have a list) and store in
            walk_geno_list.append(list(i))  # walk_geno_list
            
        self.total_allele_list = []  # will be used to store the alleles present in each site for all genotypes involved in 
                                    # an adaptive walk
        for i in range(self.n): # for each site
            allele_list = []  # will be used to store all alleles in site i
            for j in walk_geno_list:  
                allele_list.append(j[i])
                
            allele_list = list(set(allele_list))#convert allele_list into a set to exclude the repeated one and then 
                                                # back to list
            self.total_allele_list.append(allele_list) # append the alleles in site i into self.total_allele_list
                
        
        self.shared_allele_indice = [] #will be used to store the indexes of shared alleles for the genotypes in an 
                                    # adaptive walk
        self.diff_allele_indice = [] #will be used to store the indexes of different alleles for the genotypes in 
                                    # an adaptive walk
        
        for i in range(self.n):
            if len(self.total_allele_list[i])>1: # this means site i has more than one alleles-->different allele,append i into
                self.diff_allele_indice.append(i) # self.diff_allele_indice
            else: # else means site i only have one allele for all genotypes--> shared allele, append i into self.shared_allele
                self.shared_allele_indice.append(i) # # _indice
        
        self.ordered_diff_allele_indice = []  # will be used to store the different alleles accroding to the mutation order
        for i in range(len(self.walk_genotype)-1): 
            change_indice = self.compare_list(list(self.walk_genotype[i]),list(self.walk_genotype[i+1]))[0]
                 # compare the neighboring two genotypes in an adaptive walk to get the order of sites mutated
            self.ordered_diff_allele_indice.append(change_indice) #append the sites mutated into self.ordered_diff_allele_indice
        
        
        self.shared_allele = ''  # # will be used to store the shared sequences for all genotypes
        
        for i in self.shared_allele_indice: 
            self.shared_allele+=self.walk_genotype[0][i] # # get the shared sequence by all genotype (use the first genotype in 
                                                      # adaptive walk to construct, also can use other genotypes in adatptive 
                                                      # walk because this sequence is shared by all genotypes
                    
        self.genotype_subset = []  # this list will store the subsets of genotypes in adaptive walk (here subset means that 
                         # only the different alleles are present, excluding the shared alleles for all genotypes in the walk)
        self.landscape_subset = {}  # this is what we would to get, the subset of landscape which only contain the information
                                    # of different alleles
        
        for i in self.walk_genotype: # get different allele sequence for each genotype in an adaptive walk
            s = '' 
            for j in self.diff_allele_indice: # the different allele sequence in genotype i
                s+=i[j]
            
            self.genotype_subset.append(s)  # append the subsets of genotypes into self.genotype_subet (will use the different
                                         # allele sequence to represent the whole genotype in an adaptive walk)
            self.landscape_subset[s] = self.landscape[i] # update the subset of landscape (key is the different allele 
                                                         # sequence in genotype i, value is the fitness of genotype i)
            
            
      
    def creat_landscape_subset(self):
        '''Generate the subset of lanscape which contains all possible combinations in different allele sites'''
                      
        all_possible_combination = [''.join(i) for i in itertools.product('01', repeat = len(self.diff_allele_indice))]
        # generate all possible allele combinations in different allele sites
        
        already_find = list(self.landscape_subset.keys()) # this is allele combinations that have already been got
        
        need_to_find_com = []  # will be used to store the allele combinations that haven't been got

        for i in all_possible_combination: 
             if i not in already_find:  # get all allele combinations that haven't been got
                need_to_find_com.append(i)

        need_to_find_com_list = []
        for i in need_to_find_com:   # convert each of the allele combinations in need_to_find_com to list and 
            need_to_find_com_list.append(list(i))  # store in need_to_find_com_list

            
        self.shared_allele_list = list(self.shared_allele) # convert the shared sequence (str) into a list

        orginal_genotype = []  # will be used to store the whole genotype constructed by the shared sequence and different 
                                # allele combinations (will be then used to find the fitness in self.landscape)
        for i in need_to_find_com_list: # only the different allele combinations that haven't been found (i.e. in need_to_find
                                        # _com_list)  will be converted back to the whole genotype sequence
            for j in self.shared_allele_indice:
                index_j = self.shared_allele_indice.index(j)
                i.insert(j, self.shared_allele_list[index_j]) # insert the alleles in shared_seq according to i according to the
            orginal_genotype.append(''.join(i)) # indices of shared_sequence in whole genotype sequence

            
        orginal_genotype_fitness = [] #will be used to store the fitness of each construced whole genotypes in orginal_genotype

        all_geno = list(self.landscape.keys()) # first get all whole genotype sequences in self.landscape

        for i in orginal_genotype:
            if i in all_geno:  # i in all_geno--> i already in self.landcape, get the fitness of i directly
                orginal_genotype_fitness.append(self.landscape[i])
            else:   # else i not in all_geno--> i not in self.landscape, need to calculate the fitness of i
                allelic_comb_i = self.get_epistic_combination(i, self.epistic_site) # generate the epistic pairs in i 
                                                                     # according to the network (i.e. self.epistic_site)
                    
                for j in range(self.n): # check whether the all_comb_i[j] have already been generated and given a fitness in 
                                        # site j
                    if allelic_comb_i[j] not in self.allelic_each_site[j]: # if not, just append the new combination into self.
                        self.allelic_each_site[j].append(allelic_comb_i[j]) # allelic_each_site[j] and given a fitness
                        self.fitness_list[j][allelic_comb_i[j]] = random.random(1)[0]

                fitness_i = self.cal_fitness(allelic_comb_i, self.fitness_list) # calcualte the fitness of genotype i

                orginal_genotype_fitness.append(fitness_i) # append the fitness of genotype i into orginal_genotype_fitness
                self.landscape[i] = fitness_i  # also update self.landscape with genotype i

        updated_landscape_subset = dict(zip(need_to_find_com, orginal_genotype_fitness)) # the updated_landscape_subset 
                                    # (need_to_find_com stores the keys, orginal_genotype_fitness stores the fitness,
                                    # need_to_find_com has the same order with orginal_genotype, orginal_genotype has the 
                                    # same order with orginal_genotype_fitness, so need_to_find_com has the same order with
                                    # with orginal_genotype_fitness)

        self.landscape_subset.update(updated_landscape_subset) # updated the self.landscape_subset with updated_landscape_subset.
    
    
    @staticmethod
    def compare_list(list1, list2):
        '''compare the elements in the two lists and return the indexes of different elements
        parameters: 
        ----------
            two lists with the same length (list1 and list2) to be compared
            
        output: 
        -------
            list that stores the indexes of the different elements in list1 and list2'''
        
        index_list = [] 
        
        for i in range(len(list1)):
            if list1[i] != list2[i]:
                index_list.append(i)
            else:
                pass
        
        return index_list
    
    
    def get_one_mutant(self):
        '''Get the one-site mutants (i.e.genotypes corresponding to W1, W2, W3, W4 and W5 in 
        Draghi 2013 paper Equation 5~8 )(the genotypes in adaptive walks are W12, W123, W1234 and W12345)'''
        
        self.ancestor = self.genotype_subset[0]  # first get the ancestor (the first genotype in an adaptive walk, self.
                         # genotype_subset uses the different allele sequences to represents the genotypes in adaptive walk)
            
        self.ancestor_fitness = self.landscape_subset[self.ancestor]
        
        ancestor_list = list(self.ancestor) # then convert the genotype of ancestor to list
        
        self.subset_change_index = []  # will be used to store the order of sites mutated in an adaptive walk
        
        for i in range(len(self.genotype_subset)-1): 
            change_index = self.compare_list(list(self.genotype_subset[i]),list(self.genotype_subset[i+1]))[0]
                 # compare the neighboring two genotypes in an adaptive walk to get the order of sites mutated
            self.subset_change_index.append(change_index) # append the sites mutated into self.subset_change_index
            
        self.one_mutant  = []  # will be used to store the one-site mutants (i.e. genotypes of W1, W2, W3, ...)
        
        for i in self.subset_change_index: # self.subset_change_index stores the indexes of mutated sites according to the 
                                            # mutated order
            start = ancestor_list.copy()
            if ancestor_list[i] =='0':
                start[i] = '1'
            else:
                start[i] = '0'
            
            self.one_mutant.append(''.join(start)) # get the one-site mutants
            
        self.one_mutant_fitness = []    
        for j in self.one_mutant: # get the fitness of the one_site mutants (i.e. W1, W2, W3, W4 and W5)
            self.one_mutant_fitness.append(self.landscape_subset[j])
        
        self.w = []
        
        for i in self.one_mutant_fitness:
            self.w.append(i - self.ancestor_fitness)
            
            
    @staticmethod
    def get_mutated_genotype(genotype, mutate_site):
        
        '''Get the mutated genotype based on the given genotype
        parameters: 
        ----------
            genotype: str that represents the given genotype
            mutate_site: list, the list that stores the indexes of sites that has mutated
            
        output: 
        -------
            str, the mutated genotype
        '''
        geno_list = list(genotype) 
        
        mutate_list = geno_list
        
        for i in mutate_site:
            if geno_list[i]=='0':
                mutate_list[i] = '1'
            else:
                mutate_list[i] = '0'
            
        mutate_geno = ''.join(mutate_list)
            
        return mutate_geno
    
    
    def construct_pairwise_index(self):
        '''Coustruct the pairwise index (i.e. pairwise indexes corresponding to W12, W13, W14, W15, W23, W24, W25, W34, W35, 
        W45, i.e, substitution 12, 13, 14, ..., 35, 45)'''
        
        self.pairwise_site = []  # will be used to store the pairwise substitution sites (combined sites corresponding to 
                                # substitution 12, 13, 14, ..., 34, 35 and 45)
        self.pairwise_site_index = []  # will be used to store the index combinations of the pairwise substitutions 
                                    # (i.e, 12, 13, 14, ..., 34, 35 and 45 in substitution 12, 13, 14, ..., and 45)
        for i in self.subset_change_index:
            index_i = self.subset_change_index.index(i)  # get the index of i in self.subset_change_index
            for j in self.subset_change_index[(index_i+1):]:
                index_j = self.subset_change_index.index(j) # get the index of j in self.subset_change_index
                self.pairwise_site.append([i, j])  # store the pairwise substitution site combination
                self.pairwise_site_index.append([index_i, index_j]) # store thee index combinations of pairwise substitutions
                
                
    def get_pairwise_genotype(self):
        '''Generate the pairwise substitution genotypes (genotypes corresponding to W12, W13, W14, W15, ..., W34, W35 and W45)'''
        
        self.pairwise_genotype = []   # will be used to store the pairwise genotypes
        self.pairwise_genotype_fitness = [] # will be used to store the fitness of pairwise genotypes
        self.pairwise_w = [] # will be used to store the calculated w of the pairwise genotypes
        
        for i in self.pairwise_site: # i is the combined pairwise substitution sites
            constructed_genotype = self.get_mutated_genotype(self.ancestor, i)  # pairwise substitution genotype
            self.pairwise_genotype.append(constructed_genotype)
            self.pairwise_genotype_fitness.append(self.landscape_subset[constructed_genotype]) # fitness of pairwise genotype
            self.pairwise_w.append(self.landscape_subset[constructed_genotype]-self.ancestor_fitness) # w of pairwise genotype
    
    
    @staticmethod
    def plus_certain_element(list1, index_list):
        '''Plus the elements in list1 according to the indexes stored in index_list
        parameters: 
        ----------
            list1: list
            index_list: list, the list that stores the indexes of elements in list1 that would be multiplied
            
        output: 
        -------
            float, the multiplied results of elements in list1 whose indexes stored in index_list'''
        
        m = 0
        
        for i in index_list:
            m+=list1[i]
            
        return m
    
    
    def calculate_e_effect(self):
        '''Calculate the e of each pairwised genotypes (i.e. e12, e13, e14, ..., e34, e35 and e45.
        e12 = W12-W1W2, ..., e35 = W35-W3W5, e45=W45-W4W5)'''
        
        self.pairwise_e_effect = [] # will be used to store the calculated e of each pairwised genotypes
        
        for i in self.pairwise_w:
            index_i = self.pairwise_w.index(i) # get the index of i 
            corr_one_mutant_index = self.pairwise_site_index[index_i] # get the indexes corresponding one site mutates(i.e., 1
                                                                     # and 2 for W1 and W2, 3 and 4 for W3 and W4) 
            plus_one_mutant_fitness = self.plus_certain_element(self.w, corr_one_mutant_index) # get the multiply results
                                                                        # (i.e., W1*W2, W1*W3,...)
            print('I', i)
            print('PLUS', plus_one_mutant_fitness)
            self.pairwise_e_effect.append(i-plus_one_mutant_fitness) # calculate e and store in self.pairwise_e_effect
