# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:11:06 2015

@author: Richard
"""


# import the module Entrez to have access to Genbank
from Bio import Entrez


def genus_nucleotide(input_file):
    '''
    (str, str, str) -> int
    Use Biopython to return the number of nucleotide entries in Genbank corresponding to Genus
    in Genbank from date1 to date2 given in the input file
    '''
        
    # open file for reading
    infile = open(input_file, 'r')
    # create a list to store the paramters needed to resolve 
    parameters = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            parameters.append(line)
    # close file after reading
    infile.close()
    
    # get the paramter values    
    genus, date1, date2 = parameters[0], parameters[1], parameters[2]
    
    # create a handle
    handle = Entrez.esearch(db = 'nucleotide',\
    term = "{0}[Organism] AND {1}[Publication Date]:{2}[Publication Date]".format(genus, date1, date2))
    # create a record by reading the handle, in a form of a dict
    record = Entrez.read(handle) 
    # access information, including counts for record object with keywords
    
    return record['Count']
    
    

    


    
    
    
    
    
    
    