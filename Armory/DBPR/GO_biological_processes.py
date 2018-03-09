# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:45:59 2015

@author: Richard
"""

# uses Biopython to retrieve Biological Processes GO annotations from SwissProt



def GO_biological_processes(input_file):
    '''
    (file) -> list
    Given the Uniprot IDof a protein in the input file, return the list of
    Gene Ontologies corresponding to Biological Processes annotated in Swissprot
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    prot_ID = infile.readline().rstrip()
    
    # make a handle to read the swissprot infor from Expasy
    handle = ExPASy.get_sprot_raw(prot_ID)
    # create a a record object with all information about the protein ID
    record = SwissProt.read(handle)
    
    # get a list of references from other databases
    outside_ref = record.cross_references
    # get the information from the GO database
    GO = [ref for ref in outside_ref if 'GO' in ref[0]]
    # parse the list of GO terms to grab only the GO biological processes
    GO_BP = []
    for item in GO:
        if item[2].startswith('P:'):
            process = item[2][2:]
            GO_BP.append(process)
    
    # print biological processes
    for process in GO_BP:
        print(process)
    

        
    


if __name__ == '__main__':
    # import the ExPasy module
    from Bio import ExPASy
    # import the Swissprot module
    from Bio import SwissProt
    
    