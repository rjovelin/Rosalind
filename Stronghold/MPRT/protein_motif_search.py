# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:26:16 2015

@author: Richard
"""

# import module to open url
import urllib
# import module to search text with regular expression
import re

def unipro_seq_to_fasta(input_file):
    '''
    (file) -> dict
    Given UniProt Protein Database access IDs, return a dictionnary with the 
    pritein ID : protein sequence pairs
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    
    protein_names = []
    
    # read file and get the protein ID
    for line in infile:
        line = line.rstrip()
        if line != '':
            # get the protein ID and store it in a list
            protein_ID = line
            # parse the protein ID if it contains '_' because the online ID
            # is without underscores
            if '_' in protein_ID:
                protein_ID = protein_ID[:protein_ID.index('_')]
            protein_names.append(protein_ID)
    # close file
    infile.close()
    
    # create a dictionnary to store the protein ID : protein sequence pairs
    protein_sequences = {}
    
    # grab sequences from the unipro database
    for protein_ID in protein_names:
        # build url to open
        url = 'http://www.uniprot.org/uniprot/' + protein_ID + '.fasta'
        # open the page and read content
        page = urllib.request.urlopen(url)
        # read the sequence header
        header = page.readline().rstrip()
        # grab the sequence
        sequence = page.read()
        # convert the byte string to string
        sequence = sequence.decode('utf-8')
        # replace end of line characters
        sequence = sequence.replace('\n', '')
        # populate dictionnary
        protein_sequences[protein_ID] = sequence
    
    return protein_sequences        
    

def get_unipro_protein_IDs(input_file):
    '''
    (file) -> tuple
    Given a collection of unipro potein IDs in input file, return a tuple 
    with 2 lists: a list of protein names without any underscore
    and a list of unmodified protein IDs.
    Elements of each list are corresponding for a same index
    '''
    # initialize lists
    protein_names = []
    protein_IDs = []
    
    # open file for reading
    infile = open(input_file, 'r')
    
    # read file and get the protein ID and protein name
    for line in infile:
        line = line.rstrip()
        if line != '':
            # get the protein ID and store it in a list
            ID = line
            # parse the protein ID if it contains '_' because the online ID
            # is without underscores
            if '_' in ID:
                name = ID[:ID.index('_')]
            else:
                name = ID
            # populate each list
            protein_names.append(name)
            protein_IDs.append(ID)
    # close file
    infile.close()

    return protein_names, protein_IDs
    

def protein_motif_finder(protein_sequences):
    '''
    (dict) -> dict
    Given a dictionnary of string protein_name : string protein_sequence pairs,
    return a dictionnary with {string protein_name : list of lists for each
    match of the N-glycosil motif in protein sequence} pairs.
    Each inner list contains the starting index, and ending index (non-inclusive),
    and the match sequence, of the matching pattern in protein sequence
    '''

    # make list of tuples each containing seq name and string sequence
    protein = [(name, sequence) for name, sequence in protein_sequences.items()]
    
    # make a dictionnary with index of match object : match object pair
    patterns = {}
    
    # loop over the protein list,
    # search all, overlapping patterns in prot seq
    # re.match search for patterns at the begining of the sequence only
    for prot in protein:
        # loop over the sequence to find matches at the beginning of sequence
        for i in range(len(prot[1])):
            # look for matches at the begining the sequence
            motif = re.match('N[^P][ST][^P]',prot[1][i:])
            # if match, grab the start and end positions, and the pattern
            if motif != None:
                # use protein name as key
                if prot[0] in patterns:
                    patterns[prot[0]].append([i, i+ len(motif.group()), motif.group()])
                else:
                    patterns[prot[0]] = [[i, i+ len(motif.group()), motif.group()]]
    
    return patterns
    


def pattern_positions_printer(input_file):
    '''
    (file) -> None
    Given a list of UniProt Protein Database access IDs, 
    Print for each protein possessing the N-glycosylation motif:
    the protein ID on one line, followed by the list of starting indices
    of the motif on the next line.
    Use 1-base indexing
    '''
    
    # get the protein sequences for each protein ID
    protein_sequences = unipro_seq_to_fasta(input_file)
    
    # find the motif in each sequence
    patterns = protein_motif_finder(protein_sequences)
    
    # make lists containing the protein names without underscore
    # as on the uniprot website, and the protein IDs given as input
    protein_names, protein_IDs = get_unipro_protein_IDs(input_file)
    
    # loop over the protein names which are keys in protein_sequences
    for i in range(len(protein_names)):
        # check that the protein has a motif
        if protein_names[i] in patterns:
            # use protein ID as input
            print(protein_IDs[i])
            # print starting index, adjusting for 1-base indexing
            for motif in patterns[protein_names[i]]:
                start = motif[0] + 1
                print(start, end = ' ')
            print('')
                
            
    
    