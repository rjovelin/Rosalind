# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 23:25:20 2015

@author: Richard
"""

def shortest_superstring(input_file):
    '''
    (file) -> str
    Given at most 50 DNA strings in FASTA in the input file
    with length <= 1 Kb representing reads derived from the same strand
    of a single linear chromosome return the shortest superstring containing
    all the given strings (ie. = reconstructed chromosome).
    Precondtion: There exists a unique way to reconstruct the entire
    chromosome from these reads by gluing together pairs of reads that
    overlap by more than half their length.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # create dict to hold the fasta seq
    fasta = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                name = line[1:]
                fasta[name] = ''
            else:
                fasta[name] += line
    #close file
    infile.close()
    
    # find the first sequence
    # make a list of sequences
    unique_first = True
    sequences = [fasta[gene] for gene in fasta]
    # go through the sequences and find the seq with unique first half
    not_first = []
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i != j:
                if sequences[i][:len(sequences[i]) // 2] in sequences[j]:
                    not_first.append(sequences[i])
    if len(not_first) == len(sequences) - 1:
        for seq in sequences:
            if seq not in not_first:
                first_seq = seq
    else:
        print('found ambiguous first sequences')
        unique_first = False
    
    if unique_first == True:
        chromo = first_seq
        sequences.remove(first_seq)
    
    # find and add the next sequence until all the sequences have been used
    while len(sequences) != 0:
        for seq in sequences:
            if seq[: len(seq) //2] in chromo:
                new_seq = seq
                sequences.remove(new_seq)
                print(len(sequences))
                break
        stop = chromo.index(new_seq[: len(new_seq) //2])
        chromo = chromo[:stop]
        chromo += new_seq
    return chromo

    
    
    

            
        
        
        
        
        
    
    
    
    
    