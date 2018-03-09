# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 23:46:47 2015

@author: Richard
"""

def overlap_graph(input_file, outputfile):
    '''
    (file, file) -> file
    Given a collection of Patterns k-mers in input file, return the overlap graph
    Overlap(Patterns) in the form of a adjacent list (ie. node1 -> node2)
    Precondition: the sequence overlap between any overlapping substring is the same
    and is equal to length of substring - 1
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # create list to hold sequences
    sequences = [line.rstrip() for line in infile if line != '\n']
    # close file
    infile.close()

    # make a dictionnary to store the source: target node pairs
    graph = {} 
        
    # go through the list of sequences
    # for each sequence, ask if another sequence's prefix == focal sequence's suffix
    # record the focal seq in dict and add the other seq in the list of target nodes
    
    # make a set of outcoming nodes
    outcoming = set()
    # make a set of incoming nodes
    incoming  = set()
    
    # go through each sequence in the file
    for i in range(len(sequences)):
        # compare this sequence to every other sequence
        for j in range(len(sequences)):
            # do not consider auto loop
            if i != j:
                # check if the sequence's suffix is the other sequence's prefix
                if sequences[i][1:] == sequences[j][0:-1]:
                    # check if the outcoming node is in dict
                    if sequences[i] in graph:
                        # add the incoming node
                        graph[sequences[i]].append(sequences[j])
                        # add the incoming node to set
                        incoming.add(sequences[j])
                    else:
                        # add both the incoming and outcoming nodes
                        graph[sequences[i]] = [sequences[j]]
                        # add nodes to respective sets
                        incoming.add(sequences[j])
                        outcoming.add(sequences[i])
        
        
    # open file for writing
    newfile = open(outputfile, 'w')        
        
    # check tthat all nodes are used
    missing = set()
    for seq in sequences:
        if seq not in incoming.union(outcoming):
            missing.add(seq)
    # print warning if any sequence is missing
    if len(missing) != 0:
        print('warning: {0} sequences are missing'.format(len(missing)))
    else:
        # print to screen the overlapping graph (node -> node)
        for node in graph:
            if len(graph[node]) == 1:
                newfile.write(node + ' -> ' + graph[node][0] + '\n')
            else:
                for seq in graph[node]:
                    newfile.write(node + ' -> ' + seq + '\n')
    
    newfile.close()
                
                
                
                

    
    
    
    
    
    
    
    
    