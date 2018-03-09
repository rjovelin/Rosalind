# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:36:51 2015

@author: Richard
"""



def eulerian_cycle(input_file, outputfile):
    '''
    (file) -> file
    
    Given an Eulerian directed graph in the form of an adjacency list,
    return an Eulerian cycle in this graph.
    '''
    
    #open file for reading
    infile = open(input_file, 'r')
    # make a dictionnary to store the nodes
    graph = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            first_node = line[0:line.index(' ')]
            second_node = line[line.index('> ')+2:]
            if ',' in second_node:
                second_node = second_node.split(',')
                graph[first_node] = second_node
            else:
                graph[first_node] = [second_node]
    # close file
    infile.close()
    
    # count the number of edges
    edges = 0
    for node in graph:
        edges += len(graph[node])
        
    # make a list of nodes already explored  
    explored = []
    
    # form an initial cycle starting at '0'
    cycle = []
    firstnode = '0'
    # add node to cycle
    cycle.append(firstnode)
    secondnode = graph[firstnode][0]
    cycle.append(secondnode)
    explored.append(secondnode)
    # remove secondnode from graph
    graph[firstnode].remove(secondnode)
    firstnode = secondnode
    
    # iterate until the cycle has the expected number of edges
    while len(cycle) < edges:
        
        # build the cycle until its stuck and return to the first node
        while firstnode in graph and graph[firstnode] != []:
            secondnode = graph[firstnode][0]
            graph[firstnode].remove(graph[firstnode][0])
            cycle.append(secondnode)
            explored.append(secondnode)
            firstnode = secondnode
            
        # read the cycle and find the node with nodes unexplored
        # start from that node and go through cycle, populating the new cycle
        # until that node is reached. add that node also to new cycle
        # resume looking for nodes in graph until its stuck
        
        not_found = True
        new_cycle = []
        for i in range(len(cycle)):
            if not_found:
                firstnode = cycle[i]
                for node in graph[firstnode]:
                    if node not in explored:
                        firstnode = cycle[i]
                        if cycle[0] == cycle[-1]:
                            secondnode = cycle[i:-1]
                        else:
                            secondnode = cycle[i:]
                        new_cycle.extend(secondnode)
                        new_cycle.extend(cycle[:i+1])
                        cycle = new_cycle
                        not_found = False
                        break
        
    # open file for writing
    newfile = open(outputfile, 'w')
    # write content of cycle to file
    for item in cycle[:-1]:
        newfile.write(item + '->')
    newfile.write(cycle[-1])
    #close file
    newfile.close()


# some explanations for the following graph:
# 0 -> 3
#1 -> 0
# 2 -> 1,6
# 3 -> 2
# 4 -> 2
# 5 -> 4
# 6 -> 5,8
# 7 -> 9
# 8 -> 7
# 9 -> 6
                
# form initial cycle by starting at a random node in the adjacency list
# Let's start at 0: build the cycle, we get 0 -> 3 -> 2
# since "2" has more ways to go, choose a random one (start at the 1st doesnt matter)
# now cycle is 0 -> 3 -> 2 -> 1
# one more iteration from graph and the cycle is stuck at 0 -> 3 -> 2 -> 1 -> 0
# Since the cycle is shorter than expected 
# traverse through your cycle and search a node with unexplored edges:
#0 -> no
#3 -> no 
#2 -> yes
# 2 is the new starting point now.
# build a new cycle by starting at 2 and traversing through the old cycle
# until 2 is reached, keep 2 in the new cycle 2- > 1 -> 0 -> 3- > 2
# select the other unexplored node connected to 2
# keep building the new cycle until its stuck again and repeat   
              
    
    