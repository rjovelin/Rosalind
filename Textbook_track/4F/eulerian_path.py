# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:06:11 2015

@author: Richard
"""

def eulerian_path(input_file, outputfile):
    '''
    Given a directed graph that contains an Eulerian path, 
    where the graph is given in the form of an adjacency list,
    return an Eulerian path in this graph.
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # make a dictionnary to store graph
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
    
    # make a list to store the eulerian cycle
    cycle = []
     
    # find starting node and end nodes
    first_nodes = [node for node in graph]
    second_nodes = []
    for node in graph:
        for item in graph[node]:
            second_nodes.append(item)
    
    for node in first_nodes:
        # starting node outdegree > indegree
        if len(graph[node]) > second_nodes.count(node):
            firstnode = node
    for node in second_nodes:
        # end node indegree < outdegree
        if node in graph:
            if len(graph[node]) < second_nodes.count(node):
                lastnode = node
        else:
            lastnode = node
    print(firstnode, lastnode)
     
    # create a link between the lastnode and the first node to form a 
    # graph that contains an eulerian cycle
    if lastnode in graph:
        graph[lastnode].append(firstnode)
    else:
        graph[lastnode] = [firstnode]
        
    # count the number of edges
    edges = 0
    for node in graph:
        edges += len(graph[node])
    
    # form an eulerian cycle starting at the firstnode    
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
                    if not node in explored:
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
        
    
    # transform the eulerian cycle into an eulerian path from starting node
    # to end node
    euler_path = []
    end = cycle.index(lastnode)
    euler_path.extend(cycle[end+1:-1])
    euler_path.extend(cycle[:end+1])
    
    # convert list to string
    euler_path = '->'.join(euler_path)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    newfile.write(euler_path)
    newfile.close()
    


    
    
    
    
    
    