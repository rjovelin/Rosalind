# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 22:49:09 2019

@author: rjovelin
"""


def ReadGraphFromFile(PbFile):
    '''
    (file) -> dict    

    Returns a dictionary representing graph in PbFile 
    '''
    
    infile = open(PbFile)
    graph = {}
    for line in infile:
        if '->' in line:
            line = line.split('->')
            parent = int(line[0].strip())
            children = set(int(i.strip()) for i in line[1].split(',')) 
            assert parent not in graph
            graph[parent] = children
    infile.close()
    return graph



def FindUnbalancedNodes(Graph):
    '''
    (dict) - > dict
    
    Return a dictionary with unbalanced nodes and their count of incoming
    and outgoing edges. Dictionary has necesarily only 2 nodes
    '''
    
    # count incoming and outgoing edges for every node in graph
    C = {}
    
    # make a list with all nodes
    L = list(Graph.keys())
    for i in Graph.values():
        for j in i:
            L.append(j)
    L = list(set(L))
    
    for i in L:
        C[i] = {}
        C[i]['in'] = 0
        C[i]['out'] = 0
        
    # count outgoing edges
    for i in Graph:
        C[i]['out'] += len(Graph[i])
    for i in Graph:
        for j in Graph[i]:
            C[j]['in'] += 1

    # find unbalanced nodes
    nodes = [i for i in C if C[i]['in'] != C[i]['out']]
    
    D = {}
    for i in nodes:
        D[i] = C[i]
    return D
    
    
def AddEdge(Graph, Unbalanced):
    '''
    (dict, dict) -> dict, list
    
    Returns a modified and balanced Graph in which an edge in added
    between unbalanced nodes and returns a list with the balanced nodes representing
    that directed edge
    Unbalanced has necessarily 2 nodes
    '''

    # find node missing an incoming edge 
    for i in Unbalanced:
        if Unbalanced[i]['in'] < Unbalanced[i]['out']:
            Incoming = i
    # find node with 
    for i in Unbalanced:
        if i != Incoming:
            Outgoing = i
    
    # modify Graph by adding a directed edge from Outgoing to Incoming 
    if Outgoing not in Graph:
        Graph[Outgoing] = {Incoming}
    else:
        Graph[Outgoing].add(Incoming)
    
    edge = [Outgoing, Incoming]
    
    return Graph, edge    


def FromCycleToGraph(cycle, edge):
    '''
    (list, list) -> list
    
    Take a eularian cycle and the edge added to graph to obtain cycle
    and return a eularian path by breaking that edge in cycle and reformating cycle
    '''

    # break edge in cycle and reformat cycle so that it becomes a eularian path
    
    for i in range(len(cycle)):
        if cycle[i] == edge[0] and cycle[i+1] == edge[1]:
            break
    path = cycle[i+1: -1] + cycle[: i+1]
    return path





# use breadth-first-search to find all paths from f to goal
def FindCycle(tree, f, goal):
    '''
    (dict, str, str) -> list
    Take a dictionary representing with interactions among node, and return
    a list os lists of nodes representing the path betwen nodes f and goal
    '''
    # create a queue to add the discovered nodes that are not yet visited
    # keep track of the path visited to reach discovered node
    
    L = []
    
    queue = [(f, [f])]
    visited = []
    while len(queue) != 0:
        #print('queue', queue)
        #print('L', L)
        # get a node to visit. keep track of the path up to node
        (node, path) = queue.pop(0)
        #print(node, path)
        if node not in visited:
            visited.append(node)
            #print('visited', visited)
            # add node children to queue
            for child in tree[node]:
                #print('child', child)
                if child == goal:
                    path.append(child)
                    L.append(path)
                else:
                    queue.append((child, path + [child]))
    return L


def TrackVisited(Graph):
    '''
    (dict) -> dict    
    Take a dictionary representring edges between nodes in Graph
    and return a dictionary with number of ingoing and outgoing edges for each node
    '''
    
    Tag = {}
    for i in Graph:
        # count outgoing edges:
        if i not in Tag:
            Tag[i] = {}
        Tag[i]['out'] = len(Graph[i])
        Tag[i]['in'] = 0
    for i in Graph:
        for j in Graph[i]:
            if j in Tag:
                Tag[j]['in'] += 1
    return Tag
    

def SelectCycle(Paths, Tag, C):
    '''
    (list, dict, list) -> list
    
    Take a list of list cycles, a dictionary with number of unused edges for
    each node in graph and a list of previously explored cycles and return a 
    cycle not previously explored and for which all nodes have unused incoming
    and outgoing edges
    '''
    
    L = []
    
    for cycle in Paths:
        Unexplored = True
        if cycle not in C:
            for i in cycle:
                if Tag[i]['out'] <= 0 or Tag[i]['in'] <= 0:
                    Unexplored = False
            L.append((cycle, Unexplored))
    for i in L:
        if i[1] == True:
            return i[0]
            
            
def UpdateVisited(cycle, Tag):
    '''
    (list, dict) -> dict
    Update Tag dictionary by removing 1 edge from the count of available outgoing
    edges for each node in cycle
    '''
    
    # don't remove an edge twice
    for i in range(len(cycle)):
        if i == 0:
            Tag[cycle[i]]['out'] -= 1
        elif i == len(cycle)-1:
            Tag[cycle[i]]['in'] -= 1
        else:
            Tag[cycle[i]]['out'] -= 1
            Tag[cycle[i]]['in'] -= 1
        if Tag[cycle[i]]['out'] < 0 or Tag[cycle[i]]['in'] < 0:
            print(cycle[i], Tag[cycle[i]]['out'], Tag[cycle[i]]['in'], cycle)
    return Tag



def FindStartingNode(cycle, Tag):
    '''
    list, dict) -> int    
    Take a list of nodes representing a cycle in graph and a dictionary with
    number of unused incoming and outgoing edges and return a node in cycle
    with remining unused incoming and outgoing edges
    '''
    
    Newstart = None
    for i in cycle:
        # check that node has outgoing edges unexplored
        if Tag[i]['in'] >= 1 and Tag[i]['out'] >= 1:
            Newstart = i
            break
    return Newstart
            
def EulerianCycle(Graph):
    '''
    (dict) -> list
    Return the eularian cycle in Graph
    '''
    
    # make a list with all explored cycles
    C = []
    
    # keep track of available outgoing edges
    Tag = TrackVisited(Graph)

    # form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
    # use BFS to find all cycles
    
    paths = FindCycle(Graph, 0, 0)
        
    # select a cycle for which all nodes have unexplored edges
    cycle = SelectCycle(paths, Tag, C)
    # update list of explored cycles
    C.append(cycle)
    
    # update used edges
    Tag = UpdateVisited(cycle, Tag) 
        
    # make a list of nodes with unused outgoing edges    
    L = [i for i in Tag if Tag[i]['out'] > 0]
    
    #  while there are unexplored edges in Graph
    while len(L) != 0:
        # find new starting node with unexplored edges
        Newstart = FindStartingNode(cycle, Tag)    
        if Newstart == None:
            print('no newstart', 'remaining nodes', len(L))
            return None
               
        # form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
        cyclePrime = cycle[cycle.index(Newstart):-1] + cycle[:cycle.index(Newstart)+1]
        # find new cycle from Newstart node
        paths = FindCycle(Graph, Newstart, Newstart)
        walk = SelectCycle(paths, Tag, C)
        # update list of explored cycles
        C.append(walk)
        # update count of available outgoing edges
        Tag = UpdateVisited(walk, Tag)
        # form new cycle by combing cyclePrime with new cycle
        cyclePrime = cyclePrime + walk[1:]
        # Cycle ← Cycle’
        cycle = cyclePrime
        # update list of nodes with availe outgoing edges
        L = [i for i in Tag if Tag[i]['out'] > 0]
        print('L', len(L))
        
    return cycle



def FormatCycle(cycle):
    
    print('->'.join(list(map(lambda x: str(x), cycle))))

def SolvePb(PbFile):
    # read graph from file
    Graph = ReadGraphFromFile(PbFile)    
    # find unbalanced nodes
    Unbalanced = FindUnbalancedNodes(Graph)
    # Add an edge to obtain a balanced node 
    Graph, edge = AddEdge(Graph, Unbalanced)
    # Find eularian cycle
    cycle = EulerianCycle(Graph)
    # break edge in cycle to obtain eularian path
    path = FromCycleToGraph(cycle, edge)
    # format for grader
    FormatCycle(path)
