def make_tree(file):
    '''
    Read the file and get the number of nodes in the first line and the pairs of connected nodes in each subsequent line.
    Return the number of edges necessary to produce a tree with the given number of nodes
    '''

    # read the file, store the list of connections in a larger list
    # Total_edges = Total_nodes -1
    # get the number of total nodes
    # count the number of current edges
    # return the number of additional edges

    

    myfile = open(file, 'r')
    Total_nodes = int(myfile.readline().rstrip())
    connections = []
    for line in myfile:
        line = line.rstrip()
        line = line.split()
        connections.append(line)

    current_edges = len(connections)
    Total_edges = Total_nodes - 1

    missing_edges = Total_edges - current_edges

    myfile.close()


    return missing_edges
