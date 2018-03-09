def overlap_graph(k, fasta):
    '''
    (int, file) -> list
    Return the list of adjacent strings in fasta that are connected by substring of length k
    respectively as their suffix and prefix
    '''


    genome = {}
    with open(fasta, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line

    # make a list of keys
    sequences = []
    for sequence in genome:
        sequences.append(sequence)

    # for each sequence in the genome compare its suffix with the preffix of the next sequence
    # ignore edges that connect the string to itself

    directed_edges = []
    for name in sequences:
        for sequence in genome:
            if name != sequence:
                if genome[name][-k:] == genome[sequence][:k]:
                    pair = (name, sequence)
                    directed_edges.append(pair)

    return directed_edges
                    
                
