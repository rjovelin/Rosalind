def spliced_motif(fasta):
    '''
    Get the sequences S and t from the fasta file and return one collection of 
    positions of S in which symbols of t appear in S

    >>> spliced_motif('ACGTACGTGACG', 'GTA')
    3 8 10

    '''


    sequences = {}
    with open(fasta, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                sequences[line[1:]] = ""
                seq_name = line[1:]
            else:
                sequences[seq_name] += line
    
    dna = []
    for seq in sequences:
        dna.append(sequences[seq])
    if len(dna[0]) <= len(dna[1]):
        t = dna[0]
        S = dna[1]
    else:
        S = dna[0]
        t = dna[1]
        
    positions = []
    j = 0
    for i in range(len(t)):
        positions.append(S.index(t[i], j))
        j = S.index(t[i], j) + 1

    for item in positions:
        print (item + 1, end = ' ')
                         
