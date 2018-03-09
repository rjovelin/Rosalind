def ts_tvs_ratio(file):
    '''
    (file) -> float
    Return the transition/transversion ratio between 2 sequences s1 and s2 in file
    Precondition: s1 and s2 have the same length


    >>> ts_tvs_ratio('ts_tvs_example.txt')
    1.21428571429
    '''



    # convert the fasta sequences into single string sequences and store them in a dictionnary
    ID_seq = {}

    myfile = open(file, 'r')
    for line in myfile:
        line = line.rstrip()
        if line == '':
            continue
        elif line.startswith('>'):
            ID_seq[line[1:]] = ""
            seq_name = line[1:]
        else:
            ID_seq[seq_name] += line
    myfile.close()

    # make a list with the 2 sequences
    sequences = []
    for values in ID_seq.values():
        sequences.append(values)

    # retrieve the sequences
    s1 = sequences[0]
    s2 = sequences[1]

    s1 = s1.upper()
    s2 = s2.upper()

    ts = 0
    tvs = 0

    assert len(s1) == len(s2), 's1 and s2 have unequal length'

    for i in range(len(s1)):
        if s1[i] != s2[i]:
            if s1[i] == 'A' and s2[i] == 'G':
                ts += 1
            elif s1[i] == 'G' and s2[i] == 'A':
                ts += 1
            elif s1[i] == 'C' and s2[i] == 'T':
                ts += 1
            elif s1[i] == 'T' and s2[i] == 'C':
                ts += 1
            else:
                tvs += 1

    return ts/tvs



            
