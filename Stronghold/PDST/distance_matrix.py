def distance_matrix(fasta_file):
    '''
    (file) -> str
    Return the p-distance matrix between sequences in the fasta_file
    Precondition: all sequences in the fasta file have the same length
    '''

    # convert the fasta file into a dictionary
    # change the sequence names to numbers
    seq = {}
    i = 1
    infile = open(fasta_file, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                name = i
                seq[name] = ''
                i += 1
            else:
                seq[name] += line.upper()
    infile.close()



    # make an ordered list with the sequence names
    names = [seq_name for seq_name in seq]
    names.sort()

    # compute p-distance between all pairs of sequence and print result
    for i in range(len(names)):
        for j in range(len(names)):
            p_dist = 0
            for k in range(len(seq[names[i]])):
                if seq[names[i]][k] != seq[names[j]][k]:
                    p_dist += 1
            print('{0:.5f}'.format(p_dist / len(seq[names[0]])), end = ' ')
        print('')
    
                
