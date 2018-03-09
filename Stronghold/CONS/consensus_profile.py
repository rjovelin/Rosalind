def consensus_profile(fasta_file):
    '''
    (file) -> str
    Return the profile matrix corresponding to the sequences in the fasta_file
    along with the consensus sequence
    Precondition: all sequences in the fasta file have the same length
    '''

    # convert the fasta file into a dictionary
    seq = {}
    infile = open(fasta_file, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                name = line
                seq[name] = ''
            else:
                seq[name] += line.upper()
    infile.close()

    # make a dictionnary to store the nucleotide counts
    seq_names = [name for name in seq]
    first_seq = seq_names[0]

    nucleotides = {'A', 'T', 'C', 'G'}
    profile = {}
    for base in nucleotides:
        profile[base] = []
        for i in range(len(seq[first_seq])):
            profile[base].append(0)
    
    # update the number of each nucleotide at each position
    for name in seq:
        for i in range(len(seq[name])):
            if seq[name][i] == 'A':
                profile['A'][i] += 1
            elif seq[name][i] == 'C':
                profile['C'][i] += 1
            elif seq[name][i] == 'G':
                profile['G'][i] += 1
            elif seq[name][i] == 'T':
                profile['T'][i] += 1

    # compute consensus sequence
    consensus = ''
    for i in range(len(profile['A'])):
        sites = []
        sites.append(profile['A'][i])
        sites.append(profile['T'][i])
        sites.append(profile['C'][i])
        sites.append(profile['G'][i])
        consensus_base = max(sites)
        if profile['A'][i] == consensus_base:
            consensus += 'A'
        elif profile['T'][i] == consensus_base:
            consensus += 'T'
        elif profile['C'][i] == consensus_base:
            consensus += 'C'
        elif profile['G'][i] == consensus_base:
            consensus += 'G'

    # print results
    # make a list of nucleotides
    nucleo = [base for base in profile]
    nucleo.sort()
    
    print(consensus)
    for base in nucleo:
        print(base + ':', end = ' ')
        for item in profile[base]:
            print(item, end = ' ')
        print('')


            







    
    
            
