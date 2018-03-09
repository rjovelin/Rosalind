def translate_RNA_splice(file):
    '''
    (file) -> str
    File contains a dna sequence and introns in FASTA format
    Return the protein translation of the dna sequence after introns have been spliced


    >>> translate_RNA_splice('splice_example.txt')
    MVYIADKQHVASREAYGHMFKVCA*
    '''


    # convert the FASTA sequences into a single string
    # store the single string sequences and their ID in a dictionnary
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

    # make a list with the single string sequences
    single_seq = []
    for key in ID_seq:
        single_seq.append(ID_seq[key])

    # make a list with the length of single_seq
    seq_length = []
    for i in single_seq:
        seq_length.append(len(i))
    # get the index of the longest single_seq
    longest = seq_length.index(max(seq_length))
    # add each single_seq to the list of sequences starting with the longest
    sequences = []
    sequences.append(single_seq[longest])
    single_seq.remove(single_seq[longest])
    for i in single_seq:
        sequences.append(i)

 # if the sequence is already a single string and not fasta the vode below is enough       
##    sequences = []
##    for line in myfile:
##        line = line.rstrip()
##        if not line.startswith('>'):
##            sequences.append(line)


    # verify that  the first sequence in the list is the dna sequence
    for i in range(1, len(sequences)):
        assert len(sequences[0]) >= len(sequences[i]), 'dna is not the first sequence'

    # make a list of tuple that includes the start and end positions of each intron
    dna = sequences[0]
    intron_positions = []
    for i in range(1, len(sequences)):
        intron = sequences[i]
        intron_start = dna.index(intron)
        intron_end = intron_start + len(intron) -1
        intron_positions.append((intron_start, intron_end))

    # sort the list to get the introns in their order of appearance in dna
    intron_positions.sort()
    exons = []
    
    # add the first exon to the spliced_dna
    exon1 = dna[0:intron_positions[0][0]]
    spliced_dna = exon1
    # add the other exons using the intron_positions as coordinates
    for i in range(1, len(intron_positions)):
        exon = dna[intron_positions[i-1][1]+1:intron_positions[i][0]]
        spliced_dna += exon
    # add the last exon
    spliced_dna = spliced_dna + dna[intron_positions[-1][1]+1:]

    rna = spliced_dna.upper().replace('T', 'U')

    import rna_translate

    
    protein = rna_translate.rna_translate(rna)

    return protein

