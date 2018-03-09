def find_all_motifs(input_file):
    '''
    (file) -> None
    Return all the positions of string a in string b into a list,
    using index-1 based numbering

    >>>find_all_motifs('ATAT', 'GATATATGCATATACTT')
    [1, 3, 9]
    '''

    # open file for reading
    infile = open(input_file, 'r')
    seq1 = infile.readline().rstrip()
    seq2 = infile.readline().rstrip()

    if len(seq1) > len(seq2):
        b, a = seq1, seq2
    else:
        b, a = seq2, seq1

    i = 0
    pos = []
    found = True
    while found:
        i = b.find(a, i)
        if i >=0:
            pos.append(i)
            i += 1
        else:
            found = False

    for item in pos:
        print(item+1, end = ' ')
