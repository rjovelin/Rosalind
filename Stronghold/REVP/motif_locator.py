def locate_motif(seq, sub):

    '''
    (str, str) -> int
    Return all the position of substring sub in the string seq
    Precondition: length of seq >= length of sub
    
    >>> locate_motif('GATATATGCATATACTT', 'ATAT')
    2 4 10
    '''

    # if len(sub) > len(seq):
        # print ('sub is larger than seq')
    indices = []
    i = 0
    loc = 0
    
    if len(sub) > len(seq):
        print('sub should be smaller than seq')
    else:
        if sub not in seq:
            return -1
            print('sub is not in seq')
        else:
            # get the list of indices of sub in seq
            while loc != -1:
                loc = seq.find(sub, i)
                if loc != -1:
                    indices.append(loc)
                    i = loc + 1

    # get the positions of sub in seq relative to the 1st nucleotide
    positions = []      
    for item in indices:
        positions.append(item+1)

    for item in positions:
        print(item, end = ' ')
