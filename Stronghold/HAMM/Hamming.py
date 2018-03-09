def Hamming_D(s1, s2):
    '''
    (str, str) -> float
    Return the Hamming distance (the # of differences) between 2 sequences s1 and s2
    Precondition: s1 and s2 have the same length


    >>> Hamming_D('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT')
    7
    '''
    diffs = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diffs += 1
    return diffs
