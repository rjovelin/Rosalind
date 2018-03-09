def reverse_complement(dna):
    '''
    (str) -> str
    Return the reverse complement of a dna string

    >>> reverse_complement('ATGTC')
    GACAT
    '''

    dna2 = dna.upper()
    compl_dna = ''

    for base in dna2:
        if base == 'A':
            compl_dna += 'T'
        elif base == 'T':
            compl_dna += 'A'
        elif base == 'C':
            compl_dna += 'G'
        elif base == 'G':
            compl_dna += 'C'
            
    rev_compl = ''
    for base in compl_dna:
        rev_compl = base + rev_compl

    return rev_compl
        
