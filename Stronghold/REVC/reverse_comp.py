def reverse_complement(dna):
    '''
    (str) -> (str

    >>> reverse_complement('atcg')
    'CGAT'
    '''

    dna2 = dna.upper()
    dna_comp = ''
    for i in dna2:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'

    reverse_comp_dna = ''
    for i in reversed(dna_comp):
        reverse_comp_dna += i

    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
        
    return reverse_comp_dna


def rev_compl(dna):

    '''
    (str) -> str
    return the reverse complement of string dna

    >>> rev_compl('atcg')
    'cgat'
    >>> rev_compl('TTCGAT')
    'ATCGAA'
    '''

    dna2 = dna.upper()
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    reverse_comp_L = []
    for i in reversed(dna2):
        reverse_comp_L.append(complement[i])
    reverse_comp_dna = ''.join(reverse_comp_L)

    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()

    return reverse_comp_dna
        
        

	
