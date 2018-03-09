def clump_finding(S, k, L, t):
    '''
    (str, int, int, int) -> set
    Return a set of containing all k-mer substrings that form a clump
    of at least t times in the L interval within string S
    
    >>> S = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
    >>> clump_finding(S, 5, 75, 4)
    {'CGACA', 'GAAGA', 'AATGT'}

    >>> S2 = 'GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC'
    >>> clump_finding(S2, 4, 25, 3)
    {'TGCA'}

    '''

    K_mers = set()
        
    for i in range(len(S)-L):
        interval = S[i: i + L]
        for i in range(len(interval)-k):
            motif = interval[i:i+k]
            occurence = interval.count(motif)
            if occurence >= t:
                K_mers.add(motif)
                
    return K_mers
                    
