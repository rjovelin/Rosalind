def reverse_translate(peptide):
    '''
    (str) -> list
    Return a list of dna strings that all encode for peptide
    >>> reverse_translate('M')
    ['ATG']
    >>> reverse_translate('*')
    ['TAA', 'TAG', 'TGA']

    
    '''
    
    genetic_code = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
                    'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
                    'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                    'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
                    'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
                    'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                    'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                    'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                    'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                    'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                    'UAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                    'UAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                    'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
                    'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                    'UGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                    'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}


    # make a reverse genetic code to look for aa
    aa_to_codon = {}
    for codon, aa in genetic_code.items():
        if aa in aa_to_codon:
            aa_to_codon[aa].append(codon)
        else:
            aa_to_codon[aa] = [codon]
    
    # make a list of list with all possible codon for each aa in peptide
    possible_codons = []
    for aa in peptide:
        possible_codons.append(aa_to_codon[aa])
    
    # make a list of all substrings corresponding to the peptide
    # store the result of the itertools.product in a temporary list
        
    import itertools
    storage = list(itertools.product(*possible_codons))

    # make a list of substrings by adding all the codons contained in the list of tuples storage 
    # convert each rna substrings into dna
    amino = ''
    all_substrings = []
    for i in range(len(storage)):
        for j in range(len(storage[i])):
                       amino += storage[i][j]
        amino = amino.replace('U', 'T')
        all_substrings.append(amino)
        amino = ''

    return all_substrings


def all_seq_to_peptide(dna, peptide):
    '''
    (str, str) -> list
    Return the list of all possible substrings of dna in the sense and antisense that encode for peptide

    >>> all_seq_to_peptide('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')
    ['ATGGCC', 'GGCCAT', 'ATGGCC']
    '''

    subs_sense = reverse_translate(peptide)
    

    # count each substring in dna
    sense = []
    for sub in subs_sense:
        count = dna.count(sub)
        while count > 0:
            sense.append(sub)
            count -= 1

    # get the the substrings in the antisense
    import reverse_complement

    subs_antisense = []
    for sub in subs_sense:
        rev_sub = reverse_complement.reverse_complement(sub)
        subs_antisense.append(rev_sub)

    rev = reverse_complement.reverse_complement(dna)
    antisense = []
    for sub in subs_antisense:
        count = dna.count(sub)
        while count > 0:
            antisense.append(sub)
            count -= 1
    
    sense.extend(antisense)

    return sense



    
