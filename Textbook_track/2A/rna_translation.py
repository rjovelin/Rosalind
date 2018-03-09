def rna_translation(rna):
    '''
    (str) -> str
    Return a protein sequence translation of the string rna given a genetic code

    >>> rna_translation('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    MAMAPRTEINSTRING
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


    protein = ''
    for i in range(0, len(rna), 3):
        codon = rna[i: i+3]
        protein  += genetic_code[codon]

    return protein
