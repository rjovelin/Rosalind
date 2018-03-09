def rna_translate(rna):
    '''
    (str) -> str
    Translate a RNA sequence into a protein sequence according to the standard genetic code

    >>> rna_translate('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    MAMAPRTEINSTRING
    >>> rna_translate('AUGUACUAA')
    MY*
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

    RNA = rna.upper()
    protein = ''

    for i in range(0, len(RNA), 3):
        codon = RNA[i:i+3]
        protein += genetic_code[codon]

    return protein
        
    



    
