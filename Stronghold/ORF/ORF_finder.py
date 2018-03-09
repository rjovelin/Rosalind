def CDS_finder(dna):
    '''
    (str) -> list
    Return a list of ORFs found in the direct sense of dna

    >>> ORF_finder('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
    {'M', 'MGMTPRLGLESLLE', 'MTPRLGLESLLE'}
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

    dna = dna.upper()
    rna = dna.replace('T', 'U')
    CDS_list = set()

    # walk through rna
    for i in range(len(rna)):
        CDS = ''
        j = i
        # make sure the key is in the dictionnary
        if rna[j: j+3] in genetic_code and genetic_code[rna[j: j+3]] == 'M':
            # start translating from that codon until a stop codon is reached
            for j in range(j, len(rna), 3):
                if rna[j: j+3] in genetic_code and genetic_code[rna[j: j+3]] == '*':
                    CDS = CDS + genetic_code[rna[j: j+3]]
                    break
                elif rna[j: j+3] in genetic_code and genetic_code[rna[j: j+3]] != '*':
                    CDS = CDS + genetic_code[rna[j: j+3]]
            # keep only CDS starting with M and ending with *
            # but do not include the stop codon in the list of CDS
            if CDS.endswith('*'):
                CDS_list.add(CDS[:-1])

    return CDS_list

def six_frames_CDS_finder(dna):
    '''
    (str) -> list
    Return a list of ORFs found in the direct sense of dna and in its reversed complement

    >>> six_frames_CDS_finder('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
    {'MLLGSFRLIPKETLIQVAGSSPCNLS, 'M', 'MGMTPRLGLESLLE', 'MTPRLGLESLLE'}
    '''

    direct_sense = CDS_finder(dna)

    import reverse_comp

    reverse_dna = reverse_comp.reverse_complement(dna)
    anti_sense = CDS_finder(reverse_dna)

    six_frames = direct_sense.union(anti_sense)

    return six_frames
    

