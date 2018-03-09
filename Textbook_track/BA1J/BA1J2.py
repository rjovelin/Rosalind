# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:42:59 2016

@author: RJovelin
"""



def ReverseComplement(DNA):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    revcomp = ''
    for i in DNA.upper():
        revcomp = complement[i] + revcomp
    return revcomp.upper()


def FastPatternToNumber(Pattern):
    '''
    (str) -> int
    Take a string kmer and return the index of that kmer in the list of
    all possible kmers of length k ordered lexicographically
    '''
    
    # use the following relationship
    # PatternToNumber(Pattern) = 4 · PatternToNumber(Prefix(Pattern)) + SymbolToNumber(LastSymbol(Pattern))
    # where Prefix(Pattern) = Pattern minus last symbol
    # and SymBolToNumber of last symol of Pattern is defined by A = 0, C = 1, G = 2, T =3
    
    SymbolToNumber = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    if Pattern == '':
        return 0
    else:
        symbol = Pattern[-1]
        prefix = Pattern[:-1]
        return 4 * FastPatternToNumber(prefix) + SymbolToNumber[symbol]
    
    
   
def FastNumberToPattern(index, k):
    '''
    (str) -> int
    Take the index kmer in the list of all possible kmers of length k
    ordered lexicographically and the length k and return the corresponding kmer
    '''
    
    # use the following relationship
    # PatternToNumber(Pattern) = 4 · PatternToNumber(Prefix(Pattern)) + SymbolToNumber(LastSymbol(Pattern))
    # where Prefix(Pattern) = Pattern minus last symbol
    # and SymBolToNumber of last symol of Pattern is defined by A = 0, C = 1, G = 2, T =3
    
    NumberToSymbol = {0:'A', 1: 'C', 2: 'G', 3: 'T'}
    
    if k == 1:
        return NumberToSymbol[index]
    prefixIndex = index // 4
    remainder = index % 4
    symbol = NumberToSymbol[remainder]
    PrefixPattern = FastNumberToPattern(prefixIndex, k - 1)
    return PrefixPattern + symbol    



def HammingDistance(seq1, seq2):
    D = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            D += 1
    return D


def FrequentWordsMismatch(Genome, k, d):
    '''
    (str, int, int) -> set
    Return a dict of kmers with counts in Genome. Kmers in dict have a hamming
    distance > d (ie, each kmer counts its occurence or the occurence of
    similar kmers with < mismatches in Genome)
    
    >>>FrequentWordsMismatch('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
    {ATGT, ACAT}
    '''
    
    # create a frequency array
    frequency = []
    for i in range(4**k):
        frequency.append(0)
    # create a dict with pattern : counts
    counts = {}
    # initialize pattern : 0 values
    for i in range(len(frequency)):
        pattern = FastNumberToPattern(i, k)
        counts[pattern] = 0
    # loop through genome
    for i in range(len(Genome) - k + 1):
        # extract kmer
        kmer = Genome[i:i+k]
        for pattern in counts:
            if HammingDistance(kmer, pattern) <= d:
                # update pattern count
                counts[pattern] += 1
    revcomp = ReverseComplement(Genome)
    for i in range(len(revcomp) - k + 1):
        # extract kmer
        kmer = revcomp[i:i+k]
        for pattern in counts:
            if HammingDistance(kmer, pattern) <= d:
                # update pattern count
                counts[pattern] += 1
    
    # find highest count
    maxval = max([counts[pattern] for pattern in counts])
    # create set of most frequent patterns
    mostfrequent = set()
    for pattern in counts:
        if counts[pattern] == maxval:
            mostfrequent.add(pattern)
            
    return mostfrequent
    
    

    