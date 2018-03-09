# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:39:41 2016

@author: RJovelin
"""



# use this function to find the pattern of lenght k given its index in the 
# list of all possible kmers ordered lexicographically   
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