# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:40:08 2016

@author: RJovelin
"""


# use this function to get the position of pattern in the list of all possible patterns
# ordered lexicographically
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
    