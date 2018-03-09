# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:53:11 2015

@author: Richard
"""

import scipy.misc as msc
import decimal
import math


def maximum_matching(input_file):
    '''
    (file) -> int
    Given a RNA string S in input_file, return the total possible number
    of maximum matchings of basepair edges in the bonding graph of S
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # grab RNA seq
    S = ''
    for line in infile:
        if line != '\n' and not line.startswith('>'):
            line = line.rstrip()
            S += line
        
    S = S.upper()
    
    # count the number of A, U, C, G 
    A = S.count('A')
    U = S.count('U')
    C = S.count('C')
    G = S.count('G')
    
    decimal.Context(prec=10000)    
    
    
    # don't use approximation when computing factorial
    if A == U and C != G:
        if C > G:
            matching = msc.factorial(A, True) * decimal.Decimal(msc.factorial(C, True)/msc.factorial(C-G, True))
        elif C < G:
            matching = msc.factorial(A, True) * decimal.Decimal(msc.factorial(G, True)/msc.factorial(G-C, True))
    elif A != U and C != G:
        if A > U and C > G:
            matching = decimal.Decimal(msc.factorial(A, True)/msc.factorial(A-U, True)) * decimal.Decimal(msc.factorial(C, True)/msc.factorial(C-G, True))
        elif A > U and C < G:
            matching = decimal.Decimal(msc.factorial(A, True)/msc.factorial(A-U, True)) * decimal.Decimal(msc.factorial(G, True)/msc.factorial(G-C, True))
        elif A < U and C > G:
            matching = decimal.Decimal(msc.factorial(U, True)/msc.factorial(U-A, True)) * decimal.Decimal(msc.factorial(C, True)/msc.factorial(C-G, True))
        elif A < U and C < G:
            matching = decimal.Decimal(msc.factorial(U, True)/msc.factorial(U-A, True)) * decimal.Decimal(msc.factorial(G, True)/msc.factorial(G-C, True))
    elif A == U and C == G:
        matching = msc.factorial(A, True) * msc.factorial(C, True)
    elif A != U and C == G:
        if A > U:
            matching = decimal.Decimal(msc.factorial(A, True)/msc.factorial(A-U, True)) * msc.factorial(C, True)
        elif A < U:
            matching = decimal.Decimal(msc.factorial(U, True)/msc.factorial(U-A, True)) * msc.factorial(C, True)
    
    
    
    
    print(A, U, C, G)    
    print(matching)
    print(decimal.Decimal(matching))    