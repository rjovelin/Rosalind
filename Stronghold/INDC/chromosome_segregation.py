# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 01:10:13 2015

@author: Richard
"""


from scipy import misc
import decimal
decimal.Context(prec = 100)



def independent_segregation(input_file):
    '''
    (file) -> list
    Given a positive integer n ≤ 50 return an array A of length 2n in which
    A[k] represents the common logarithm of the probability that two
    diploid siblings share at least k of their 2n chromosomes
    Precondition: assume no recombination    
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    # get n
    n = int(infile.readline().rstrip())
    # close file
    infile.close()
    
    # the probability of k successes in n trials is given by the binomial
    # distribution P(k) = C(n,k) * p^k * q^(n-k)
    
    p = 0.5
    N = 2*n
    
    for k in range(1, N+1):
        P = 0
        for i in range(k, N+1):
           P += misc.comb(N, i) * p**i * (1-p)**(N-i)
        Q = round(math.log10(P), 3)
        if Q == 0:
            print('0.000', end = ' ' )
        else:
            Q = str(Q)
            print(str(Q[:Q.index('.')+4]), end = ' ')
        
        
 
