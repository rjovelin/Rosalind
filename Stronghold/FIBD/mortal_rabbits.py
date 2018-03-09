# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 20:25:53 2015

@author: Richard
"""

# use bottom-up dynamic programming to calculate the number of rabbits after n months    
def mortal_rabbits(input_file):
    '''
    (file) -> int
    Given positive intergers n and m in input file, return the total number
    of pairs of rabbits that will remain after the n-th month if all rabbits
    live for m months. Precondition: we start with 1 pair, and it takes 1
    month for rabbit to reach age of reproduction 
    '''
    
    # create a dictionnary to store the values of n: 1 to n
    fibo = {}
    # initialize fibo with the 2 first values
    
    # initialize at fibo[0] because fibo[i - (m+1)] returns key 0
    fibo[0] = 1
    fibo[1] = 1
    # compute the number of rabbits for each value before n
    for i in range(2, n):
        # the recurrency relationaship changes dpending on the value of i
        if i < m:
            fibo[i] = fibo[i-1] + fibo[i-2]
        elif i == m:
            fibo[i] = fibo[i-1] + fibo[i-2] - 1
        else:
            fibo[i] = fibo[i-1] + fibo[i-2] - fibo[i - (m+1)]
    # return the value at the nth mobth = n-1
    return fibo[n-1]
 
 
 
 
 
 
