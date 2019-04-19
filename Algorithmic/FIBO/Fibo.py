# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


# implement Fibonacci Numbers


def Fibo(n):
    '''
    (int) -> int
    Return the Fibonacci number of n
    '''
    
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        return Fibo(n-1) + Fibo(n-2)
    