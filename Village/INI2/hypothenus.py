# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 13:44:59 2015

@author: Richard
"""

def hypothenus(input_file):
    '''
    (file) -> int
    Get the 2 positive integers a and b from the input file
    and return the integer corresponding to the square of the hypotenuse
    of the right triangle whose legs have lengths a and b
    '''
    
    # open file for reading
    infile = open(input_file, 'r')
    line = infile.readline().rstrip().split()
    a = int(line[0])
    b = int(line[1])
    
    # close file
    infile.close()
    
    return a**2 + b**2
    
    
    
    