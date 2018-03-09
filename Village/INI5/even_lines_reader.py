# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:04:31 2015

@author: Richard
"""

def read_even_numbered_lines(input_file, output_file):
    '''
    (file) -> file
    Read the text (<= 1000 lines) in input file
    and return a file containing all the even-numbered lines
    from the original file. Assume 1-based numbering of lines
    '''
    
    # open files
    infile = open(input_file, 'r')
    newfile = open(output_file, 'w')
    i = 1    
    for line in infile:
        if i % 2 == 0:
            newfile.write(line)
        i += 1
    
    # close files
    infile.close()
    newfile.close()
        
    
    
    
    