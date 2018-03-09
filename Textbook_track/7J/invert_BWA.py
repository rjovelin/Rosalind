# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 22:38:01 2015

@author: Richard
"""



# Reconstruct a string from its Burrows-Wheeler transform
def invert_BWA(input_file):
    '''
    (file) -> str
    Given a string with a single "$" sign in input_file which is
    the Wheeler-Burrows Transform of a string Text, return 
    the string Text such that BWT(Text) = Transform.
    '''
    # open file for reading
    infile = open(input_file, 'r')
    # get Transform into a list
    transform = [x for x in infile.readline().rstrip()]
    # close file
    infile.close()

    # copy Transform to initiate the IBW matrix
    IBWA = transform[:]
       
    # sort last column of matrix gives 1st column
    IBWA.sort()
    
    # the number of iteration is length of Tranform -1
    n = len(transform) -1
    
    # iterate through cycles of adding elements of Trannsform to elements
    # in matrix and then sorting matrix
    while n != 0:
        # add elements of IBWA to elements in transform
        for i in range(len(transform)):
            IBWA[i] = transform[i] + IBWA[i]
    
        IBWA.sort()    
        n -= 1
    
    # find the word  starting with "$"
    for word in IBWA:
        if word.startswith('$'):
            text = word
            break
    
    # by definition, text ends with '$' so move '$' from start to end
    text = text[1:] + text[0]
    
    return text
            
    
    