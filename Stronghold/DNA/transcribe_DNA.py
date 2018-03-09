# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 23:14:48 2015

@author: Richard
"""

def transcribe_DNA(DNA):
    '''
    (str) -> str
    Return the RNA sequence transcribed from the DNA sequence
    '''
    
    DNA = DNA.upper()
    return DNA.replace('T', 'U')