# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 14:03:12 2015

@author: Richard
"""



def most_frequent_element(L):
    '''
    (list) -> element of list
    Return the most frequent element from list L
    >>> most_frequent_element([1, 2, 5, 6, 1, 1, 5, 2, 1, 8, 9, 10])
    1
    >>> most_frequent_element(['a', 'a', 'b', 1, 1, 'a', 'c', [1, 3]])
    'a'
    '''
    # set up the maximum frequency
    frequency = 0
    
    # go through the list, count each element, update the frequency
    # and update the most frequent element
    for item in L:
        if L.count(item) > frequency:
            most_frequent = item
            frequency = L.count(item)
            
    return most_frequent

def spectral_convolution(spectrum):
    '''
    (list) -> list
    Compute the spectral convolution of spectrum and return a list of elements
    in the convolution
    >>> spectral_convolution([0, 137, 186, 323])
    [137, 186, 49, 323, 186, 137]
    '''
    
    # create a dictionnary to hold the elements on spectum convolution
    convolution = {}
    # initialize i
    i = 1
    # compute the elements of spectrum convolution as absolute value
    # of each 2 elements that are different: ie: compute a matrix of differences
    while i != len(spectrum):
        # initialize value of i
        convolution[i] = []
        for k in range(i):
            convolution[i].append(abs(spectrum[i] - spectrum[k]))
        i += 1
    # create  a list to store the elements of the spectral convolution
    spectral = []
    # add all elements to list
    for i in convolution:
        for item in convolution[i]:
            # add only item > 0
            if item > 0:
                spectral.append(item)
    
    return spectral


def multiplicity_spectral_convolution(spectrum):
    '''
    (list) -> list
    Return a list of elements from the convolution of spectrum ordered in
    decreasing order of their multiplicity. Return each element the number
    of times it appears in the convolution
    >>> multiplicity_spectral_convolution([0, 137, 186, 323])
    [137 137 186 186 323 49]
    '''
    
    # compute the convolution of spectrum
    convolution = spectral_convolution(spectrum)
    
    # make a dictionnary with element : count pairs
    multiplicity = {}
    for item in convolution:
        multiplicity[item] = convolution.count(item)
    
    # create a list of tuples (count, item) for item in multiplicity
    counts = [(count, item) for item, count in multiplicity.items()]
    # reverse sort the list
    counts.sort()
    counts.reverse()
    
    # create a list to store all elements by decreasing order
    ordered = []
    for pair in counts:
        j = pair[0]
        while j != 0:
            ordered.append(pair[1])
            j -= 1
    
    return ordered
        
    

if __name__ == '__main__':
    filename = input('enter filename: ')
    infile = open(filename, 'r')
    spectrum = infile.readline().rstrip().split()
    infile.close()
    for i in range(len(spectrum)):
        spectrum[i] = int(spectrum[i])
    multiplicity = multiplicity_spectral_convolution(spectrum)
    newfile = open('spectral_multiplicity.txt', 'w')
    for item in multiplicity:
        newfile.write(str(item) + ' ')
    newfile.close()
    




# code for printing a matrix
#i = 0
#while i != len(spectrum):
#        print(i, end  = ' ')
#        for k in range(i):
#            print(abs(spectrum[i] - spectrum[k]), end = ' ')
#        print('\n')
#        i += 1