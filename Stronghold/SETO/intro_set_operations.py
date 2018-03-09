# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:58:47 2015

@author: Richard
"""

def intro_set_operations(input_file, outputfile):
    '''
    (file) -> file
    Grab a positive interget n and (n≤20,000)
    and two subsets A and B of the set of {1,2,…,n} from the input file
    Save to file 6 sets: union of A and B, intersection of A and B, A−B,
    B−A, complement of A and complement of B
    (where set complements are taken with respect to {1,2,…,n})    
    '''
    # open file for reading
    infile = open(input_file, 'r')
    # grab n, A and B from file object
    n = int(infile.readline().rstrip())
    s1 = infile.readline().rstrip()
    s2 = infile.readline().rstrip()
    
    # get rid of the '{' and '}' symbols from each string
    s1 = s1[1:-1]
    s2 = s2[1:-1]
    
    # make lists with elements of each string
    s1 = s1.split(',')
    s2 = s2.split(',')
    
    # remove trailing space for each string element of both lists
    for i in range(len(s1)):
        s1[i] = s1[i].rstrip()
    for i in range(len(s2)):
        s2[i] = s2[i].rstrip()
    
    # make sets containing integers using both lists
    A = set(int(i) for i in s1)
    B = set(int(i) for i in s2)
    
    # make master set {1, 2, ...n}
    U = set(i for i in range(1, n+1))
    
    # compute the different sets and assign to lists
    AunB = list(A.union(B))
    AintersB = list(A.intersection(B))
    Aonly = list(A - B)
    Bonly = list(B - A)
    cA = list(U - A)
    cB = list(U - B)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # create list of lists to loop over
    all_sets = [AunB, AintersB, Aonly, Bonly, cA, cB]
    
    # loop over the lists and write to file to content of each list
    for set_item in all_sets:
        newfile.write('{')
        for item in set_item[:-1]:
            newfile.write(str(item) + ', ')
        newfile.write(str(set_item[-1]) + '}' + '\n')
    
    # close file after writing
    newfile.close()
    
    
    
    
    
    
    
    
    
    
    