# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:35:46 2015

@author: Richard
"""

from math import sqrt

# An atom is defined by its three coordinates (x, y, z) in 3D space.
# Determine the distances between any atoms
def atomic_distance(input_file):
    '''
    (file) -> list
    Given five atom's coordinates in space in input file, return the list of
    distances between all pairs of distinct atoms in any order with a
    decimal precision of 2.

Sample Dataset

17.426 -32.764 65.278
2.109 -38.295 41.517
-7.758 -42.568 33.470
-5.866 -3.238 13.986
-4.720 -38.377 -4.862
Sample Output

28.806 41.738 63.601 73.767 13.430 45.283 46.879 43.932 38.680 39.891    
    
    '''
    # The distance between two points in space is given by Pythagora's theorem:
    # d = square_root((x2 - x1)**2 + (y2 -y1)**2 + (z2-z1)**2)
    
    # open file for reading
    infile = open(input_file, 'r')
    # make a list to store the atoms' coordinates
    atoms = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            for i in range(len(line)):
                line[i] = float(line[i])
            atoms.append(line)
    
    # create list to store the distances
    distances = []
    
    # loop over the list of atoms, compute d between any 2 different atoms
    for i in range(len(atoms)-1):
        for j in range(i+1, len(atoms)):
            # assign atoms to be compared
            atom1 = atoms[i]
            atom2 = atoms[j]
            # compute the distance between atom1 and atom2
            d = sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)
            distances.append(d)
    
    # print distances in expected format
    for diff in distances:
        print(round(diff, 3), end = ' ')
    
