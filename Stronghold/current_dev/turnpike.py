# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 11:00:10 2015

@author: Richard
"""
# import necessary modules
import math
from copy import deepcopy

# function to resolve the tunrpike problem
def turnpike(input_file):
    '''
    (file) -> list
    Given a multiset L containing C(n, 2) positive integers for some
    positive integer n, return a set X containing n nonnegative integers
    such that ΔX=L. ΔX is defined as the collection of all positive
    differences between elements of X
    '''
    
    
    # open file for reading
    infile = open(input_file, 'r')
    # get list of positive differences
    L = infile.readline().rstrip().split()
    for i in range(len(L)):
        L[i] = int(L[i])
    # close file
    infile.close()
    
    # get the number of element in set X
    # resolve quadratic equation 
    m = len(L)
    n = int((1+math.sqrt((1 -4*(-2*m)))/2))
    
    # make a list of n negatif elements that will be changed
    # the final list must have only positive elements
    X = [-1] * n
    
    # first element in list is 0
    X[0] = 0
    # last element is the highest number in L
    X[-1] = max(L)
    # remove highest number from L
    L.remove(max(L))
    
    # now the problem consists in placing points in list X such that
    # distances between points are elements of list L
    
    # use pointers to indicate the location in X of the next element
    # we know the position of the first and last elements so the next element
    # to place will be at 2nd position or 2nd to last
    i = 1
    j = -2
    
    # until no -1 in X:
    # take largest number, check if Xn - largest in L
    # can assign lasgest to X[j] or Xn-largest to X[i]: check if 
    # try one, make deepcopy of L and X, see if a solution exists
        # stop if solution doesn't exist and 
    
    while -1 in X and len(L) != 0:
        largest = max(L)
        
        if X[-1]-largest in L:
            X1 = deepcopy(X)
            L1 = deepcopy(L)
            
            
                   
            
            
            highest_position = True
        if largest - X[-1] in L:
            lowest_position = True
        # test solution with largest value
        while highest_position:
            k = j
            X1[k] = largest
            # remove Xn-largest and largest - X0
            L1.remove(X1[-1] - largest)
            L1.remove(largest - X1[0])
            # reassign k
            k -= 1
            newlargest = max(L1)
        if X[-1] - largest in L:
            X2[j] = largest


# take the largest distance:
    # 2 possible solution: X[-2] = largest or X[1] = X[-1] - largest
    # check that X[-1] - largest is in L and other condition
    # assign largest to X[-2]
    # remove x[-1] - largest and largest - x[0]
    
    # assign X[-1] - largest to X[1]
    # remove x[-1] - largest and largest - x[0]
    
    # evaluate both possibilities.
    # for each possibility, evaluate the 2 possibilities with the new largest
    # if both possibilities are correct, then backtrack and choose the X[-1] in previous step
    # if one possibility is notcorrect, then back tracj and choose tthe assignmenent accordingly in the previous step


The next step is not obvious. Since 7 is the largest value in D , either x 4=7 or x 2=3. If
x 4=7, then the distances x 6−7=3, x 5−7=1, must also be present in D . A quick check shows that
indeed they are. On the other hand, if we set x 2=3, then 3−x 1=3 and x 5−3=5 must be present in
D. These distances are also in D , so we have no guidance on which choice to make. Thus we try
one, and see if it leads to a solution. If it turns out that it doesn’t, we can come back and try the
other. Trying the first choice, we set x 4=7, which leaves



At this point, we have x 1=0, x 4=7, x 5=8, and x 6=10. Now the largest distance is 6, so either
x 3=6 or x 2=4. But if x 3=6, then x 4−x 3=1, which is impossible since 1 is not in D . On the other
hand, if x 2=4, then x 2−x 0=4, and x 5−x 2=4. This is also impossible, since 4 only appears once in
D . Thus, this line of reasoning leaves no solution, so we backtrack.
Since x 4=7 failed to produce a solution, we try x 2=3. If this also fails, then we give up and
report no solution. We now have
x 1=0 x 2=3 x 5=8 x 6=10
D = { 1, 2, 2, 3, 3, 4, 5, 5, 6 }.
Once again, we have to choose between x 4=6 and x 3=4. x 3=4 is impossible because D only
has one occurrence of 4, and two would be implied by this choice. x 4=6 is possible, so we obtain
x 1=0 x 2=3 x 4=6 x 5=8 x 6=10
D = { 1, 2, 3, 5, 5}.
The only remaining choice is to assign x 3=5; this works because it leaves D empty, and so we
have a solution.
x 1=0 x 2=3 x 3=5 x 4=6 x 5=8 x 6=10
D = { }.
Fig. 0.1 shows a decision tree representing the actions taken to arrive at the solution.
Instead of labelling the branches, we’ve placed the label in the branches’ destination node. A
node with an asterisk indicates that the points chosen are inconsistent with the given distances;
nodes with two asterisks have only impossible nodes as children, and thus represent an incorrect
path.