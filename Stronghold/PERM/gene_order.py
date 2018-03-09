def gene_order(I):
    '''
    (int) -> tuple

    Return a tuple with the number of permutations of length I and the list
    containing all the permuations in tuples

    >>> gene_order(3)
    (6, [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)])
    '''


    # generate a list of  the range of I
    # use numbers instead of indices: add 1, and start at 1
    #import itertools
    # use the permutations methods
    # return the length of the permuatios and the permuations

    seq = list(range(1, I + 1))
    import itertools
    permutations = list(itertools.permutations(seq))
    return len(permutations), permutations



