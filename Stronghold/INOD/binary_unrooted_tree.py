def binary_unrooted_tree(n):
    '''
    Return the number of internal nodes in an unrooted binary tree given the number of leaves n
    

    >>> binary_unrooted_tree(3)
    1
    >>> binary_unrooted_tree(4)
    2
    '''
    # in an unrooted binary tree, leaves have degree 1 and internal nodes have degree 3
    # the number of internal nodes is
    # N_leaves = N_internal_nodes + 2

    return n - 2
