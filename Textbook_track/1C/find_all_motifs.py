def find_all_motifs(a, b):
    '''
    (str, str) -> list
    Return all the indices of string a in string b into a list

    >>>find_all_motifs('ATAT', 'GATATATGCATATACTT')
    [1, 3, 9]
    '''
    i = 0
    pos = []
    found = True
    while found:
        i = b.find(a, i)
        if i >=0:
            pos.append(i)
            i += 1
        else:
            found = False

    return pos
