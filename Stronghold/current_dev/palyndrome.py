def find_palyndrome(S, min_L, max_L):
    '''
    (str, int, int) -> None
    Print the positions and length of all palyndromes in S that have length between min)L and max_L

    >>> find_palyndrome('TCAATGCATGCGGGTCTATATGCAT', 4, 12)
    4 6
    5 4
    6 6
    7 4
    17 4
    18 4
    20 6
    21 4
    '''

    from reverse_comp import reverse_complement

    positions = []

    while max_L >= min_L:
        for i in range(len(S)):
            if S[i: i + int(max_L / 2)] == reverse_complement(S[i + int(max_L / 2) : i + max_L]):
                pos = (i + 1, max_L)
                positions.append(pos)
        max_L -= 2

    for pos in positions:
        print(pos[0], pos[1])
            
            
            
            
            
    
