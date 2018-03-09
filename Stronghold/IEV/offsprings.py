def expected_offsprings(offspring):
    '''
    Return the number of expected offsprings with the dominant phenotype in the next generation
    given that each couple has exactly 2 offsprings
    '''


    myfile = open(offspring, 'r')
    expected = 0
    line = myfile.readline()
    line = line.rstrip()
    line = line.split()
    myfile.close()

    p = [1, 1, 1, 3/4, 1/2, 0]

    total = 0
    for i in range(6):
        total += (2 * int(line[i]) * p[i])

    return total
    
    
    
