def protein_weight(protein):
    '''
    (str) -> float
    Return the protein weight as the sum of the weight of each amino acid

    >>> protein_weight('SKADYEK')
    821.392
    '''

    # generate a dictionnary of amino acid weight from the monoisotopic table
    mass_table = open('monoisotopic.txt', 'r')
    aa = []
    weight = []
    for line in mass_table:
        line = line.split()
        aa.append(line[0])
        weight.append(line[1])
    mass_table.close

    mass_aa = {}
    for i in range(len(aa)):
        mass_aa[aa[i]] = float(weight[i])

    mass_protein = 0
    for amino_acid in protein:
        mass_protein += mass_aa[amino_acid]

    return mass_protein


    











 
