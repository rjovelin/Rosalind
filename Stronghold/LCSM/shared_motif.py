def common_substring(fasta):
    '''
    (file) -> str
    Return the longest common substring to all string sequences in the fasta file
    If multiple longest common substring exists, return any of them
    '''

    # make a dictionnary with the fasta sequences
    sequences = convert_fasta(fasta)

    # go through a sequence and find all common substrings
    # check all substrings adding a nucleotide from left to right
    names = []
    for seq in sequences:
        names.append(seq)

    
    commons = set()
    substring = ''

    for i in range(len(sequences[names[0]])):
        substring = sequences[names[0]][i]
        for j in range(i + 1, len(sequences[names[0]])):
            substring += sequences[names[0]][j] 
            total = 0
            for seq in sequences:
                if substring in sequences[seq]:
                    total += 1
            if total == len (sequences):
                commons.add(substring)
       

    sizes = []
    for item in commons:
        sizes.append(len(item))

    longest_commons = []
    for item in commons:
        if len(item) == max(sizes):
            longest_commons.append(item)
    if longest_commons != []:
        return longest_commons[0]

def convert_fasta(fasta):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    sequence ID as key and single string sequence as value
    '''
    
    genome = {}
    with open(fasta, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line
    return genome
