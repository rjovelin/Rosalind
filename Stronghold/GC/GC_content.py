def GC_calculate(s):
    '''
    (str) -> float
    Return the GC content of sequence s
    '''
    S = s.upper()
    G, C = S.count('G'), S.count('C')
    # multiply by 100 to get GC %
    return (G + C)/len(s) * 100 

def GC_content(file_name):

    '''
    (file) -> file
    Print in a new_file the ID of the sequence in file_name that has the highest GC content
    along the value of the GC content
    '''
       
    ID_seq = {}

    with open(file_name, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                ID_seq[line[1:]] = ""
                seq_name = line[1:]
            else:
                ID_seq[seq_name] += line

    names = []
    GCs = []
    for key, value in ID_seq.items():
        names.append(key)
        GCs.append(GC_calculate(value)) 
    
                    
    result_file = open('GC_result.txt', 'w')

    i = GCs.index(max(GCs))
    result_file.write(names[i] + '\n')
    result_file.write(str(max(GCs)))
    result_file.close()

    
