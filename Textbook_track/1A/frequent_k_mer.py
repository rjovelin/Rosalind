def frequent_k_mer(seq, k):
    '''
    (str, int) -> str
    Return the number of k-mers substring of seq of length k such that k is maximized

    >>> frequent_k_mer('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
    CATG GCAT
    '''
    # make a list of all k-mers of length k substring of seq
    k_mer_all = []
    for i in range(0, len(seq)-k):
        k_mer = seq[i: i+k]
        k_mer_all.append(k_mer)

    # count each k-mer and store the count in a dictionnary
    k_mer_to_count = {}
    for item in k_mer_all:
        if item in k_mer_to_count:
            k_mer_to_count[item] += 1
        else:
            k_mer_to_count[item] = 1
            
    # reverse the dictionnary to access the list of k-mers by their count
    count_to_k_mer = {}
    for key, value in k_mer_to_count.items():
        if value in count_to_k_mer:
            count_to_k_mer[value].append(key)
        else:
            count_to_k_mer[value] = [key]

    # get the highest number of occurence
    longest_k = max(count_to_k_mer)

    #return the list of k-mer with the highest occurence
    return count_to_k_mer[longest_k]
