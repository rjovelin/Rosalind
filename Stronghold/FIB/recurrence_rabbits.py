# use a recursive function to calculte the number of rabbits after n months
def rabbits(n, k):
    '''
    (int, int) -> int
    Return the number of rabbit pairs present after n months if each
    reproduction-age pair of rabbits produce a litter of k rabbits pairs
    Precondition: we start with 1 pair, and it takes 1 month for rabbit to
    reach age of reproduction (ie. the number of offpsrings = N rabbits present 2 months before * k)
    '''
    
    # intitialise the recursive function
    if n == 1 or n == 2:
        return 1
    else:
        return rabbits(n-1, k) + (k * (rabbits(n-2, k)))

# use bottom-up dynamic programming to calculate the number of rabbits after n months    
def smart_rabbits(n, k):
    '''
    (int, int) -> int
    Return the number of rabbit pairs present after n months if each
    reproduction-age pair of rabbits produce a litter of k rabbits pairs
    Precondition: we start with 1 pair, and it takes 1 month for rabbit to
    reach age of reproduction (ie. the number of offpsrings = N rabbits present 2 months before * k)
    '''
    # create a dictionnary to store the values of n: 1 to n
    fibo = {}
    # initialize fibo with the 2 first values
    fibo[1] = 1
    fibo[2] = 1
    # compute the number of rabbits for each value before n
    for n in range(3, n+1):
        fibo[n] = fibo[n-1] + fibo[n-2] * k
    return fibo[n]