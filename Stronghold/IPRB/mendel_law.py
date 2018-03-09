def mendel_law(k, m, n):
    '''
    (int, int, int) -> float
    Return the probability that 2 invididuals mating at random drawn from
    the pool of k dominant homozygotes, m heterozygotes and n recessive heterozygotes
    will produce an offspring with a dominant allele

    >>> mendel_law(2, 2, 2)
    0.78333
    '''

    # draw a tree to compute the probability of choosing each parent
    # and then adjusting each cross for the probability of having an offspring with a dominant allele

    T = k + m + n
    P = (k*(k-1))/(T*(T-1)) + (2*k*n)/(T*(T-1)) + (2*k*m)/(T*(T-1)) + (3/4 * (m*(m-1))/(T*(T-1))) + (m*n)/(T*(T-1))

    return P
