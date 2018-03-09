# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:23:43 2015

@author: Richard
"""

# Find the minimum number of coins needed to make change
def optimum_change(input_file):
    '''
    (file) -> int
    Given an integer money and an array Coins of positive integers, 
    return the minimum number of coins with denominations Coins
    that changes money.
    example: money = 40, coins = [1,5,10,20,25,50]
    return 2    
    '''
    
    
    # use dynamic programming to find the optimum solution
    
    # precompute the number of coins to return for all values before money
    
    # at each value m, compare the minimum of mininumcoins(money -coin) + 1
    # for all coin
    
    # open file for reading
    infile = open(input_file)
    money = int(infile.readline().rstrip())
    coins = infile.readline().rstrip().split(',')
    # close file
    infile.close()
    # convert str to int
    for i in range(len(coins)):
        coins[i] = int(coins[i])
    
    
    # make a dictionnary to hold the values of minimumcoins for each value
    # of m up to money
    minNumcoins = {}
    
    # initialize dictionnary
    minNumcoins[0] = 0
    
    # loop over each value of m up to money
    # add 1 to money to include money in range
    for m in range(1, money + 1):
        # assign an very large value to minNumcoins[m] so it can be changed
        minNumcoins[m] = 1000000000
        # take the minimum of minNumcoins[m] and nimNumcoins(m - coin)
        for j in range(len(coins)):
            # check that m - coin >= 0
            if m >= coins[j]:
                # compare minNumcoins[m] and minNumcoins[m - coin]
                # eventually update minNumcoins[m]
                if minNumcoins[m - coins[j]] + 1 < minNumcoins[m]:
                    minNumcoins[m] = minNumcoins[m - coins[j]] + 1
    return minNumcoins[money]
                    
                    
        
    
    
    
    
    