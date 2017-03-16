# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.




import itertools
import sys



# Store as global variable since its values must be constant.
keys = 'ACGT'


# Converts single digit to single symbol.
def NumberToSymbol(Value):
    if Value == 0:
        return 'A'
    elif Value == 1:
        return 'C'
    elif Value == 2:
        return 'G'
    elif Value == 3:
        return 'T'
    else:
        return 'Error'

# Final defnition of function converting numbers to patterns.
# Uses recursion.
def NumberToPattern(Number, k):
    if k == 1:
        return NumberToSymbol(Number)

    Prefix = Number // 4
    Remainder = Number % 4

    Pattern = NumberToSymbol(Remainder)

    PrefixPattern = NumberToPattern(Prefix, k-1)

    return PrefixPattern + Pattern
            




'''
We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi. For example, CGAAT and CGGAC have two mismatches.
The number of mismatches between strings p and q is called the Hamming distance between these strings and is denoted
HammingDistance(p, q).

Hamming Distance Problem: Compute the Hamming distance between two strings.
     Input: Two strings of equal length.
     Output: The Hamming distance between these strings.
'''

def HammingDistance(p, q):
    if len(p) != len(q):
        return 'Invalid input'

    hamming_dist = 0 # Initialize to zero

    for i in range(len(p)):
        if p[i] != q[i]:
            hamming_dist += 1
    #print hamming_dist
    return hamming_dist





##    Instead of generating all possible motifs in DNA we instead generate k-mers in each
##    string in DNA. We then compute the minimum Hamming Distance between any k-mer 
##    and Text.
##
##    We perform this computation for each string in DNA then add together the resulting
##    minimums.

##    INPUT: string Pattern, list of strings DNA
##    OUTPUT: d(DNA, Pattern)
def DistanceBetweenPatternAndStrings(DNA, Pattern):
    k = len(Pattern)
    distance = 0

    for Text in DNA:
        HammingDistanceMax = sys.maxsize

        for i in range(len(Text) - k + 1):
            test_string = Text[i:i+k]
                
            condition = HammingDistance(Pattern, test_string)
            if HammingDistanceMax > condition:
                HammingDistanceMax = condition
        distance = distance + HammingDistanceMax

    return distance




##     INPUT: Int k, list of strings DNA.
##     OUTPUT: k-mer Pattern that minimizes d(Pattern, DNA) among all k-mers Pattern.
##             (Return the first one if there are multiple.)
##    
def MedianString(DNA, k):
    distance = sys.maxsize  # Initialize distance to "infinity"
    Median = ''

    for val in range(4**k):
        Pattern = NumberToPattern(val, k)
        test = DistanceBetweenPatternAndStrings(DNA, Pattern)

        if distance > test:
            distance = test
            Median = Pattern

    return Median


