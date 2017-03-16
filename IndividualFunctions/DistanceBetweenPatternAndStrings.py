# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.




import itertools
import sys
import numpy

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]



def HammingDistance(p, q):
    if len(p) != len(q):
        return 'Invalid input'

    hamming_dist = 0 # Initialize to zero

    for i in range(len(p)):
        if p[i] != q[i]:
            hamming_dist += 1
    #print hamming_dist
    return hamming_dist



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
