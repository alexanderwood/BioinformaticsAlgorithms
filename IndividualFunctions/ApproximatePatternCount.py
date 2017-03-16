# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.



# The following is a variation on the HammingDistance function which quits running its calculations once a bound
# on the hamming distance has been reached. So for instance if the bound is d and p and q have a hamming distance of
# anything greater than d, it will simply return d+1.

def HammingDistance_bound(p, q, d):
    if len(p) != len(q):
        return 'Invalid input'

    hamming_dist = 0 # Initialize to zero

    for i in range(len(p)):
        if p[i] == q[i]:
            continue
        elif p[i] != q[i] and hamming_dist <= d:
            hamming_dist += 1
            continue
        else:
            break
    #print hamming_dist
    return hamming_dist


'''
We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern'
of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA
box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
     Input: Strings Pattern and Text along with an integer d.
     Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

Code Challenge: Solve the Approximate Pattern Matching Problem.
'''

def ApproximatePatternMatching(Pattern, Text, d):
    # Create a list to store indexes of approximate matches.
    approx_match = []

    for i in range(len(Text) - len(Pattern) + 1):
        temp = HammingDistance_bound(Pattern, Text[i: i+len(Pattern)], d)

        if temp <= d:
            approx_match.append(str(i))

    return approx_match




'''
Our goal now is to modify our previous algorithm for the Frequent Words Problem in order to find DnaA boxes by identifying
frequent k-mers, possibly with mismatches. Given strings Text and Pattern as well as an integer d, we define
Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches. For example,
Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at most one mismatch:
AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.
'''

def Count(Pattern, Text, d):
    total_count = len(ApproximatePatternMatching(Pattern, Text, d))
    return total_count
