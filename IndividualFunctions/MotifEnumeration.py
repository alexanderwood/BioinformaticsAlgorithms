# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.




import itertools
import sys


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
Code Challenge: Implement MotifEnumeration (reproduced below).
     Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna.
'''
# A brute-force algorithm for finding regulatory motifs in DNA strands.

# Notes that a regulatory (k,d)-motif is a k-mer which appears in every
# string from some collection of strings with at most d mismatches.
# It appears in scattered regions within the DNA strand.

# A brute-force algorithm like the one which follows is highly
# unlikely to be of any real use but is an interesting learning exercise.

# First we create the following function:

# INPUT: A k-mer Text, integer d.
# OUTPUT: Every possible k-mer which is a d-mismatch of Text.
def GenerateMismatches(Text, d):
    TEXT = list(Text)

    for position in itertools.combinations(range(len(Text)), d): # Every way to pick d letters from Text
        for key in itertools.product(keys, repeat=d): # Every way to pick d chars from keys
            Flag = False # Escape flag

            ans = dict(zip(position, key)) # Create a dictionary which associates the
                                           # string position with the key you wish to overwrite.
            #print(ans)
            yield ''.join([TEXT[p] if p not in position else ans[p] for p in range(len(Text))])

            
#ans = GenerateMismatches('AAA', 2)
#for an in ans:
#    print(an, end=' ')
    
# INPUT: Integers k,d; List of strings DNA
# OUTPUT: List Motifs of all (k,d)-motifs within DNA
def MotifEnumeration(DNA, k, d):
    Patterns = [] # Initialize an empty list to store patterns.

    for Strand in DNA:                          # For each strand....           
        for i in range(len(Strand) - k + 1):    # For each k-mer in each strand...
            Pattern = Strand[i:i+k]
            
            Mismatches = list(set(GenerateMismatches(Pattern, d)))

            for Mismatch in Mismatches:         # For each d-mismatch....
                FoundFlags = []
                for i in range(len(DNA)):
                    FoundFlags.append(0)
                        
                for index in range(len(DNA)):
                    dna = DNA[index]
                    
                    temp = ApproximatePatternMatching(Mismatch, dna, d)

                    if temp: # If it appears in the strand
                        FoundFlags[index] += 1
                if 0 not in FoundFlags:
                    Patterns.append(Mismatch)          
                
    return list(set(Patterns))



