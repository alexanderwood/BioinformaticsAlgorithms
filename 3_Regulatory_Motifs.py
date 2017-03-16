'''
    Following is my code from week 3 of Coursera's Bioinformatics 1 course.
    Copyright 2016 Alexander Wood.
    
    Includes the following functions:
    Motif Count
    Profile Matrix
    Consensus Vector
    Brute-Force Algorithm for Finding Regulatory motifs
    Faster Methods for finding Regulatory Motifs
'''



# Uses modules from previous weeks.

import 2_PatternsAndFrequencies as m2 # Module from week 2
import itertools
import numpy
import sys


# Store as global variable since its values must be constant.
keys = 'ACGT'


##
##    Next we create the Count matrix, a 4 x k matrix which
##    stores the number of times a nucleotide appears in each column.
##
##    For example, input:
##    AAA
##    ACG
##    AGA
##    ATG
##
##    Yields:
##    Count = [ [4, 1, 2],   # of times A appears in each col
##              [0, 1, 0],   # of times C appears in each col
##              [0, 1, 2],   # of times G appears in each col
##              [0, 1, 0] ]  # of times T appears in each col
##

# INPUT: Motif matrix
# OUTPUT: Count matrix
def Count(Motif):
    Count = []

    # Initialize Count as a 4 x k matrix filled with zeros
    for i in range(len(keys)):
        temp_row = []
        for j in range(len(Motif[0])):
            temp_row.append(0)
        Count.append(temp_row)
    
    # If the argument given is the file path name, call Motif() and
    # generate its Motif matrix.
    if type(Motif) == str:
        Motif = Motifs(Motif)

    print(Motif)
    
    for j in range(len(Motif[0])):         # Vary over columns...
        for i in range(len(Motif)):        # Vary over rows...
            for k in range(len(keys)):     # Vary over nucleotides...
                print('(i,j)', i, j)
                if Motif[i][j] == keys[k]: 
                    Count[k][j] += 1

    return Count   



                

##    Now we wish to compute the Profile of the Motif matrix, which is
##    just the count matrix converted to percentages. Thus we wish to divide
##    each entry in the count matrix by the total number of rows. The element
##    Profile[i][j] represents the percentage of nucleotides in column j which
##    take on value keys[i].
##    For instance:
##        INPUT:
##    2 0 2
##    4 4 5
##    1 1 1
##    3 5 2
##
##        OUTPUT:
##    .2  0 .2
##    .4 .4 .5
##    .1 .1 .1
##    .3 .5 .2


# INPUT: Count matrix
# OUTPUT: Profile matrix
def Profile(Count):
    N = list(numpy.sum(Count, 0))
    N = N[0]

    for i in range(len(Count)):
        for j in range(len(Count[i])):
            Count[i][j] /= N

    Profile = Count
    return Profile


# INPUT: Profile matrix
# OUTPUT: Consensus vector where Consensus[j] yields the most common letter in column j.
def Consensus(Profile):
    Consensus = []

    print('Prof:', Profile)

    for j in range(len(Profile[0])): #Column
        temp = []
        for i in range(len(Profile)): #Rows
            print('(i,j)', i, j)
            temp.append(Profile[i][j])

        letter_max = max(set(temp), key=temp.count)
        index_max = keys.index(letter_max)

        Consensus.append(keys[index_max])

    return Consensus










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
                    
                    temp = mod2.ApproximatePatternMatching(Mismatch, dna, d)

                    if temp: # If it appears in the strand
                        FoundFlags[index] += 1
                if 0 not in FoundFlags:
                    Patterns.append(Mismatch)          
                
    return list(set(Patterns))


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
                
            condition = m2.HammingDistance(Pattern, test_string)
            if HammingDistanceMax > condition:
                HammingDistanceMax = condition
        distance = distance + HammingDistanceMax

    return distance





##
##     INPUT: Int k, list of strings DNA.
##     OUTPUT: k-mer Pattern that minimizes d(Pattern, DNA) among all k-mers Pattern.
##             (Return the first one if there are multiple.)
##    
def MedianString(DNA, k):
    distance = sys.maxsize  # Initialize distance to "infinity"
    Median = ''

    for val in range(4**k):
        Pattern = m2.NumberToPattern(val, k)
        test = DistanceBetweenPatternAndStrings(DNA, Pattern)

        if distance > test:
            distance = test
            Median = Pattern

    return Median



##  INPUT: A profile matrix Profile, a DNA pattern Pattern
##  OUTPUT: The probability that Profile generates Pattern.

def ProbabilityGivenProfile(Profile, Pattern):
    Probability = 1 # Initialize to 1

    #print('pattern:', Pattern)
    for j in range(len(Pattern)):
        i = keys.index(Pattern[j])
        Probability *= Profile[i][j]

    return Probability



##  Profile-most k-mer problem.
##  INPUT: String Text, int k, 4xk matrix Profile
##  OUTPUT: the k-mer that was most likely to have been generated by Profile among all k-mers in Text
def ProfileMost_kmer(Text, Profile, k):   
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i+k]
        Probability = ProbabilityGivenProfile(Profile, Pattern)

        if i == 0:
            PMK = Pattern
            ProbPMK = Probability

        elif Probability >= ProbPMK:
            PMK = Pattern
            ProbPMK = Probability

    return PMK
    



# INPUT: Motifs
##  OUTPUT: Score(motifs), the total number of mismatches.
def Score(Motifs):
    Consen = Consensus(Motifs)
    Score = 0

    for j in range(len(Motifs[0])): # Column
        for i in range(len(Motifs)): # Rows
            if Motifs[i][j] != Consen[j]:
                Score += 1

    return Score



##  INPUT: Integers k, t and list of strings DNA
##  OUTPUT: Collection of strings BestMotifs.
##  This algorithm runs WITHOUT pseudocounts, vastly damaging its reliability.
##  We improve on this algorithm after.
def GreedyMotifSearch_noPC(DNA, k, t):
    BestMotifs = [] # BestMotifs will be a motif matrix constructed from the first kmer in each string in DNA

    for dna in DNA:
        BestMotifs.append(dna[:k])

    FirstString = DNA[0]
    
    for i in range(len(FirstString) - k + 1):
        Motif = []
        Motif.append(FirstString[i:i+k])   # So, for each k-mer in DNA[0]...
        ProfMatrix = Profile(Count(Motif)) # ...Build profile from this k-mer

        for l in range(1, t):
            NextMotif = ProfileMost_kmer(DNA[l], ProfMatrix, k)
            Motif.append(NextMotif)
            ProfMatrix = Profile(Count(Motif))

        if Score(Motif) < Score(BestMotifs):
            BestMotifs = Motif

    return BestMotifs


##  INPUT: Integers k, t and list of strings DNA
##  OUTPUT: Collection of strings BestMotifs.
##  This algorithm runs WITH pseudocounts
def GreedyMotifSearch(DNA, k, t):
    BestMotifs = [] # BestMotifs will be a motif matrix constructed from the first kmer in each string in DNA

    for dna in DNA:
        BestMotifs.append(dna[:k])

    FirstString = DNA[0]
    
    for i in range(len(FirstString) - k + 1):
        Motif = []
        Motif.append(FirstString[i:i+k])   # So, for each k-mer in DNA[0]...
        CountMatrix = Count(Motif)               # ... apply LaPlace's rule to implement Pseudocounts....
        for row in range(len(CountMatrix)):
            for col in range(len(CountMatrix[0])):
                CountMatrix[row][col] += 1
        ProfMatrix = Profile(CountMatrix)             # ...Build profile from this k-mer

        for l in range(1, t):
            NextMotif = ProfileMost_kmer(DNA[l], ProfMatrix, k)
            Motif.append(NextMotif)
            CountMatrix = Count(Motif)
            for row in range(len(CountMatrix)):
                for col in range(len(CountMatrix[0])):
                    CountMatrix[row][col] += 1
            ProfMatrix = Profile(CountMatrix)

        if Score(Motif) < Score(BestMotifs):
            BestMotifs = Motif

    return BestMotifs



    





    
