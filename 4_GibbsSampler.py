'''
    My solutions for the coding exercises in week 4 of Coursera's Bioinformatics 1.
    Copyright 2016 Alexander Wood.
    
    Includes:
    Randomized Motif Search 
    Gibbs Sampler
'''


import random
import 3_Regulatory_Motifs as m3
import sys
import 2_PatternsAndFrequencies as me
import numpy


##  INPUT: Collection of strings DNA, int k, int t
##  OUTPUT: Collection of strings BestMotifs, with pseudocounts,
##  Resulting from applying RandomizedMotifSearch 1000 times.
def RandomizedMotifSearch(DNA, k, t):
    Motifs = [] # Initialize random motifs
    i = 1##REMOVE

    # Randomly select a k-mer from each string & store it as Motifs.
    for dna in DNA:
        r = random.randint(0, len(DNA[0]) - k)
        Motifs.append(dna[r:r+k])

    BestMotifs = Motifs # Initialize the BestMotifs to the random Motifs.

    BestMotifs = ['GTC', 'CCC', 'ATA', 'GCT'] ##REMOVE



    while True:
        # Apply Laplace's Rule of Succession to form Profile.
        CountMatrix = m3.Count(Motifs)
        
        for row in range(len(CountMatrix)):
            for col in range(len(CountMatrix[0])):
                CountMatrix[row][col] += 1
                
        Prof_Matrix = m3.Profile(CountMatrix)

        # Now, construct the Motifs matrix to be the collection of
        # profile-most probable k-mers under Prof_Matrix
        Motifs = []
        for dna in DNA:
            PMK = m3.ProfileMost_kmer(dna, Prof_Matrix, k)
            Motifs.append(PMK)

        if m3.Score(Motifs) < m3.Score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs

        print(i, BestMotifs)##REMOVE
        i += 1##REMOVE




def RandomizedMotifSearch_CallAndRepeat(DNA, k, t, N):
    most_common_matrix = []
    score = sys.maxsize
    ent = sys.maxsize
        
    for n in range(N):
        answers = RandomizedMotifSearch(DNA, k, t)

        score_temp = m3.Score(answers)

        if score_temp < score:
            score = score_temp
            final_score = answers
            
    return final_score


        


##  INPUT: Profile matrix, and a DNA string
##  OUTPUT: A profile-randomly generated k-mer in DNAString
def ProfileRandomlyGeneratedKmer(Prof_Matrix, DNAString, k):
    Probabilities = []
    Strings = []

    for i in range(len(DNAString) - k + 1):
        Strings.append(DNAString[i:i+k])
    
    for string in Strings:
        probability = m3.ProbabilityGivenProfile(Prof_Matrix, string)
        Probabilities.append(probability)
    
    SUM = sum(Probabilities)
    weight_matrix = []
    for probability in Probabilities:
        weight_matrix.append(probability/SUM)

    VALUE = numpy.random.choice(Strings, p = weight_matrix)

    return VALUE



def GibbsSampler(DNA, k, t, N):
    Motifs = []
    # Randomly select k-mers Motifs in each string from DNA
    for dna in DNA:
        n = random.randint(0, len(DNA[0]) - k)
        Motifs.append(dna[n:n+k])

    # Initialize BestMotifs to Motifs
    BestMotifs = Motifs.copy()


    for j in range(1, 5*N+1):
        i = random.randint(0, t-1)

        # Profile = profile matrix constructed from all strings
        # in Motif EXCEPT Motif[i]
        Motifs_Cut = Motifs[:i] + Motifs[i+1:]
        
        # Apply Laplace's Rule of Succession to form Profile from Motifs_Cut.       
        CountMatrix = m3.Count(Motifs_Cut)
        for row in range(len(CountMatrix)):
            for col in range(len(CountMatrix[0])):
                CountMatrix[row][col] += 1
        Prof_Matrix = m3.Profile(CountMatrix)

        I = ProfileRandomlyGeneratedKmer(Prof_Matrix, DNA[i], k)

        Motifs[i] = I

        if m3.Score(Motifs) < m3.Score(BestMotifs):
            BestMotifs = Motifs.copy()
            SCORE = m3.Score(BestMotifs)

    return BestMotifs


