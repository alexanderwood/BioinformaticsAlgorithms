# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.






import random
import sys
import numpy


# Store as global variable since its values must be constant.
keys = 'ACGT'


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]





##  INPUT: A profile matrix Profile, a DNA pattern Pattern
##  OUTPUT: The probability that Profile generates Pattern.

def ProbabilityGivenProfile(Profile, Pattern):
    Probability = 1 # Initialize to 1
    #print(Profile)

    for letter_index in range(len(Pattern)):
        letter = Pattern[letter_index]
        i = keys.index(letter)
        Probability *= Profile[i][letter_index]

    return Probability




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


    for j in range(len(Motif[0])):         # Vary over columns...
        for i in range(len(Motif)):        # Vary over rows...
            for k in range(len(keys)):     # Vary over nucleotides...
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
def Consensus(Profile_Mat):
    Consensus = []
    
    for j in range(len(Profile_Mat[0])): #Column
        temp = []
        for i in range(len(Profile_Mat)): #Rows
            temp.append(Profile_Mat[i][j])

        letter_max = max(set(temp), key=temp.count)
        index_max = keys.index(letter_max)

        Consensus.append(keys[index_max])

    return Consensus



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




##  INPUT: Profile matrix, and a DNA string
##  OUTPUT: A profile-randomly generated k-mer in DNAString
def ProfileRandomlyGeneratedKmer(Prof_Matrix, DNAString, k):
    Probabilities = []
    Strings = []

    for i in range(len(DNAString) - k + 1):
        Strings.append(DNAString[i:i+k])
    
    for string in Strings:
        probability = ProbabilityGivenProfile(Prof_Matrix, string)
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

    for j in range(1, 50*N+1):
        i = random.randint(0, t-1)

        # Profile = profile matrix constructed from all strings
        # in Motif EXCEPT Motif[i]
        Motifs_Cut = Motifs[:i] + Motifs[i+1:]
        
        # Apply Laplace's Rule of Succession to form Profile from Motifs_Cut.       
        CountMatrix = Count(Motifs_Cut)
        for row in range(len(CountMatrix)):
            for col in range(len(CountMatrix[0])):
                CountMatrix[row][col] += 1
        Prof_Matrix = Profile(CountMatrix)

        I = ProfileRandomlyGeneratedKmer(Prof_Matrix, DNA[i], k)

        Motifs[i] = I

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs.copy()
            SCORE = Score(BestMotifs)

    return BestMotifs

