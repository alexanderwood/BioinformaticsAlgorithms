##  The following is a module containing functions which perform computations
##  related to finding the origin of replication in DNA. Based on Week 1 of
##  the Coursera Bioinformatics Specialization Course 1.

##  Written by Alexander Wood October 2016


#   ----------------------------------------------------------------------------------------    #


##  PatternCount Algorithm
##  INPUT: String DNA, string Pattern
##  OUTPUT: The number of occurences of Pattern within DNA

def PatternCount(Text, Pattern):
    # Initialize count
    count = 0

    for i in range(len(Text) - len(Pattern) + 1):
        if Pattern == Text[i:i+len(Pattern)]:
            count += 1

    # Return the number of times count appears in the text
    return count


#   ----------------------------------------------------------------------------------------    #


##  FrequentWords Algorithm: Inefficient Version
##  INPUT: String DNA, int k
##  OUTPUT: set of all k-words (subwords of DNA of length k) which appear most frequently

def FrequentWords(Text, k):
    FrequentPatterns = set() # Empty set to store the frequent patterns in!
    Count = [] # Empty list to store each count for subword i

    for i in range(len(Text) - k + 1): # Create array of the # of occurences
        Pattern = Text[i:i+k]
        Count.append(PatternCount(Text, Pattern))

    maxCount = max(Count)

    for i in range(len(Count)):
        if Count[i] == maxCount:
            FrequentPatterns.add(Text[i:i+k])

    return FrequentPatterns


#   ----------------------------------------------------------------------------------------    #


##  PatternToNumber Algorithm
##  Input: Some DNA string Pattern
##  Output: A numerical value in one-to-one corresp. with Pattern


##  SymbolToNumber Algorithm
##  Used in PatternToNumber
##  INPUT: A, C, G, T
##  OUTPUT: 0, 1, 2, 3
def SymbolToNumber(Value):
    if Value == 'A':
        return 0
    elif Value == 'C':
        return 1
    elif Value == 'G':
        return 2
    elif Value == 'T':
        return 3
    else:
        return 'Error'


def PatternToNumber(Pattern):
    if Pattern == '':
        return 0

    LastSymbol = Pattern[len(Pattern) - 1]
    Prefix = Pattern[0:len(Pattern)-1]

    return 4*PatternToNumber(Prefix) + SymbolToNumber(LastSymbol)


#   ----------------------------------------------------------------------------------------    #


##  NumberToPattern Algorithm
##  Input: Integer Number, integer k
##  Output: The unique string Pattern of length k in one-to-one correspondence with Number

##  NumberToSymbol Algorithm
##  Used in NumberToPattern
##  OUTPUT: 0, 1, 2, 3
##  INPUT: A, C, G, T
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


def NumberToPattern(Number, k):
    if k == 1:
        return NumberToSymbol(Number)

    Prefix = Number // 4
    Remainder = Number % 4

    Pattern = NumberToSymbol(Remainder)

    PrefixPattern = NumberToPattern(Prefix, k-1)

    return PrefixPattern + Pattern


#   ----------------------------------------------------------------------------------------    #


##  ComputingFrequencies algorithm, Slow Version
##  INPUT: String Text, integer k
##  OUTPUT: Frequency array, which contains the number of times each possible pattern
##          of lencth k appears in Text.
##          We use the functions converting between number and lexicographic order
##          of DNA sequences above in order to find frequent patterns.
def ComputingFrequencies(Text, k):
    FrequencyArray = dict()

    # Initialize a dictionary with 0 for the count for each word.
    for i in range(4**k):
        FrequencyArray[i] = 0

    # Add one to the count each time you encounter a pattern.
    for i in range(len(Text) - (k-1)):
        Pattern = Text[i:i+k]
        j = PatternToNumber(Pattern)
        FrequencyArray[j] += 1

    return FrequencyArray


#   ----------------------------------------------------------------------------------------    #


##  FasterFrequentWords Algorithm
##  INPUT: String DNA, int k
##  OUTPUT: set of all k-words (subwords of DNA of length k) which appear most frequently

def FasterFrequentWords(Text, k):
    FrequentPatterns = set()
    FrequencyArray = ComputingFrequencies(Text, k)

    MaxIndex = max(FrequencyArray, key=FrequencyArray.get)
    MaxCount = FrequencyArray[MaxIndex]

    for i in range(4**k - 1):
        if FrequencyArray[i] == MaxCount:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.add(Pattern)

    return FrequentPatterns


#   ----------------------------------------------------------------------------------------    #


##  FasterFrequentWords Algorithm, lower bound variant
##  INPUT: String DNA, int k, int t
##  OUTPUT: set of all k-words (subwords of DNA of length k) which appear at least t times.

def FasterFrequentWords_t(Text, k, t):
    FrequentPatterns = set()
    FrequencyArray = ComputingFrequencies(Text, k)

    for i in range(4**k - 1):
        if FrequencyArray[i] >= t:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.add(Pattern)

    return FrequentPatterns


#   ----------------------------------------------------------------------------------------    #


##  ReverseComplement algorithm
##  INPUT: A DNA strand Text
##  OUTPUT: The reverse complement of Text (recall that A<->T and G<->C.)

def ReverseComplement(Text):
    ReverseText = '' # Initialize empty string
    
    for i in range(len(Text)):
        b = len(Text) - (i+1) # Backwards index!
        
        if Text[b] == 'A':
            ReverseText += 'T'
        elif Text[b] == 'T':
            ReverseText += 'A'
        elif Text[b] == 'G':
            ReverseText += 'C'
        elif Text[b] == 'C':
            ReverseText += 'G'
        else:
            return 'Error'

    return ReverseText


#   ----------------------------------------------------------------------------------------    #


##  PatternMatching algorithm
##  INPUT: DNA string Genome, substring Pattern
##  OUTPUT: List of indexes in Geneome where Pattern begins

def PatternMatching(Genome, Pattern):
    indexes = []
    
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i: i+len(Pattern)] == Pattern:
            indexes.append(i)

    return indexes


#   ----------------------------------------------------------------------------------------    #


##  ClumpFinding algorithm
##  INPUT: String Genome, int k, int L, int t
##  OUTPUT: All distinct k-mers forming (L,t)-clumps in Genome.

def ClumpFinding(Genome, k, L, t):
    clumps = set()
    
    for i in range(len(Genome) - L + 1):
        L_Window = Genome[i:i+L]

        words = FasterFrequentWords_t(L_Window, k, t)

        for word in words:
            clumps.add(word)

    return clumps


#   ----------------------------------------------------------------------------------------    #


##  A faster ClumpFinding algorithm
##  INPUT: String Genome, int k, int L, int t
##  OUTPUT: All distinct k-mers forming (L,t)-clumps in Genome.

def BetterClumpFinding(Genome, k, L, t):
    clumps = dict()
    FrequentPatterns = set()

    for i in range(4**k): # Inititalize all clump vals to zero
        clumps[i] = 0

    Text = Genome[0:L]

    FrequencyArray = ComputingFrequencies(Text, k)

    for i in range(4**k):
        if FrequencyArray[i] >= t:
            clumps[i] += 1

    for i in range(1,len(Genome) - L+1):
        FirstPattern = Genome[i-1:i-1+k]
        
        index = PatternToNumber(FirstPattern)

        FrequencyArray[index] = FrequencyArray[index] - 1

        LastPattern = Genome[i+L-k:i+L]

        index = PatternToNumber(LastPattern)

        FrequencyArray[index] = FrequencyArray[index] + 1

        if FrequencyArray[index] >= t:
            clumps[index] += 1

    for i in range(4**k):
        if clumps[i] >= 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.add(Pattern)
    return FrequentPatterns
            
