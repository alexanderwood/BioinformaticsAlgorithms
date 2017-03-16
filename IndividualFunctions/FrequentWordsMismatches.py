# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.




# Converts single digit to single symbol.
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
    


# Final, recursive, attempt:
def PatternToNumber(Pattern):
    if Pattern == '':
        return 0

    LastSymbol = Pattern[len(Pattern) - 1]
    Prefix = Pattern[0:len(Pattern)-1]

    return 4*PatternToNumber(Prefix) + SymbolToNumber(LastSymbol)




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




# Is called recursively until it returns a list of all d-neighbors. 
# i is current length of neighbor candidates.
def Neighbors_Test(Pattern, Neighbors, i, d):
    neighbors = []
    
    if d == 0:
        return [Pattern]
      
    nucs = ['A', 'C', 'G', 'T']
    
    candidates = []           
    for Neighbor in Neighbors:
        for nuc in nucs:
            candidates.append(Neighbor + nuc)

    new_candidates = []

    for candidate in candidates:
        if (HammingDistance(candidate, Pattern[:i+1]) < d) and (len(candidate) < len(Pattern)):
            new_candidates.append(candidate)


        elif (HammingDistance(candidate, Pattern[:i+1]) < d) and (len(candidate) == len(Pattern)):
            neighbors.append(candidate)

        elif HammingDistance(candidate, Pattern[:i+1]) == d:
            neighbors.append(candidate + Pattern[i+1:])

        else:
            continue

    if i < len(Pattern):
        more_neighbors = Neighbors_Test(Pattern, new_candidates, i+1, d)
        neighbors += more_neighbors

    elif i == len(Pattern):
        even_more = []
        for cand in new_candidates:
            for nuc in nucs:
                even_more.append(cand + nuc)
        neighbors += even_more

    return neighbors
 
# This function initializes the call to Neighbors_Test, which will run recursively until
# all d-neighbors have been found. 
def Neighbors(Pattern, d):
    neighbors_init = ['A', 'G', 'C', 'T']

    return Neighbors_Test(Pattern, neighbors_init, 1, d)










'''
A most frequent k-mer with up to d mismatches in Text is simply a string
Pattern maximizing Countd(Text, Pattern) among all k-mers. Note that
Pattern does not need to actually appear as a substring of Text; for
example, as we already saw, AAAAA is the most frequent 5-mer with 1 mismatch
in AACAAGCTGATAAACATTTAAAGAG, even though it does not appear exactly in this
string. Keep this in mind while solving the following problem.

Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
     Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
     Output: All most frequent k-mers with up to d mismatches in Text.
'''

def FrequentWordsMismatch(Text, k, d):
    FrequentPatterns = [] # Initialize an empty set to store the frequent patterns found

    close = []
    frequency_array = []
    
    for i in range(4**k):
        close.append(0)
        frequency_array.append(0)

    for i in range(len(Text) - k + 1):
        neighborhood = Neighbors(Text[i:i+k], d)

        for neighbor in neighborhood:
            index = PatternToNumber(neighbor)
            close[index] += 1

    for i in range(4**k):
        if close[i] > 0:
            Pattern = NumberToPattern(i, k)
            frequency_array[i] = len(ApproximatePatternMatching(Pattern, Text, d))

    max_count = max(frequency_array)

    for i in range(4**k):
        if frequency_array[i] == max_count:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)

    return FrequentPatterns

