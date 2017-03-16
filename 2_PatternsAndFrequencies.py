'''
    Module for algorithms from the week 2 work on the Bioinformatics course.
    Solutions Copyright Alexander Wood 2016.
    
    Includes:
    Hamming Distance
    Nearest d-Neighbors
    Frequent Words
    & more.
'''

# PatternCount algorithm from section 1.2 of
# Bioinformatics course.

# The input is a text string of DNA, along with
# some pattern you wish to count the number of
# occurences of within that string.

# The output is the number of occurences of the
# pattern within the string.
def PatternCount(Text, Pattern):
    # Initialize count
    count = 0

    for i in range(len(Text) - len(Pattern) + 1):
        if Pattern == Text[i:i+len(Pattern)]:
            count += 1

    # Return the number of times count appears in the text
    return count




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
            







# Recall that skew(i, Genome) = # of occureences of G - # occurences of C
# on the first i nucleotides in the genome.
def skew(i, Genome):
    # Store the values of skew in a list, Skew
    # So, skew_list[j] is the value of skew at index j
    # start with skew_list = [0] since the balance of G and C is initally zero.
    skew_list = []

    # To keep track of the skew values as we run through the string
    count = 0
    
    for nuc in range(i):
        if (Genome[nuc] == 'A') or (Genome[nuc] == 'T'):
            skew_list.append(count)

        elif Genome[nuc] == 'G':
            count += 1
            skew_list.append(count)

        elif Genome[nuc] == 'C':
            count -= 1
            skew_list.append(count)

        elif Genome[nuc] == '\n':
            continue

        else:
            print('An error occurred')

#    for val in range(len(skew_list)):
#        print(skew_list[val], end=' ')

    return skew_list

'''           
Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
     Input: A DNA string Genome.
     Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
'''

# This algorithm returns a list of the indexes of minimum values in a list.
def MinListIndex(Some_List):
    # Find the minimum
    minimum = min(Some_List)

    # Create a list to store the indexes.
    min_index = []

    for i in range(len(Some_List)):
        if Some_List[i] == minimum:
            min_index.append(i)

    return min_index
        

def MinimumSkew(Genome):
    # Get the list of skew funciton values
    skew_list = skew(len(Genome), Genome)
    
    # Get the list of indexes of occurences of min value.
    min_index = MinListIndex(skew_list)

    for val in range(len(min_index)):
        print(min_index[val], end=' ')
    #return min_index 

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



'''
A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing Countd(Text, Pattern) among all k-mers.

Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
     Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
     Output: All most frequent k-mers with up to d mismatches in Text.

One way to solve the above problem is to generate all 4k k-mers Pattern, compute ApproximatePatternCount(Text, Pattern, d)
for each k-mer Pattern, and then output k-mers with the maximum number of approximate occurrences. This is an inefficient
approach in practice, since many of the 4k k-mers that this method analyzes should not be considered because neither they
nor their mutated versions (with up to d mismatches) appear in Text. Check out Charging Station: Solving the Frequent Words
with Mismatches Problem to learn about a better approach that avoids analyzing such hopeless k-mers.
'''

# We need a method of analyzing the k-mers in the Text without having to test every possible k-mer as this would be massively
# ineffecient. Therefore, we first generate the d-neighborhood of a string, ie, all strings which have at most d mismatches with
# the string.
def ImmediateNeighbors(Pattern): # ie, differ by one char 
    nucs = ['A', 'C', 'G', 'T']
    
    neighborhood = [Pattern] # Initialize a list

    for i in range(len(Pattern)):
        symbol = Pattern[i]
        nucs.remove(symbol)
        tempnucs = list(nucs)
        nucs.append(symbol)
        for nuc in tempnucs:
            neighbor = Pattern[0:i] + nuc + Pattern[i+1:]
            neighborhood.append(neighbor)
    return neighborhood


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



'''
Now we wish to perform the same task as above, but with also the reverse complement.
We now redefine the Frequent Words Problem to account for both mismatches and reverse
complements. Recall that Patternrc refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent
k-mers (with mismatches and reverse complements) in a string.
      Input: A DNA string Text as well as integers k and d.
      Output: All k-mers Pattern maximizing the sum
              Countd(Text, Pattern)+ Countd(Text, Pattern_rc)
              over all possible k-mers.
'''
def FrequentWordsMismatch_ReverseComplement(Text, k, d):
    FrequentPatterns = [] # Initialize an empty set to store the frequent patterns found

    close = [] # zero array, where 
    frequency_array = [] #initialize array to store frequency of each k-word
    
    for i in range(4**k):
        close.append(0)
        frequency_array.append(0)

    # First, find all d-neighbors of Text[i:i+k]
    for i in range(len(Text) - k + 1): 
        neighborhood = Neighbors(Text[i:i+k], d)

        # Mark these k-mers in the array, so that you only test k-mers which appear at least once.
        for neighbor in neighborhood: 
            index = PatternToNumber(neighbor)
            close[index] += 1

    for i in range(4**k): # ie, For each possible (forward) k-mer:
        if close[i] > 0: # ie, If the k-mer appears in original Text at least once
            
            Pattern = NumberToPattern(i, k) # Convert the index to a pattern.
            Pattern_Reverse = ReverseComplement(Pattern) # Find the reverse complement of this pattern
            
            # Set the frequency equal to the number of times an approximate match occurs for the original pattern
            frequency_array[i] = len(ApproximatePatternMatching(Pattern, Text, d))

            # Also add the number of times an approximate match for the reverse complement appears
            frequency_array[i] += len(ApproximatePatternMatching(Pattern_Reverse, Text, d))


    max_count = max(frequency_array)

    for i in range(4**k):
        if frequency_array[i] == max_count:
            Pattern = NumberToPattern(i, k)
            ReversePattern = ReverseComplement(Pattern)
            FrequentPatterns.append(Pattern)
            FrequentPatterns.append(ReversePattern)

    FrequentPatterns = list(set(FrequentPatterns))

    return FrequentPatterns

    



