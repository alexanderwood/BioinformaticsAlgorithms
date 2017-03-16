# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.




'''           
Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
     Input: A DNA string Genome.
     Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
'''


# Recall that skew(i, Genome) = # of occureences of G - # occurences of C
# on the first i nucleotides in the genome.
def skew(i, Genome):
    # Store the values of skew in a list, Skew
    # So, skew_list[j] is the value of skew at index j
    # start with skew_list = [0] since the balance of G and C is initally zero.
    skew_list = [0]

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
