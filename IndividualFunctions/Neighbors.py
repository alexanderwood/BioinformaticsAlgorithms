# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.






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

#Is called recursively until it returns a list of all d-neighbors. 
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


