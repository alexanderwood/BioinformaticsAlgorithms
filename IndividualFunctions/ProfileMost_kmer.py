# Copyright Alexander Wood 2016.
# Solutions for Coursera's Bioinformatics 1 Course.



        

Temp = []
for i in range(2, len(tokens)):
    Temp.append(float(tokens[i]))

n = len(Temp)//4
Profile = list(chunks(Temp,n))




# Store as global variable since its values must be constant.
keys = 'ACGT'



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


