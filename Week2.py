import numpy as np
from Week1 import PatternCount
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)],symbol)
    return array
print(SymbolArray("AAAAGGGG","A"))
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2],symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array
def FindSkew(Genome):
    skew = np.zeros(len(Genome)+1)
    for i in range(len(Genome)):
        if Genome[i] == "A" or Genome[i] == "T":
            skew[i+1] = skew[i]
        elif Genome[i] == "C":
            skew[i+1] = skew[i]-1
        else:
            skew[i+1] = skew[i]+1
    return skew

def SkewArray(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == "A" or Genome[i] == "T":
            skew.append(skew[i])
        elif Genome[i] == "C":
            skew.append(skew[i]-1)
        else:
            skew.append(skew[i]+1)
    return skew
def FindOri(Genome):
    skew = FindSkew(Genome)
    mn = min(skew)
    result =[]
    for i in range(len(skew)):
        if skew[i] == mn:
            result.append(i)
    return result
def HammingDistance(p, q):
    # your code here
    dist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            dist = dist +1
    return dist

def ApproximatePatternMatching(Text,Pattern,d):
    positions = []
    n_t = len(Text)
    n_p = len(Pattern)
    for i in range(n_t-n_p+1):
        if HammingDistance(Text[i:i+n_p],Pattern) <= d:
            positions.append(i)
    return positions
def ApproximatePatternCount(Text,Pattern,d):
    positions = []
    n_t = len(Text)
    n_p = len(Pattern)
    for i in range(n_t - n_p + 1):
        if HammingDistance(Text[i:i + n_p], Pattern) <= d:
            positions.append(i)
    return len(positions)
text = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
print(FindOri(text))

p = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
q = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"

print(HammingDistance(p,q))