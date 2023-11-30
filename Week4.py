from Week2 import HammingDistance
from Week3 import ProfileMostProbableKmer,Score,Consensus,Count,Pr
import numpy as np
import random
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {'A': np.ones(k).tolist(), 'T': np.ones(k).tolist(), 'C': np.ones(k).tolist(), 'G': np.ones(k).tolist()}
    for i in range(k):
        for j in range(t):
            count[Motifs[j][i]][i] = count[Motifs[j][i]][i] + 1
    return count
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    count = CountWithPseudocounts(Motifs)
    profile={'A':(np.array(count['A'])/(t+4)).tolist(),'T':(np.array(count['T'])/(t+4)).tolist(),'C':(np.array(count['C'])/(t+4)).tolist(),'G':(np.array(count['G'])/(t+4)).tolist()}
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
def Motifs(Profile, Dna):
    Motifs_out = []
    k = len(Profile)
    t = len(Dna)
    for i in range(t):
        Motifs_out.append(ProfileMostProbableKmer(Dna[i],k,Profile))
    return Motifs_out
def RandomMotifs(Dna, k, t):
    k_mers= []
    m = len(Dna[0])
    for i in range(t):
        idx = random.randint(0,m-k)
        k_mers.append(Dna[i][idx:idx+k])
    return k_mers
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna,k,t)
    BestMotif = M
    score = Score(BestMotif)
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile,Dna)
        if Score(M) < Score(BestMotif):
            BestMotif = M
            score.append(Score(BestMotif))
        else:
            return BestMotif
def Normalize(Probabilities):
    p_sum = sum(Probabilities.values())
    for k,p in Probabilities.items():
        Probabilities[k] = p/p_sum
    return Probabilities
def WeightedDie(Probabilities):
    sort_p = list(Probabilities.values())
    sort_p.sort()
    temp_p = np.array(sort_p)
    r = random.uniform(0,1)
    flag = 0
    for i in range(len(sort_p)):
        if r> temp_p[i]:
            temp_p[i+1] = temp_p[i+1] + temp_p[i]
        else:
            m = sort_p[i]
            flag =1
            break
    if flag == 0:
        m = sort_p[-1]
    for k,p in Probabilities.items():
        if p == m:
            return k
P = {"AATC":0.23,"AAAC":0.5,"ATTC":0.27}
print(WeightedDie(P))