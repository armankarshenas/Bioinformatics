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
def Motifs(Profile, Dna,k):
    Motifs_out = []
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
    flag = 1
    while flag == 1:
        print(score)
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile,Dna)
        if Score(M) < Score(BestMotif):
            BestMotif = M
            score = Score(BestMotif)
        else:
            flag = 0
    return BestMotif
def Normalize(Probabilities):
    p_sum = sum(Probabilities.values())
    for k,p in Probabilities.items():
        Probabilities[k] = p/p_sum
    return Probabilities
def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    # Iterate through k-mers and their probabilities
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer
def ProfileGeneratedString(Text,profile,k):
    n = len(Text)
    prob = {}
    for i in range(0,n-k+1):
        prob[Text[i:i+k]] = Pr(Text[i:i+k],profile)
    prob = Normalize(prob)
    return WeightedDie(prob)
def GibbsSampler(Dna,k,t,N):
    motif = RandomMotifs(Dna,k,t)
    Bestmotif = motif
    for i in range(N):
        p = random.randint(0,t-1)
        motif.pop(p)
        profile = ProfileWithPseudocounts(motif)
        m = ProfileGeneratedString(Dna[p],profile,k)
        motif.insert(p,m)
        if Score(motif) < Score(Bestmotif):
            Bestmotif = motif
    return Bestmotif
P = {"A":[1/8, 2/8, 3/8],"C":[1/8,1/8,1/8],"T":[3/8,2/8,3/8],"G":[4/8,3/8,1/8]}
Dna = ["TGACGTTC", "TAAGAGTT", "GGACGAAA","CTGTTCGC"]
for i in range(len(Dna)):
    print(ProfileMostProbableKmer(Dna[i],3,P))

pr = np.array([0.15,0.6,0.225,0.225,0.3])
print(pr/np.sum(pr))
