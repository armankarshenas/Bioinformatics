import numpy as np
from Week2 import HammingDistance
def Count(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {'A':np.zeros(k).tolist(),'T':np.zeros(k).tolist(),'C':np.zeros(k).tolist(),'G':np.zeros(k).tolist()}
    for i in range(k):
        for j in range(t):
            count[Motifs[j][i]][i] = count[Motifs[j][i]][i] + 1
    return count
def Profile(Motifs):
    # insert your code here
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {'A':np.zeros(k).tolist(),'T':np.zeros(k).tolist(),'C':np.zeros(k).tolist(),'G':np.zeros(k).tolist()}
    for i in range(k):
        for j in range(t):
            profile[Motifs[j][i]][i] = profile[Motifs[j][i]][i] + 1/t
    return profile
def Consensus(Motifs):
    count = Count(Motifs)
    cons_seq = ""
    k = len(Motifs[0])
    for i in range(k):
        local_count = [count['A'][i],count['T'][i],count['C'][i],count['G'][i]]
        m = max(local_count)
        if count['A'][i] == m:
            cons_seq = cons_seq + "A"
        elif count['T'][i] == m:
            cons_seq = cons_seq + "T"
        elif count['C'][i] == m:
            cons_seq = cons_seq + "C"
        else:
            cons_seq = cons_seq + "G"
    return cons_seq
def Score(Motifs):
    cons_seq = Consensus(Motifs)
    t = len(Motifs)
    score = 0
    for i in range(t):
        score = score + HammingDistance(cons_seq,Motifs[i])
    return score
def Pr(Text,Profile):
    p = 1
    for i in range(len(Text)):
        p = p*Profile[Text[i]][i]
    return p
def  ProfileMostProbableKmer(Text, k, Profile):
    mp = Pr(Text[0:k],Profile)
    mfp = Text[0:k]
    for i in range(1,len(Text)-k+1):
        p = Pr(Text[i:i+k],Profile)
        if p>mp:
            mfp = Text[i:i+k]
            mp = p
    return mfp

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ComputeEntropy(Motifs):
    P = Profile(Motifs)
    score = 0
    for i in range(len(Motifs[0])):
        score = score - P['A'][i]*np.log2(P['A'][i])-P['T'][i]*np.log2(P['T'][i]) - P['C'][i]*np.log2(P['C'][i]) - P['G'][i]*np.log2(P['G'][i])
    return score
