# This is a python script for week 1 of the bioinformatics course


def PatternCount(Text,Pattern):
    count = 0
    n_p = len(Pattern)
    n_t = len(Text)
    for i in range(n_t-n_p+1):
        if Text[i:i+n_p]==Pattern:
            count = count +1
    return count

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n - k + 1):
        freq[Text[i:i + k]] = freq[Text[i:i + k]] + 1
    return freq


def FrequentWords(Text,k):
    words = []
    freq = FrequencyMap(Text,k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words
def Reverse_array(Text):
    rev_text = []
    for i in range(len(Text)):
        rev_text.append(Text[len(Text)-1-i])
    return rev_text
def Reverse(Text):
    rev_text = ""
    for char in Text:
        rev_text = char + rev_text
    return rev_text
def Complement(Text):
    comp_text = ""
    comp = {'A':'T','C':'G','T':'A','G':'C'}
    for char in Text:
        comp_text = comp_text + comp[char]
    return comp_text
def Reverse_complement(Text):
    rev_comp_text = ""
    rev_comp_text = Reverse(Text)
    rev_comp_text = Complement(rev_comp_text)
    return rev_comp_text
def PatternMatching(Pattern, Genome):
    position = []
    n_p = len(Pattern)
    n_g = len(Genome)
    for i in range(n_g-n_p+1):
        if Genome[i:i+n_p] == Pattern:
            position.append(i)
    return position
ori = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
Pattern = "TGATCA"
print(PatternCount(ori,Pattern))
print(Reverse_complement("ATTCG"))
print(FrequentWords(ori,10))