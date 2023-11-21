import numpy as np

Text = "AACTGCGTGGGGCCCATA"
Pattern = "GG"

n_Text = len(Text)
n_p = len(Pattern)
sheila = 0
for i in range(n_Text-n_p+1):
    if Text[i:i+n_p] == Pattern:
        sheila = sheila +1
print(sheila)

numbers = np.array([1,2,3,4])
print(numbers)
