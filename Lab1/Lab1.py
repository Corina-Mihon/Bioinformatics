# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:13:38 2026

@author: Corina
"""

import matplotlib.pyplot as plt

frequency = {
    "A" : 0,
    "C" : 0,
    "G" : 0,
    "T" : 0
    }

S = "ATCGACGCTAACGTA"

for ch in S:
    if ch == "A":
        frequency["A"] +=1
    elif  ch == "C":
        frequency["C"] +=1 
    elif  ch == "G":
        frequency["G"] +=1 
    elif  ch == "T":
        frequency["T"] +=1 

print(frequency)


def freq (S) :
    combinations_2 = {}
    for c in range(len(S)-1):
        comb_2 = S[c : c + 2]
        if comb_2 in combinations_2:
            combinations_2[comb_2] +=1
        else:
            combinations_2[comb_2] = 1
    
    plt.figure()
    plt.bar(combinations_2.keys(), combinations_2.values(), color="blue")
    plt.show()
    
    combinations_3 = {}
    for c in range(len(S)-1):
        comb_3 = S[c : c + 3]
        if comb_3 in combinations_3:
            combinations_3[comb_3] +=1
        else:
            combinations_3[comb_3] = 1
            
    plt.figure()  
    plt.bar(combinations_3.keys(), combinations_3.values(), color="blue")
    plt.show()
        
freq(S)