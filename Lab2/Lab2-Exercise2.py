# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:03:31 2026

@author: Corina
"""

import math

S= "ACGCGTGCCA"

frequency = {
    "A" : 0,
    "C" : 0,
    "G" : 0,
    "T" : 0
    }


Na= 0.05

def tempMelting(S):
    for ch in S:
        if ch == "A":
            frequency["A"] +=1
        elif  ch == "C":
            frequency["C"] +=1 
        elif  ch == "G":
            frequency["G"] +=1 
        elif  ch == "T":
            frequency["T"] +=1 
    #First formula
    Tm = 4 * (frequency["G"] + frequency["C"]) + 2 * (frequency["A"] + frequency["T"])
    print("With first formula: " + str(Tm))
    
    #Second formula
    Tm_alternative = 81.5 + 16.6 * math.log10(Na) + 41 * ((frequency["G"] + frequency["C"]) / len(S)) - 600 / len(S)
    print("With second formula: " + str(Tm_alternative))
    
tempMelting(S)