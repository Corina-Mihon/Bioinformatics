# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:34:54 2026

@author: Corina
"""

S = "ATCGCGTA"

frequency = {
    "A" : 0,
    "C" : 0,
    "G" : 0,
    "T" : 0
    }


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
    
    Tm = 4 * (frequency["G"] + frequency["C"]) + 2 * (frequency["A"] + frequency["T"])
    print(Tm)
    
tempMelting(S)