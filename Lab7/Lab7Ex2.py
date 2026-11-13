# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 15:33:41 2026

"""

import random

random.seed(42)

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def generate_dna(n):
    return "".join(random.choices("ACGT", k=n))

repeats = ["ATGCCA", "CGTATC", "GACCTG"]


cases = {
    "Normal": [(20, 65), (110, 160), (210, 275)],
    "Overlapping": [(20, 100), (65, 150), (115, 200)],
    "Nested": [(20, 220), (60, 175), (95, 140)]
}

def detect(dna, tir):
    right_tir = reverse_complement(tir)
    results = []

    for i in range(len(dna) - len(tir) + 1):
        if dna[i:i+len(tir)] != tir:
            continue

        # At least 10 bases between the inverted repeats
        for j in range(i + len(tir) + 10,
                       len(dna) - len(tir) + 1):

            if dna[j:j+len(tir)] == right_tir:
                results.append((i + 1, j + len(tir)))

    return results


for case, positions in cases.items():
    dna = list(generate_dna(320))

    for (start, end), tir in zip(positions, repeats):
        right_tir = reverse_complement(tir)

        dna[start:start+len(tir)] = tir
        dna[end-len(tir)+1:end+1] = right_tir

    dna = "".join(dna)

    print("\nCASE:", case)
    print("DNA:", dna)

    for k, tir in enumerate(repeats, 1):
        for start, end in detect(dna, tir):
            print(f"TE{k}: Start={start}, End={end}, Length={end-start+1}")