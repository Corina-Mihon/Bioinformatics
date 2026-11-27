# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:14:58 2026

@author: Corina
"""

import gzip
import math
import csv

# Cele 7 secvente din tabel
patterns = [
    "ATTCCGAAGTG",
    "TACGCAATCGA",
    "CGTGCTTACAG",
    "GTAATCTGAGC",
    "TGATCAGGCTA",
    "GACCTAAGGAT",
    "TCGATATGCCG"
]


def read_fasta(path):
    sequences = []

    with open(path, "rb") as f:
        compressed = f.read(2) == b"\x1f\x8b"

    opener = gzip.open if compressed else open

    with opener(path, "rt", encoding="utf-8") as f:
        sequence = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence:
                    sequences.append("".join(sequence))
                    sequence = []
            else:
                sequence.append(line.upper().replace("U", "T"))

        if sequence:
            sequences.append("".join(sequence))

    return sequences


# Calculate probabilities
probabilities = []

for i in range(11):
    column = {}

    for base in "ACGT":
        count = sum(p[i] == base for p in patterns)
        column[base] = count / 7

    probabilities.append(column)


# Calculate log-odds scores
scores = []

for column in probabilities:
    score_column = {}

    for base in "ACGT":
        p = column[base]

        if p == 0:
            score_column[base] = -20
        else:
            score_column[base] = math.log(p / 0.25)

    scores.append(score_column)


# Read DNA from FASTA
path = input("FASTA file path: ").strip().strip('"')
sequences = read_fasta(path)


with open("scores.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Sequence", "Window", "Start", "End", "DNA", "Score"])

    for seq_no, dna in enumerate(sequences, 1):
        print(f"\nDNA sequence {seq_no}")

        best_score = float("-inf")
        best_position = 0

        for i in range(len(dna) - 10):
            window = dna[i:i + 11]

            score = sum(
                scores[j].get(base, -20)
                for j, base in enumerate(window)
            )

            writer.writerow([
                seq_no, f"S{i+1}", i+1, i+11,
                window, round(score, 2)
            ])

            if i < 20:
                print(f"S{i+1} = {score:.2f}")

            if score > best_score:
                best_score = score
                best_position = i + 1

        if best_position:
            print("Best score:", round(best_score, 2))
            print("Best position:", best_position)

print("\nAll window scores saved in scores.csv")