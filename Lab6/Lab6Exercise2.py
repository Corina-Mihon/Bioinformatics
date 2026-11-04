# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:32:03 2026

@author: Corina
"""

import gzip

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


def reverse_complement(seq):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(pairs[b] for b in seq[::-1])


def detect_transposons(dna):
    tir_length = 6
    dr_length = 4
    min_core = 10
    max_core = 150

    count = 0

    for start in range(len(dna)):
        dr = dna[start:start + dr_length]
        tir = dna[start + dr_length:start + dr_length + tir_length]

        if len(dr) != dr_length or len(tir) != tir_length:
            continue

        if any(b not in "ACGT" for b in dr + tir):
            continue

        for core_length in range(min_core, max_core + 1):
            right = start + dr_length + tir_length + core_length

            right_tir = dna[right:right + tir_length]
            right_dr = dna[right + tir_length:right + tir_length + dr_length]

            if (right_tir == reverse_complement(tir)
                    and right_dr == dr):

                count += 1
                print(f"Candidate {count}: Start={start + dr_length + 1}, "
                      f"End={right + tir_length}")
                print(f"  TIR: {tir}, DR: {dr}, Core: {core_length} bp")

    print("Total candidates:", count)

sequences = read_fasta(r"D:\Scoala\UNIVERSITATE\AnulIV-sem1\Bioinformatics\Lab5\GCF_002194405.1_ViralProj390412_genomic.fna.gz")

for i, dna in enumerate(sequences):
    print(f"\nSequence {i+1} ({len(dna)} bases)")
    detect_transposons(dna)