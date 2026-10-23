# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:16:52 2026

@author: Corina
"""


import random
import gzip

file_path = r"GCF_003575385.1_ASM357538v1_genomic.fna.gz"
output_path = r"reconstructed.fasta"

SAMPLE_COUNT = 20000
MIN_LENGTH = 100
MAX_LENGTH = 150
MIN_OVERLAP = 20


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



def generate_samples(genome, sample_count=20000,
                     min_length=100, max_length=150):

    samples = []
    position = 0

    while position < len(genome) and len(samples) < sample_count:
        length = random.randint(min_length, max_length)
        sample = genome[position:position + length]

        if len(sample) < min_length:
            break

        samples.append(sample)
        position += length

    random.shuffle(samples)

    return samples


def find_overlap(a, b):
    max_length = min(len(a), len(b))

    for length in range(max_length, MIN_OVERLAP - 1, -1):
        if a.endswith(b[:length]):
            return length

    return 0


def assemble_genome(samples):
    remaining = samples.copy()
    assembled = remaining.pop(0)

    while remaining:
        best_overlap = 0
        best_index = -1
        best_fragment = ""

        for i, fragment in enumerate(remaining):
            overlap = find_overlap(assembled, fragment)

            if overlap > best_overlap:
                best_overlap = overlap
                best_index = i
                best_fragment = fragment

        if best_index == -1:
            assembled += "\n"
            assembled += remaining.pop(0)
        else:
            assembled += best_fragment[best_overlap:]
            remaining.pop(best_index)

    return assembled


def main():
    print("Reading original genome...")
    sequences = read_fasta(file_path)
    genome = "".join(sequences)

    print("Original genome length:", len(genome))
    samples = generate_samples(genome)
    reconstructed = assemble_genome(samples)

    with open(output_path, "w", encoding="utf-8") as file:
        file.write(">reconstructed_genome\n")
        file.write(reconstructed.replace("\n", "") + "\n")

    print("\nReconstruction finished.")
    print("Reconstructed length:", len(reconstructed.replace("\n", "")))
    print("Saved to:", output_path)

    reconstructed_clean = reconstructed.replace("\n", "")

    if reconstructed_clean == genome:
        print("The reconstructed genome matches the original exactly.")
    elif reconstructed_clean in genome:
        print("The reconstructed sequence is a substring of the original.")
    else:
        print("The reconstructed sequence differs from the original.")
        print("This can happen because of incorrect or ambiguous overlaps.")


if __name__ == "__main__":
    main()
