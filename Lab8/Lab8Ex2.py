# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:48:28 2026

@author: Corina
"""

import matplotlib.pyplot as plt

dna = "CGGACTGATCTATCTA" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"

window_size = 30


# GC Content
def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def information_content(seq):
    n = len(seq)
    total = 0

    for u in range(1, n):
        matches = 0

        for i in range(n - u):
            if seq[i] == seq[i + u]:
                matches += 1

        total += matches / ((n - u) * 100)

    return 100 - total / (n - 1)


# Sliding window
positions = []
gc_values = []
ic_values = []

for i in range(len(dna) - window_size + 1):
    window = dna[i:i + window_size]

    positions.append(i + 1)
    gc_values.append(gc_content(window))
    ic_values.append(information_content(window))


# Center of weight (using low-value regions as weights)
def center_of_weight(positions, values):
    weights = [max(values) - v for v in values]

    if sum(weights) == 0:
        return sum(positions) / len(positions)

    return sum(p * w for p, w in zip(positions, weights)) / sum(weights)


gc_center = center_of_weight(positions, gc_values)
ic_center = center_of_weight(positions, ic_values)


# Results
print("DNA length:", len(dna))
print("Number of windows:", len(positions))
print("GC Content:", round(gc_content(dna), 2))
print("Average IC:", round(sum(ic_values) / len(ic_values), 2))

print("GC center:", round(gc_center, 2))
print("IC center:", round(ic_center, 2))


# Plot 1 - GC Content
plt.figure(figsize=(10, 4))
plt.plot(positions, gc_values, color="blue")
plt.axvline(gc_center, color="red", linestyle="--", label="GC center")
plt.xlabel("Window position")
plt.ylabel("GC (%)")
plt.title("GC Content")
plt.legend()
plt.grid()
plt.show()


# Plot 2 - Information Content
plt.figure(figsize=(10, 4))
plt.plot(positions, ic_values, color="green")
plt.axvline(ic_center, color="red", linestyle="--", label="IC center")
plt.xlabel("Window position")
plt.ylabel("Information Content")
plt.title("Information Content")
plt.legend()
plt.grid()
plt.show()


# Plot 3 - Centers of weight
plt.figure(figsize=(10, 3))
plt.scatter(gc_center, 1, color="blue", s=100, label="GC Center")
plt.scatter(ic_center, 2, color="green", s=100, label="IC Center")

plt.yticks([1, 2], ["GC", "IC"])
plt.xlabel("DNA position")
plt.title("Centers of Weight")
plt.legend()
plt.grid()
plt.show()