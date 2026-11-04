import random

random.seed(42)

def dna(n):
    return ''.join(random.choice("ACGT") for _ in range(n))

def reverse_complement(s):
    return ''.join({"A": "T", "T": "A", "C": "G", "G": "C"}[b]
                   for b in s[::-1])

def make_transposon():
    dr = dna(4)
    tir = dna(6)
    core = dna(20)

    return dr + tir + core + reverse_complement(tir) + dr

sequence = list(dna(320))

transposons = [make_transposon() for _ in range(4)]

starts = [30, 100, 114, 128]

for te, start in zip(transposons, starts):
    sequence[start:start + len(te)] = te

sequence = ''.join(sequence)

def detect(te):
    return max(
        range(len(sequence) - len(te) + 1),
        key=lambda i: sum(
            a == b for a, b in zip(te, sequence[i:i + len(te)])
        )
    )

print("DNA sequence:")
print(sequence)

print("\nTransposons:")
for i, te in enumerate(transposons):
    print(f"TE{i+1}: {te}")

print("\nDetected positions:")
for i, te in enumerate(transposons):
    found = detect(te)
    print(f"TE{i+1}: Start = {found+1}, End = {found+len(te)}")

print("\nOverlaps:")
for i in range(4):
    for j in range(i + 1, 4):
        left = max(starts[i], starts[j])
        right = min(starts[i] + 40, starts[j] + 40)

        if left < right:
            print(f"TE{i+1} and TE{j+1}: {left+1}-{right}")
            