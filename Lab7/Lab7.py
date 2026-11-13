import gzip
import matplotlib.pyplot as plt

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2)
}


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


def fragments(dna, cuts):
    points = [0] + sorted(set(cuts)) + [len(dna)]
    return [points[i+1] - points[i] for i in range(len(points)-1)]


def draw_gel(lanes):
    fig, ax = plt.subplots()

    for x, (name, sizes) in enumerate(lanes.items()):
        for size in sizes:
            ax.hlines(size, x-0.3, x+0.3, color="blue", linewidth=2)

    ax.set_yscale("log")
    ax.set_xticks(range(len(lanes)), list(lanes))
    ax.set_ylabel("Fragment length (bp)")
    ax.set_title("Electrophoresis Gel")

    plt.show()


# Read DNA
path = input("FASTA file path: ").strip().strip('"')
dna = read_fasta(path)[0]

print("DNA length:", len(dna))

# Choose enzymes
names = list(enzymes)

for i, name in enumerate(names):
    print(i+1, name)

n = int(input("How many enzymes (1-5)? "))
choices = list(map(int, input("Choose enzyme numbers: ").split()))

if not 1 <= n <= 5 or len(choices) != n or len(set(choices)) != n or any(c < 1 or c > 5 for c in choices):
    raise ValueError("Choose exactly n different enzymes from 1 to 5.")

lanes = {}
all_cuts = set()

# Restriction digestion
for choice in choices:
    name = names[choice-1]
    site, offset = enzymes[name]

    cuts = []

    for i in range(len(dna) - len(site) + 1):
        if dna[i:i+len(site)] == site:
            cut = i + offset
            if 0 < cut < len(dna):
                cuts.append(cut)

    all_cuts.update(cuts)
    lengths = fragments(dna, cuts)
    lanes[name] = lengths

    print("\n", name)
    print("Number of cuts:", len(cuts))
    print("Cut positions:", cuts)
    print("Fragment lengths:", lengths)

# Combined digestion
if n > 1:
    combined_cuts = sorted(all_cuts)
    combined_lengths = fragments(dna, combined_cuts)
    lanes["Combined"] = combined_lengths

    print("\nCombined digestion")
    print("Number of cuts:", len(combined_cuts))
    print("Cut positions:", combined_cuts)
    print("Fragment lengths:", combined_lengths)

# Electrophoresis gel
draw_gel(lanes)