
from collections import Counter
from google import genai
import gzip
import json
import os

file_path = r"GCF_002194405.1_ViralProj390412_genomic.fna.gz"
MODEL = "gemini-3.5-flash-lite"

def read_genome(path):
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


def count_codons(sequences):
    counts = Counter()

    for sequence in sequences:
        for i in range(len(sequence) - 2):
            codon = sequence[i:i + 3]

            if all(base in "ACGT" for base in codon):
                counts[codon] += 1

    return counts


def get_top_10(counts):
    total = sum(counts.values())

    if total == 0:
        raise ValueError("No valid triplets found.")

    top_10 = []

    print("\nTOP 10 MOST FREQUENT CODONS")
    print("-" * 45)

    for codon, count in counts.most_common(10):
        frequency = round(count / total * 100, 2)

        print(f"{codon}: {count} ({frequency:.2f}%)")

        top_10.append({
            "codon": codon,
            "frequency": frequency
        })

    return top_10

def analyze_with_gemini(top_10):
    client = genai.Client()

    prompt = """
You are an AI assistant for a bioinformatics university project.

You receive the top 10 most frequent nucleotide
triplets found in a viral genome.

Your task is to identify common foods derived from
organisms whose DNA is likely to contain these triplets.

Instructions:
1. Analyze the provided top 10 triplets.
2. Suggest 10 common foods from plants, animals,
   or fungi whose genomes are likely to contain
   these nucleotide triplets.
3. For each food, explain briefly why its source
   organism is expected to contain these triplets.
4. Do not claim that these foods are dangerous.
5. Clearly state that the results are theoretical
   examples, not verified genomic matches.
6. Explain that sharing triplets with a virus
   does not make a food unsafe to eat.

Respond using this format:

FOODS WITH POTENTIALLY MATCHING CODONS:

1. Food name - explanation
2. Food name - explanation
...

TOP 10 CODONS ANALYZED:
...

IMPORTANT BIOLOGICAL NOTE:
...

Here are the top 10 triplets:
"""

    prompt += json.dumps(top_10, indent=2)

    response = client.models.generate_content(
        model=MODEL,
        contents=prompt
    )

    return response.text


def main():
    print("VIRAL GENOME ANALYZER")
    print("=" * 45)

    if not os.path.isfile(file_path):
        print("File not found:", file_path)
        return

    print("\nReading genome...")

    sequences = read_genome(file_path)

    if not sequences:
        print("No DNA sequences found.")
        return

    print("Sequences:", len(sequences))
    print("Nucleotides:",
          sum(len(s) for s in sequences))

    counts = count_codons(sequences)
    top_10 = get_top_10(counts)

    print("\nConnecting to Gemini...")

    try:
        result = analyze_with_gemini(top_10)

        print("\nAI DIETARY ANALYSIS")
        print("=" * 45)
        print(result)

    except Exception as error:
        print("Gemini API error:", error)


if __name__ == "__main__":
    main()
