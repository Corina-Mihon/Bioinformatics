import numpy as np
import matplotlib.pyplot as plt

s1 = "ACCGTGAAGCCAATAC"
s2 = "AGCGTGCAGCCAATAC"

match = 1
mismatch = -1
gap = -2

# Initialize matrix
matrix = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)

for i in range(len(s1) + 1):
    matrix[i, 0] = i * gap

for j in range(len(s2) + 1):
    matrix[0, j] = j * gap

# Fill matrix
for i in range(1, len(s1) + 1):
    for j in range(1, len(s2) + 1):

        diagonal = matrix[i-1, j-1] + (
            match if s1[i-1] == s2[j-1] else mismatch
        )
        up = matrix[i-1, j] + gap
        left = matrix[i, j-1] + gap

        matrix[i, j] = max(diagonal, up, left)

# Traceback
i, j = len(s1), len(s2)
a1, a2 = "", ""

while i > 0 or j > 0:
    if i > 0 and j > 0 and matrix[i, j] == matrix[i-1, j-1] + (
        match if s1[i-1] == s2[j-1] else mismatch
    ):
        a1 = s1[i-1] + a1
        a2 = s2[j-1] + a2
        i -= 1
        j -= 1

    elif i > 0 and matrix[i, j] == matrix[i-1, j] + gap:
        a1 = s1[i-1] + a1
        a2 = "-" + a2
        i -= 1

    else:
        a1 = "-" + a1
        a2 = s2[j-1] + a2
        j -= 1

# Plot matrix
fig, ax = plt.subplots(figsize=(12, 10))

ax.imshow(matrix, cmap="coolwarm")

ax.set_xticks(range(len(s2) + 1))
ax.set_yticks(range(len(s1) + 1))

ax.set_xticklabels(["-"] + list(s2))
ax.set_yticklabels(["-"] + list(s1))

ax.set_xlabel("Sequence 2")
ax.set_ylabel("Sequence 1")
ax.set_title("Needleman-Wunsch Scoring Matrix")

# Display scores inside cells
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        ax.text(j, i, str(matrix[i, j]),
                ha="center", va="center", fontsize=8)

ax.set_xticks(np.arange(-0.5, len(s2) + 1, 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(s1) + 1, 1), minor=True)
ax.grid(which="minor", color="black", linewidth=0.5)
ax.tick_params(which="minor", bottom=False, left=False)

plt.tight_layout()
plt.show()

print("Alignment:")
print(a1)
print("".join("|" if x == y else " " for x, y in zip(a1, a2)))
print(a2)
print("Score:", matrix[-1, -1])