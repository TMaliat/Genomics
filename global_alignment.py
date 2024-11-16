import numpy as np

# Parameters
match = -100
mismatch = -95
gap = -5

#mismatch = -5
#gap = -2
#
# Sequences 
seqB = "CTCGCAGC"
seqA = "CAGGCAGT"
  
  

# Initialize score matrix
rows = len(seqA) + 1
cols = len(seqB) + 1
scoring_matrix = np.zeros((rows, cols), dtype=int)

# Fill the first row and column with gap penalties
for i in range(1, rows):
    scoring_matrix[i][0] = scoring_matrix[i-1][0] + gap
for j in range(1, cols):
    scoring_matrix[0][j] = scoring_matrix[0][j-1] + gap

# Populate scoring matrix
for i in range(1, rows):
    for j in range(1, cols):
        if seqA[i-1] == seqB[j-1]:
            score_diag = scoring_matrix[i-1][j-1] + match
        else:
            score_diag = scoring_matrix[i-1][j-1] + mismatch
        score_up = scoring_matrix[i-1][j] + gap
        score_left = scoring_matrix[i][j-1] + gap
        scoring_matrix[i][j] = max(score_diag, score_up, score_left)

# Print score matrix
print("Score Matrix:")
print(scoring_matrix)

# Backtracking to find the aligned sequences and traceback values
alignA, alignB, match_line = "", "", ""
traceback_values = []  
i, j = len(seqA), len(seqB)

while i > 0 or j > 0:
    current_score = scoring_matrix[i][j]
    traceback_values.append(current_score)  # Add the current score to the traceback list
    if i > 0 and j > 0 and ((seqA[i-1] == seqB[j-1] and current_score == scoring_matrix[i-1][j-1] + match) or
                            (seqA[i-1] != seqB[j-1] and current_score == scoring_matrix[i-1][j-1] + mismatch)):
        alignA = seqA[i-1] + alignA
        alignB = seqB[j-1] + alignB
        match_line = ("|" if seqA[i-1] == seqB[j-1] else " ") + match_line
        i -= 1
        j -= 1
    elif i > 0 and current_score == scoring_matrix[i-1][j] + gap:
        alignA = seqA[i-1] + alignA
        alignB = "-" + alignB
        match_line = " " + match_line
        i -= 1
    else:
        alignA = "-" + alignA
        alignB = seqB[j-1] + alignB
        match_line = " " + match_line
        j -= 1

# Calculate alignment score
alignment_score = scoring_matrix[len(seqA)][len(seqB)]

# Print aligned sequences and score in the desired format
print("\nAligned Sequences:")
print("Sequence 1:", alignA)
print("           ", match_line)
print("Sequence 2:", alignB)
print("\nAlignment Score:", alignment_score)

# Print traceback path (from bottom-right to top-left)
print("\nTraceback Path (bottom-right to top-left):")
print(" -> ".join(map(str, traceback_values)))
