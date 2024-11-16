import numpy as np

# Parameters for scoring
match_score = 10
mismatch_penalty = -2
gap_penalty = -5

# Sequences
seqA = "ACACACTA"
seqB = "AGCACACA"

# Initialize the scoring matrix with zeros
rows = len(seqA) + 1
cols = len(seqB) + 1
scoring_matrix = np.zeros((rows, cols), dtype=int)

# Variables to track the maximum score and its position
max_score = 0
max_pos = (0, 0)

# Fill the scoring matrix
for i in range(1, rows):
    for j in range(1, cols):
        if seqA[i-1] == seqB[j-1]:
            score_diag = scoring_matrix[i-1][j-1] + match_score
        else:
            score_diag = scoring_matrix[i-1][j-1] + mismatch_penalty
        score_up = scoring_matrix[i-1][j] + gap_penalty
        score_left = scoring_matrix[i][j-1] + gap_penalty
        scoring_matrix[i][j] = max(0, score_diag, score_up, score_left)  # Local alignment includes zero as a minimum

        # Track the maximum score in the matrix
        if scoring_matrix[i][j] > max_score:
            max_score = scoring_matrix[i][j]
            max_pos = (i, j)

# Print score matrix
print("Score Matrix:")
print(scoring_matrix)

# Backtracking to find the optimal local alignment
i, j = max_pos
longest_match_line = ""

# Keep track of the longest contiguous match substring
temp_match_line = ""
while i > 0 and j > 0 and scoring_matrix[i][j] != 0:
    current_score = scoring_matrix[i][j]
    if i > 0 and j > 0 and ((seqA[i-1] == seqB[j-1] and current_score == scoring_matrix[i-1][j-1] + match_score) or
                            (seqA[i-1] != seqB[j-1] and current_score == scoring_matrix[i-1][j-1] + mismatch_penalty)):
        if seqA[i-1] == seqB[j-1]:  # Only add to the match line if there's an exact match
            temp_match_line = seqA[i-1] + temp_match_line
        else:
            # If there's a mismatch, finalize the current contiguous match substring
            if len(temp_match_line) > len(longest_match_line):
                longest_match_line = temp_match_line
            temp_match_line = ""
        i -= 1
        j -= 1
    elif i > 0 and current_score == scoring_matrix[i-1][j] + gap_penalty:
        # Finalize the current contiguous match substring if a gap is encountered
        if len(temp_match_line) > len(longest_match_line):
            longest_match_line = temp_match_line
        temp_match_line = ""
        i -= 1
    else:
        # Finalize the current contiguous match substring if a gap is encountered
        if len(temp_match_line) > len(longest_match_line):
            longest_match_line = temp_match_line
        temp_match_line = ""
        j -= 1

# Check at the end of backtracking if the last temp_match_line is the longest
if len(temp_match_line) > len(longest_match_line):
    longest_match_line = temp_match_line

# Print the longest contiguous matching substring without gaps
print("\nLongest contiguous matching substring without gaps:", longest_match_line)
