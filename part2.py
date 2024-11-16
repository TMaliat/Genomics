from collections import Counter

# Provided DNA sequences
sequences = [
    "actgatcgacgctatctggctatccaGgtacTtaggtcctctgtgcgaatctatgcgtttccaaccat",
    "agtactggtgtccgtatagtCcAtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc",
    "aaacgtTAgtgcaccctctttcttcgtggctctcctgtaagagggctgatgtataagacgaaaatttt",
    "agcctccgatgtaagtcatagctgtccgtgataccctgccacccctattacatcttacgtCcAtataca",
    "ctgttatacaaccctagtaggcggggtatgcgttttggtcgtcgtacgctcgatcgttaCcgtacgGc"
]

l = 8  # Length of l-mers to find

# Function to compute the match count between two strings (instead of Hamming distance)
def match_count(str1, str2):
    return sum(1 for x, y in zip(str1, str2) if x == y)

# Function to extract all possible l-mers from a sequence
def get_l_mers(sequence, l):
    return [sequence[i:i + l] for i in range(len(sequence) - l + 1)]

# Function to compute the consensus string
def get_consensus(group):
    # Find consensus string by getting the most common base at each position
    return "".join(Counter(column).most_common(1)[0][0] for column in zip(*group))

# Open the file to write iteration outputs
with open('one.txt', 'w') as f:
    # Step 1: Initial Guessing
    groups = []
    f.write("\nStep 1: Initial Guessing\n")

    # Iterate over all l-mers in the first sequence
    for i in range(len(sequences[0]) - l + 1):  # Iterate over all l-mers in the first sequence
        initial_l_mer = sequences[0][i:i + l]
        group = [initial_l_mer]
        f.write(f"\nStarting with initial l-mer from Seq 1, position {i}: {initial_l_mer}\n")

        # Find the closest l-mer in the second sequence
        closest_l_mer = initial_l_mer
        total_match_count = 0  # Sum of match counts for the group
        
        for j in range(1, len(sequences)):  # Find closest l-mers from remaining sequences
            scores = [match_count(closest_l_mer, candidate) for candidate in get_l_mers(sequences[j], l)]
            min_index = scores.index(max(scores))  # Max match count
            closest_l_mer = get_l_mers(sequences[j], l)[min_index]
            group.append(closest_l_mer)
            match_cnt = max(scores)
            total_match_count += match_cnt
            f.write(f"  Closest l-mer from Seq {j+1}: {closest_l_mer} with match count {match_cnt}\n")
        
        # Get consensus motif for the group
        consensus_motif = get_consensus(group)
        
        # After calculating the total match counts and consensus motif, print the result
        f.write(f"  Formed group: {group}, Consensus motif: {consensus_motif}, Total match count: {total_match_count}\n")
        groups.append((group, consensus_motif, total_match_count))

    # Step 2: Find the group with the highest match count and calculate the new consensus
    highest_score_group = max(groups, key=lambda x: x[2])
    highest_consensus_motif = highest_score_group[1]
    f.write(f"\nHighest score group: Consensus motif: {highest_consensus_motif}, Score: {highest_score_group[2]}\n")

    # Step 3: Iteratively update consensus until no further improvement is possible
    previous_score = highest_score_group[2]
    current_consensus_motif = highest_consensus_motif
    iteration = 0

    # Continue until no improvement in the score can be achieved
    while True:
        iteration += 1
        f.write(f"\nIteration {iteration} - Current consensus motif: {current_consensus_motif}\n")
        
        updated_groups = []
        new_score = 0
        
        for i, seq in enumerate(sequences):
            f.write(f"Sequence {i+1}: {seq}\n")
            distances = []
            
            # Compare with all l-mers in the sequence
            for l_mer in get_l_mers(seq, l):
                cnt = match_count(current_consensus_motif, l_mer)
                distances.append((l_mer, cnt))
                f.write(f"  l-mer: {l_mer}, Match count with consensus: {cnt}\n")
            
            # Now, find the best l-mer considering the match count
            best_l_mer = max(distances, key=lambda x: x[1])  # Max match count
            updated_groups.append(best_l_mer[0])
            new_score += best_l_mer[1]
            f.write(f"  Best l-mer from Seq {i+1}: {best_l_mer[0]} with match count {best_l_mer[1]}\n")
        
        # Recalculate the new consensus motif based on updated l-mers, considering all matches
        non_zero_l_mers = [l_mer for l_mer in updated_groups if match_count(current_consensus_motif, l_mer) > 0]
        new_consensus_motif = get_consensus(non_zero_l_mers)
        f.write(f"New consensus motif: {new_consensus_motif}, Total match count: {new_score}\n")
        
        # If the new score is not better than the previous score, stop the iteration
        if new_score <= previous_score:
            f.write("No further improvement in score. Stopping iteration.\n")
            f.write(f"Final consensus motif: {current_consensus_motif}\n")
            f.write(f"Best score: {previous_score}\n")
            break

        # Update the current consensus motif and previous score for the next iteration
        current_consensus_motif = new_consensus_motif
        previous_score = new_score

# Print the best consensus and its score in the terminal
print(f"\nFinal Consensus Motif: {current_consensus_motif}")
print(f"Best Score (based on matches): {previous_score}")
