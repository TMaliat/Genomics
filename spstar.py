from collections import Counter

# Provided DNA sequences
sequences = [
    "cctgatagacgctatctggctatccaGgtacTtaggtcctctgtgcgaatctatgcgtttccaaccat",
"agtactggtgtacatttgatCcAtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc",
"aaacgtTAgtgcaccctctttcttcgtggctctggccaacgagggctgatgtataagacgaaaatttt",
"agcctccgatgtaagtcatagctgtaactattacctgccacccctattacatcttacgtCcAtataca",
"ctgttatacaacgcgtcatggcggggtatgcgttttggtcgtcgtacgctcgatcgttaCcgtacgGc"
]

l = 8  # Length of l-mers to find

# Function to compute the Hamming distance between two strings
def hamming_distance(str1, str2):
    return sum(1 for x, y in zip(str1, str2) if x != y)

# Function to extract all possible l-mers from a sequence
def get_l_mers(sequence, l):
    return [sequence[i:i + l] for i in range(len(sequence) - l + 1)]

# Function to compute the consensus string
def get_consensus(group):
    # Find consensus string by getting the most common base at each position
    return "".join(Counter(column).most_common(1)[0][0] for column in zip(*group))

# Step 1: Initial Guessing
groups = []
print("\nStep 1: Initial Guessing")

# Iterate over all l-mers in the first sequence
for i in range(len(sequences[0]) - l + 1):  # Iterate over all l-mers in the first sequence
    initial_l_mer = sequences[0][i:i + l]
    group = [initial_l_mer]
    print(f"\nStarting with initial l-mer from Seq 1, position {i}: {initial_l_mer}")

    # Find the closest l-mer in the second sequence
    closest_l_mer = initial_l_mer
    total_hamming_distance = 0  # Sum of Hamming distances for the group
    
    for j in range(1, len(sequences)):  # Find closest l-mers from remaining sequences
        scores = [hamming_distance(closest_l_mer, candidate) for candidate in get_l_mers(sequences[j], l)]
        min_index = scores.index(min(scores))
        closest_l_mer = get_l_mers(sequences[j], l)[min_index]
        group.append(closest_l_mer)
        hamming_dist = min(scores)
        total_hamming_distance += hamming_dist
        print(f"  Closest l-mer from Seq {j+1}: {closest_l_mer} with Hamming distance {hamming_dist}")
    
    # Get consensus motif for the group
    consensus_motif = get_consensus(group)
    
    # After calculating the total Hamming distances and consensus motif, print the result
    print(f"  Formed group: {group}, Consensus motif: {consensus_motif}, Total Hamming distance: {total_hamming_distance}")
    groups.append((group, consensus_motif, total_hamming_distance))

# Step 2: Find the group with the lowest Hamming distance score and calculate the new consensus
lowest_score_group = min(groups, key=lambda x: x[2])
lowest_consensus_motif = lowest_score_group[1]
print(f"\nLowest score group: Consensus motif: {lowest_consensus_motif}, Score: {lowest_score_group[2]}")

# Step 3: Iteratively update consensus until no further improvement is possible
previous_score = lowest_score_group[2]
current_consensus_motif = lowest_consensus_motif
iteration = 0

# Continue until no improvement in the score can be achieved
while True:
    iteration += 1
    print(f"\nIteration {iteration} - Current consensus motif: {current_consensus_motif}")
    
    updated_groups = []
    new_score = 0
    
    for i, seq in enumerate(sequences):
        print(f"Sequence {i+1}: {seq}")
        distances = []
        
        # Compare with all l-mers in the sequence
        for l_mer in get_l_mers(seq, l):
            dist = hamming_distance(current_consensus_motif, l_mer)
            distances.append((l_mer, dist))
            print(f"  l-mer: {l_mer}, Hamming distance with consensus: {dist}")
        
        # Now, find the best l-mer considering the Hamming distance, including zero distances
        # The best l-mer is still considered for the internal selection, even if it has distance 0
        best_l_mer = min(distances, key=lambda x: x[1])  # Consider all, including distance 0
        updated_groups.append(best_l_mer[0])
        new_score += best_l_mer[1]
        print(f"  Best l-mer from Seq {i+1}: {best_l_mer[0]} with Hamming distance {best_l_mer[1]}")
    
    # Recalculate the new consensus motif based on updated l-mers, but exclude zero-distance l-mers
    non_zero_l_mers = [l_mer for l_mer in updated_groups if hamming_distance(current_consensus_motif, l_mer) > 0]
    new_consensus_motif = get_consensus(non_zero_l_mers)
    print(f"New consensus motif: {new_consensus_motif}, Total Hamming distance: {new_score}")
    
    # If the new score is not better than the previous score, stop the iteration
    if new_score >= previous_score:
        print("No further improvement in score. Stopping iteration.")
        print(f"Final consensus motif: {current_consensus_motif}")
        print(f"Best score: {previous_score}")
        break

    # Update the current consensus motif and previous score for the next iteration
    current_consensus_motif = new_consensus_motif
    previous_score = new_score
