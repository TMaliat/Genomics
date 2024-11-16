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
    return "".join(Counter(column).most_common(1)[0][0] for column in zip(*group))

# Set to store valid d values
valid_d_values = set()

# Open the file for writing the iteration-based outputs
with open('three.txt', 'w') as f:
    # Iterate over all possible values of d (0 to l)
    for d in range(l + 1):
        f.write(f"\nTrying d = {d} (max allowed Hamming distance)\n")

        groups = []
        previous_score = float('inf')  # Set initial score to infinity
        iteration = 0
        consensus_motif = None
        
        # Step 1: Initial Guessing for d = current mutation threshold
        for i in range(len(sequences[0]) - l + 1):
            initial_l_mer = sequences[0][i:i + l]
            group = [initial_l_mer]

            # Find the closest l-mer in the other sequences
            for j in range(1, len(sequences)):
                scores = [hamming_distance(initial_l_mer, candidate) for candidate in get_l_mers(sequences[j], l)]
                
                # Filter l-mers that have a Hamming distance ≤ d
                valid_l_mers = [get_l_mers(sequences[j], l)[index] for index, score in enumerate(scores) if score <= d]
                
                if valid_l_mers:
                    closest_l_mer = valid_l_mers[0]  # Choose the first valid l-mer
                    group.append(closest_l_mer)
            
            # Compute consensus motif for the group
            consensus_motif = get_consensus(group)
            groups.append((group, consensus_motif))

        # Step 2: Refining the Consensus Motif for d = current mutation threshold
        while True:
            iteration += 1
            f.write(f"\nIteration {iteration} - Refining consensus motif\n")  # Write to the file

            new_groups = []
            new_score = 0
            
            # Recalculate groups based on closest l-mers with Hamming distance ≤ d
            for i, seq in enumerate(sequences):
                f.write(f"Sequence {i+1}: {seq}\n")  # Write sequence info to the file
                distances = []
                
                # Compare with all l-mers in the sequence
                for l_mer in get_l_mers(seq, l):
                    dist = hamming_distance(consensus_motif, l_mer)
                    distances.append((l_mer, dist))
                    f.write(f"  l-mer: {l_mer}, Hamming distance with consensus: {dist}\n")  # Write to the file
                
                # Now, find the best l-mer considering the Hamming distance, including zero distances
                valid_distances = [(l_mer, dist) for l_mer, dist in distances if dist <= d]
                
                # If there are valid l-mers, choose the best one (lowest Hamming distance)
                if valid_distances:
                    best_l_mer = min(valid_distances, key=lambda x: x[1])  # Consider the l-mer with the lowest distance
                    new_groups.append(best_l_mer[0])
                    new_score += best_l_mer[1]
                    f.write(f"  Best l-mer from Seq {i+1}: {best_l_mer[0]} with Hamming distance {best_l_mer[1]}\n")  # Write to the file
            
            # Recalculate the new consensus motif based on updated l-mers
            non_zero_l_mers = [l_mer for l_mer in new_groups if hamming_distance(consensus_motif, l_mer) > 0]
            new_consensus_motif = get_consensus(non_zero_l_mers)
            f.write(f"New consensus motif: {new_consensus_motif}, Total Hamming distance: {new_score}\n")  # Write to the file
            
            # If the new score is not better than the previous score, stop the iteration
            if new_score >= previous_score:
                f.write("No further improvement in score. Stopping iteration.\n")  # Write to the file
                break

            # Update the current consensus motif and previous score for the next iteration
            consensus_motif = new_consensus_motif
            previous_score = new_score
        
        # Check if a valid consensus motif was found
        if consensus_motif and previous_score < float('inf'):
            valid_d_values.add(d)
            print(f"Valid motif found for d = {d}: {consensus_motif}, Score: {previous_score}")

# Print the highest valid d value with its motif and score
if valid_d_values:
    highest_d = max(valid_d_values)
    print(f"\nHighest valid d value: {highest_d}")
    print(f"Motif: {consensus_motif}, Score: {previous_score}")
else:
    print("\nNo valid motifs found for any d value.")
