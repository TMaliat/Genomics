import matplotlib.pyplot as plt
import networkx as nx
import random

# Sample sequences
sequences = [
    "GTGCTGCACGGCTCAGTATAGCATTTACCCTTCCATCTTCAGATCCTGAA",
    "ACGCTGCACGGCTCAGTGCGGTGCTTACCCTCCCATCTTCAGATCCTGAA",
    "GTGCTGCACGGCTCGGCGCAGCATTTACCCTCCCATCTTCAGATCCTATC",
    "GTATCACACGACTCAGCGCAGCATTTGCCCTCCCGTCTTCAGATCCTAAA",
    "GTATCACATAGCTCAGCGCAGCATTTGCCCTCCCGTCTTCAGATCTAAAA"
]
# Hamming distance function
def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Format for joined clusters
def strFormat(i, j):
    return f"({i} {j})"

# Initialize distance matrix
n = len(sequences)
distance = {f'S{i+1}': {} for i in range(n)}
for i in range(n):
    for j in range(n):
        distance[f'S{i+1}'][f'S{j+1}'] = hamming(sequences[i], sequences[j]) if i != j else 0

# UPGMA clustering loop
joins = 0
adjacency_list = {}
while joins < n - 1:
    min_dist = float('inf')
    key1, key2 = "", ""

    # Find the minimum distance in the matrix
    for row in distance.keys():
        for col in distance[row]:
            if distance[row][col] is not None and distance[row][col] < min_dist and row != col:
                key1, key2 = row, col
                min_dist = distance[row][col]

    # Print current iteration and joining step
    print(f"Iteration {joins + 1}:")
    print(f"Joining {key1} and {key2} with distance {min_dist}")

    # Create new key for the merged cluster
    new_key = strFormat(key1, key2)
    distance[new_key] = {}

    # Update the distance matrix with the new merged cluster
    for row in list(distance.keys()):
        if row != key1 and row != key2:
            val1 = distance[key1].get(row, None)
            val2 = distance[key2].get(row, None)
            if val1 is not None and val2 is not None:
                avg_dist = 0.5 * (val1 + val2)
                distance[new_key][row] = avg_dist
                distance[row][new_key] = avg_dist

    # Remove old clusters from the distance matrix
    del distance[key1]
    del distance[key2]
    for row in distance.keys():
        distance[row].pop(key1, None)
        distance[row].pop(key2, None)

    # Print the updated distance matrix
    print("Updated Distance Matrix:")
    for row, dist_vals in distance.items():
        print(f"{row}: {dist_vals}")
    print("\n")

    # Store adjacency for visualization
    adjacency_list[new_key] = [key1, key2]
    joins += 1

# Final root node for the tree
root_key = f"Root: {new_key}"
adjacency_list[root_key] = adjacency_list.pop(new_key)

# Print the final adjacency list
print("\nAdjacency List for Tree:")
print(adjacency_list)

# Visualization code
def hierarchy_pos(G, root=None, width=1., vert_gap=0.3, vert_loc=0, xcenter=0.5):
    """Generate positions for a hierarchical tree."""
    if not nx.is_tree(G):
        raise TypeError('Cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        root = next(iter(nx.topological_sort(G)))

    def _hierarchy_pos(G, root, width=1., vert_gap=0.3, vert_loc=0, xcenter=0.5, pos=None, parent=None):
        if pos is None:
            pos = {root: (xcenter, vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)
        if children:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G, child, width=dx, vert_gap=vert_gap,
                                     vert_loc=vert_loc + vert_gap, xcenter=nextx,
                                     pos=pos, parent=root)
        return pos

    return _hierarchy_pos(G, root, width, vert_gap, vert_loc)

# Draw the final tree with root at the bottom
def draw_tree(adjacency_list, root):
    G = nx.Graph()
    for node in adjacency_list:
        G.add_node(node)
    for node, neighbors in adjacency_list.items():
        for neighbor in neighbors:
            G.add_edge(node, neighbor)

    plt.figure(figsize=(10, 8))
    pos = hierarchy_pos(G, root, vert_loc=1, vert_gap=0.15)  # Root at bottom
    nx.draw(G, pos, with_labels=True, node_size=1100, node_color='#6f5193', font_size=10, font_color='black', font_weight='bold', width=1.2)
    plt.title('Phylogenetic Tree with Root at Bottom (Upside Down)')
    plt.show()

# Call draw function with the final root
draw_tree(adjacency_list, root_key)
