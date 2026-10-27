# Project 06: Neighbor Joining 

## Description 
This project builds an unrooted phylogenetic tree, using neighbor joining, from a FASTA file of HIV-1 reverse transcriptase protein sequences. It reads the sequences, then runs the Smith-Waterman algorithm on each pair to get similarity scores, then converts those scores into a distance matrix. Neighbor joining uses that distance matrix to repeatedly join the best neighbor pair (the pair with the smallest `D*` score) until the full tree is built. It outputs a tree in Newick format and a plot of that tree. HIV lineages evolve at different rates due to rapid mutation and drug pressure, so neighbor joining is a better fit than UPGMA for this data because it does not assume all lineages evolve at the same rate. 

## Documentation and Usage
### Setup
Create and activate a conda environment with the required packages:
```bash
conda create -n module06 -c bioconda -c conda-forge python numpy matplotlib biopython ipykernel
conda activate module06
```

### Running the project
1. Open `project06/project06.ipynb` in VS Code or Jupyter and select the `module06` environment as the kernel
2. Make sure the sequence file is at `project06/data/lafayette_SARS_RT.fasta`. The notebook reads it using the relative path `data/lafayette_SARS_RT.fasta`.
3. Run all cells from top to bottom. The final cell reads the sequences, builds the distance matrix, runs neighbor joining, and plots the tree

### Functions
- **`read_fasta(filename)`**: Takes the path to a FASTA file and returns a dictionary mapping each sequence ID to its sequence. The ID is the last word of each header line since the full headers contain commas and brackets that would break the Newick format
- **`cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap)`**: Helper reused from Project 05. Returns the Smith-Waterman score for one cell of the scoring matrix (never below 0). The traceback direction was removed since only the score is needed
- **`smith_waterman(sequence1, sequence2, match_score=2.0, mismatch_penalty=-1.0, gap_penalty=-1.0)`**: Returns a normalized similarity score between 0 and 1: the best local alignment score divided by the average of the two sequences' self-scores (each sequence's length x match score).
- **`build_distance_matrix(sequences)`**: Takes the sequence dictionary and returns a tuple of the n x n distance matrix and the list of IDs in matrix order. Each similarity is clamped to [0.000001, 1] and converted to a distance with -ln
- **`neighbor_joining(distance_matrix, labels)`**: Takes the distance matrix and IDs and returns the tree as a Newick string. Negative branch lengths are set to 0, with the full pairwise distance given to the sister branch
- **`plot_tree(tree)`**: Provided with the template. Draws the Newick tree with Biopython.


### Output
The notebook diplays a plot of the tree. The plot shows a root on the left since `Bio.Phylo` draws every tree as if it had a root on the left, even though neighbor joining never chose one, so the left to right order of splits means nothing. The y-axis is labeled with numbers, but that's just how the plot spaces out the labels, and it has no meaning. Only which sequences group together, and the horizontal branch lengths should be interpreted.

## Pseudocode
### read_fasta
```
sequence = empty dict
current_id = None
for each line in file:
    strip whitespace
    if line starts with ">":
        current_id = the last word of the header line
        sequence[current_id] = empty string
    else if line is not empty:
        append line to sequence[current_id]
return sequence
```

### smith_waterman
``` 
raw local alignment score (adapted from project 05, traceback not required here)
m, n = lengths of sequence1, sequence2
create (m+1) x (n+1) score matrix of zeros (float)
highest_score = 0
for i from 1 to m:
    for j from 1 to n:
        score = cal_score(matrix, sequence1, sequence2, i, j, match_score, mismatch_penalty, gap_penalty)
        matrix[i][j] = score
        if score > highest_score:
            highest_score = score

self-scores: best possible score for each sequence 
every position matches itself so length * match_score
self1 = m * match_score
self2 = n * match_score

normalize against the average self-score of this pair
if either self-score is 0 (empty sequence)
    return 0.0
return highest_score / ((self1 + self2) / 2)
```                

### build_distance_matrix
```
ids come from the dictionary's keys
create a matrix of size nxn where n is the number of ids, filled with zeros
for i from 0 to n-1
    for j from i+1 to n-1 to only visit each pair exactly once
        seq_i = sequences[ids[i]]
        seq_j = sequences[ids[j]]
        score = smith_waterman(seq_i, seq_j)
        keep score between 0.000001 and 1 to avoid -ln(0) = infinity or negative distance scores
        distance = -ln(score)
        matrix[i][j] = distance
        matrix[j][i] = distance
return matrix, ids
```

### neighbor_joining
```
make a copy of the distance matrix
make a copy of the labels list
n = number of labels
while n > 2:
    for each row i:
        r[i] = add together all values in row i
    for each pair (i, j) where i < j:
            D*[i][j] = (n-2) * D[i][j] - r[i] - r[j]
    find the pair (i,j) with the smallest D*

    arm_i = (D[i][j] + (r[i] - r[j]) / (n-2)) / 2
    arm_j = D[i][j] - arm_i
    if arm_i < 0:
        arm_i = 0
        arm_j = D[i][j]
    if arm_j < 0:
        arm_j = 0
        arm_i = D[i][j]
    for each node x that is not i or j:
        distance from new junction to x = (D[i][x] + D[j][x] - D[i][j]) / 2
    new_label = f"({labels[i]}:{arm_i},{labels[j]}:{arm_j})"
    remove rows and columns i and j from the matrix
    add a new row and column for the junction (distances from the step above)
    remove labels[i] and labels[j] from the labels list
    add new_label to the labels list in the same position as the new row
    n = n-1
final = f"({labels[0]},{labels[1]}:{D[0][1]});"
return final 
```

### Formulas

Notation: `D[i][j]` is the distance between nodes i and j, `n` is the number of nodes currently in the matrix, and `r[i]` is the sum of row i (how far node i is from all other nodes).

**Neighbor-joining matrix D\*** (reconstructed from the Background clue; matches the Q-matrix in the project notebook's Neighbor-Joining Algorithm cell)

`D*[i][j] = (n - 2) * D[i][j] - r[i] - r[j]`

How I reconstructed it: two sequences can look close only because both are close to everything (short branches near the middle of the tree). The row sum `r[i]` measures how far node i is from everyone. Subtracting `r[i]` and `r[j]` removes that bias, so pairs that are close to each other *relative to how close they are to everyone else* get the lowest score. `D[i][j]` is a single distance while each row sum adds up many, so `D[i][j]` is multiplied by `(n - 2)` to put them on the same footing. Because `n` and the row sums change after every join, D\* is recalculated each round.

Comparison with Saitou & Nei (1987): the original paper chooses the pair whose join gives the smallest total branch length for the tree. Studier & Keppler (1988) showed this simplifies to the D\* formula above, so both choose the same pair. The difference is in how the criterion is written, not which pair it selects. [Check this against the paper before submitting.]

**Branch lengths** (source: project notebook, Neighbor-Joining Algorithm cell)
- `arm_i = (D[i][j] + (r[i] - r[j]) / (n - 2)) / 2`
- `arm_j = D[i][j] - arm_i`

Start from an even split of `D[i][j]`, then shift length toward whichever node is farther from everyone (larger row sum).

**Distance from the new node u to each remaining node x** (source: project notebook, Neighbor-Joining Algorithm cell)
- `D[u][x] = (D[i][x] + D[j][x] - D[i][j]) / 2`

The paths from i to x and from j to x both pass through u. Adding them counts `D[u][x]` twice plus both arms (`D[i][j]`), so subtract `D[i][j]` and halve.

**Similarity to distance**
- `similarity = SW_score(a, b) / ((self_a + self_b) / 2)`, where `self_a = len(a) * match_score`
- similarity is clamped to [0.000001, 1], then `distance = -ln(similarity)`

Normalizing by the pair's average self-score corrects for sequence length (based on the normalization idea in Feng & Doolittle, 1987).

### Design decisions

**Tree representation:** the tree is stored as a list of labels that stays in the same order as the matrix rows. Each label starts as a sequence ID. When two nodes are joined, the new node's label is a Newick piece that holds both children and their branch lengths, e.g. `(A:1.0,C:10.0)`. That label is everything the new node needs to remember about the two nodes it joined.

**Iterative vs. recursive:** Newick is naturally recursive, since each internal node contains subtrees. I built it iteratively instead: each join wraps the two labels in parentheses, so the nesting builds up as the loop runs. The alternative would be to store a tree of node objects and write the Newick string recursively at the end.

**Ties:** if more than one pair has the same smallest D\*, the first one found in the loop is used. When three nodes remain, all D\* values are always tied, and any choice gives the same unrooted tree.

**Stopping:** the loop runs while more than 2 nodes remain, so it does n − 2 joins. The last two nodes are connected with one branch of length `D[0][1]`.

**Negative branch lengths:** these come from the method, not a bug. Real sequence distances are not perfectly additive, so the branch-length formula can overshoot. When an arm is negative, it is set to 0 and the other arm gets the full `D[i][j]`, so the distance between the two joined nodes is kept. On the HIV RT data, this happened 2 times out of 18 joins.

**Testing:** before using the HIV data, I ran `neighbor_joining` on a 4-taxon matrix with a known additive tree (A and C on one side, B and D on the other, with outer branches of 1 and 10 and a middle branch of 1). Raw distance alone would wrongly pair A with B. The output was `((A:1.0,C:10.0),(B:1.0,D:10.0):1.0);`, which recovers the correct pairs and every branch length exactly, with no negative branches. I also checked that the final tree for 20 sequences has 19 commas (one per join, plus the final connection).


## Sucesses

## Struggles

## Reflections 

## Appendix (AI usage, citations, etc.)
