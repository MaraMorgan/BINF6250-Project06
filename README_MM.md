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
1. Open `project06.ipynb` in VS Code or Jupyter and select the `module06` environment as the kernel
2. Download the data from the Canvas assignment page and place it at `data/lafayette_SARS_RT.fasta` next to the notebook. The data is not included in this repository.
3. Run all cells from top to bottom. The final cell reads the sequences, builds the distance matrix, runs neighbor joining, and plots the tree


#### Dependencies: 
python 3.12
##### imports
- typing: List, Tuple, Dict
- io: StringIO
- numpy
- matplotlib.pyplot
- Bio: Phylo, SeqIO


### Functions
- **`read_fasta(filename)`**: Takes the path to a FASTA file and returns a dictionary mapping each sequence ID to its sequence. The sequence ID is the last word on each header line. 
- **`cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap)`**: Helper reused from Project 05. Returns the Smith-Waterman score for one cell of the scoring matrix (never below 0). The traceback direction was removed since only the score is needed
- **`smith_waterman(sequence1, sequence2, match_score=2.0, mismatch_penalty=-1.0, gap_penalty=-1.0)`**: Returns a normalized similarity score between 0 and 1: the best local alignment score divided by the average of the two sequences' self-scores (each sequence's length x match score). Because all sequences in this dataset are the same length, normalizing it by the average self-score gives the same result as other options such as the geometric mean (the square root of the product of the self-scores)

- **`build_distance_matrix(sequences)`**: Takes the sequence dictionary and returns a tuple of the n x n distance matrix and the list of IDs in matrix order. Each similarity is clamped to [0.000001, 1] and converted to a distance with -ln
- **`neighbor_joining(distance_matrix, labels)`**: Takes the distance matrix and IDs and returns the tree as a Newick string. Negative branch lengths are set to 0, with the full pairwise distance given to the sister branch
- **`plot_tree(tree)`**: Provided with the template. Draws the Newick tree with Biopython.


### Output
The notebook displays a plot of the tree. The plot shows a root on the left since `Bio.Phylo` draws every tree as if it had a root on the left, even though neighbor joining never chose one, so the left to right order of splits means nothing. The y-axis is labeled with numbers, but that's just how the plot spaces out the labels, and it has no meaning. Only which sequences group together, and the horizontal branch lengths should be interpreted.

## Pseudocode
### read_fasta
```
parse the FASTA file with Biopython's SeqIO
for each record:
    id = last word of the header line (the short ID)
    sequence = the record's sequence as a plain string
return dictionary of {id: sequence}
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

Comparison with Saitou & Nei (1987): Saitou & Nei describe the method as choosing, at each step, the pair of neighbors that minimizes the total branch length of the tree, starting from a star-shaped tree. The `D*` formula used here is the formulation given by Studier & Keppler (1988), who present it as an alternative to Saitou & Nei's formulation that runs in O(N^3) time (equations 1–2). Studier & Keppler also note that their distance-update formula differs from the one in Saitou & Nei and is the correct one when the distances are treelike (equation 3), and that Saitou & Nei's proof that the method recovers the correct tree from treelike data was incorrect, which they correct.

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

**Why -ln** -ln flips the scale so identical sequences have distance 0 and less similar sequences have larger distances. As sequences diverge, similarity decreases fractionally, not by a fixed amount, so log makes distances grow roughly in proportion to evolutionary change. This keeps distances closer to additive along a tree. For HIV RT sequences in this project, which are very similar, -ln(sim) and 1 - sim give nearly the same values; this choice matters more for distantly related sequences.

Normalizing by the pair's average self-score corrects for sequence length (based on the normalization idea in Feng & Doolittle, 1987).

### Design decisions

**Tree representation:** the tree is stored as a list of labels that stays in the same order as the matrix rows. Each label starts as a sequence ID. When two nodes are joined, the new node's label is a Newick piece that holds both children and their branch lengths, e.g. `(A:1.0,C:10.0)`. That label is everything the new node needs to remember about the two nodes it joined.

**Iterative vs. recursive:** Newick is naturally recursive, since each internal node contains subtrees. I built it iteratively instead: each join wraps the two labels in parentheses, so the nesting builds up as the loop runs. The alternative would be to store a tree of node objects and write the Newick string recursively at the end.

**Ties:** if more than one pair has the same smallest D\*, the first one found in the loop is used. When three nodes remain, all D\* values are always tied, and any choice gives the same unrooted tree.

**Stopping:** the loop runs while more than 2 nodes remain, so it does n − 2 joins. The last two nodes are connected with one branch of length `D[0][1]`.

**Negative branch lengths:** these come from the method, not a bug. Real sequence distances are not perfectly additive, so the branch-length formula can overshoot. When an arm is negative, it is set to 0 and the other arm gets the full `D[i][j]`, so the distance between the two joined nodes is kept. On the HIV RT data, this happened 2 times out of 18 joins.

**Testing:** before using the HIV data, I ran `neighbor_joining` on a 4-taxon matrix with a known additive tree (A and C on one side, B and D on the other, with outer branches of 1 and 10 and a middle branch of 1). Raw distance alone would wrongly pair A with B. The output was `((A:1.0,C:10.0),(B:1.0,D:10.0):1.0);`, which recovers the correct pairs and every branch length exactly, with no negative branches. I also checked that the final tree for 20 sequences has 19 commas (one per join, plus the final connection). The docstring examples provided with the template do not match our results exactly. They list 19 sequences, while the provided file contains 20, so our Newick string has 19 commas instead of 18. The example distances (0.004 and 0.012) also differ from ours (0.007, 0.047, and 0.040 for the first three sequences). Because all the sequences are the same length, the choice of normalization does not affect the distances, so the examples were most likely generated from a different version of the data. 


## Successes
- One of our biggest successes was how we worked together to finalize the project. We all approached certain parts differently, but when we compared our code, we found that a lot of it was actually very similar. We ended up combining parts of different implementations based on what made the most sense rather than just choosing one person's work.
- Being able to build on our Smith-Waterman code from Project 05 saved us a lot of time. We already understood the alignment process, so we could spend more time figuring out how to normalize the scores, build the distance matrix, and use it for Neighbor-Joining.
- Our final group meeting was also really productive. We went through the different implementations, compared our trees, and talked through decisions like which distance transformation to use and how to handle negative branch lengths. We were able to explain our reasoning, question each other's choices, and agree on a final version without having to start everything over.

## Struggles
- One of our biggest struggles was figuring out which formulas to use, especially for normalizing the Smith-Waterman scores and converting similarity into distance. There were different ways to approach both, and it wasn't always clear which would be most appropriate for our sequences. We spent quite a bit of time comparing the formulas and discussing how they would affect the final distance matrix and tree.
- Neighbor-Joining was probably the most challenging part of the project. There were a lot of calculations to keep track of, and understanding how the distance matrix changes after joining two nodes took some time. We also had to figure out how to handle negative branch lengths and when to stop the algorithm, since different implementations handled these differently.
- Finding a time when everyone could meet for long enough continues to be a challenge, especially with everyone's different schedules. I think we did a great job making it work though, whether that meant working on things individually beforehand or making the most of the time we did have together.

## Reflections 
### Group Leader: Katelyn
Going into this project, I had a basic understanding of how the neighbor-joining algorithm works, but the Canvas module was a very helpful step-by-step that allowed me to fully refresh and clarify all the parts. Working out an example by hand helped me a lot to understand how to go about writing the pseudocode. I don't have as much experience using NumPy arrays, so the hardest part of this for me was understanding how to implement the pseudocode in python. I think that starting with a quick refresher of the NumPy functions themselves before I start would help me work more efficiently.
Our group worked well together, and we decided to independently attempt to write the pseudocode before meeting to discuss it. This was successful because it allowed us to all identify the parts that we didn't understand and all get the chance to try it. 
In the future, I think it would be interesting to consider that our Smith-Waterman scoring gives every mismatch the same penalty. For amino acids, a substitution matrix like BLOSUM62 would be more biologically accurate, since substitutions between chemically similar amino acids are more common and should be penalized less. This could be an interesting addition if this project were to be extended. 

### Other Members: Aamna
I think one of the biggest things I took away from this project was becoming more comfortable with the idea that there isn't always one correct way to approach a problem. When I was working through the functions and when we came together to compare our code, I came across various ways to do the same thing and we had made quite a few different decisions, but a lot of our implementations were actually doing the same thing. I think previously I would have been more concerned about figuring out which approach was the "right" one, or even jumped to doubt my code over anothers but this time I felt more comfortable looking at why each approach worked, deciding what made the most sense for our final code and most importantly to me confidently discuss what I wrote and why.

I also felt like I understood the algorithms much better this time around. Being able to reuse Smith-Waterman from Project 05 helped because I wasn't starting completely from scratch. The Neighbor-Joining calculations still took me a while to understand, especially when it came to updating the distance matrix and calculating branch lengths, but I felt like I could actually follow what the code was doing rather than just trying to get it to run.

Something I really appreciated about this project was that I had more time to play around with testing. I used the walkthrough Marcus posted throughout the project, and it was really helpful to have something I could keep referring back to as I worked through the different steps. I was able to test things as I went, compare my results, and go back to figure out why something wasn't working the way I expected. I think getting to test made a big difference in how well I understood my own code, rather than just getting it to work and moving on.

I still feel like I'm slower at coding than I'd like to be, even when I understand the concepts. But I think I'm getting better at working through problems without feeling like everything needs to be perfect from the beginning. Overall, I feel more confident discussing my code with other people, questioning my decisions, and explaining why I approached something a certain way. I think that's been one of the biggest improvements for me over these last few projects.

### Other Members: Mara
Going into this project, I had only ever heard of Neighbor Joining as a vague concept in prior discussions of general phylogeny, so this was my first real introduction to it. The lecture was very confusing to me, but once provided the walkthrough, the process became much clearer, 
and adjusting the Smith-Waterman algorithm to return normalized scores was a very simple task. Building the distance matrix using NumPy and indexing methods was straight forward as well, thanks to previous projects.
I was able to create a logic flow for the Neighbor-Joining algorithm by following the walkthrough, but the actual implementation and use of recursion ended up being more complicated. 
A point of confusion for me following the walkthrough was why they chose the final node joining at 3 nodes left instead of 2. My groupmates explained to me that stop condition of 3 nodes left is valid, but makes for a less detailed tree.  

My groupmates and I were able to talk through the formulas we used and the decisions we made in order to come up with a combined approach for our final submission. I had originally used the "distance = 1-similarity" formula, but we ended up using the natural log, as it is more biologically relevant. 
I do feel that I am getting better at reading other people's code and understanding why they chose the method they did. 


## Appendix (AI usage, citations, etc.)
### References:
Feng, D.-F., & Doolittle, R. F. (1987). Progressive sequence alignment as a prerequisite to correct phylogenetic trees. Journal of Molecular Evolution, 25(4), 351–360. https://doi.org/10.1007/BF02603120

Gascuel, O., & Steel, M. (2006). Neighbor-joining revealed. Molecular Biology and Evolution, 23(11), 1997–2000. https://doi.org/10.1093/molbev/msl072

Kosakovsky Pond, S. L., Weaver, S., Leigh Brown, A. J., & Wertheim, J. O. (2018). HIV-TRACE (TRAnsmission Cluster Engine): A tool for large scale molecular epidemiology of HIV-1 and other rapidly evolving pathogens. Molecular Biology and Evolution, 35(7), 1812–1819. https://doi.org/10.1093/molbev/msy016

Rhee, S.-Y., Gonzales, M. J., Kantor, R., Betts, B. J., Ravela, J., & Shafer, R. W. (2003). Human immunodeficiency virus reverse transcriptase and protease sequence database. Nucleic Acids Research, 31(1), 298–303. https://doi.org/10.1093/nar/gkg100

Saitou, N., & Nei, M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution, 4(4), 406–425. https://doi.org/10.1093/oxfordjournals.molbev.a040454

Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. Journal of Molecular Biology, 147(1), 195–197. https://doi.org/10.1016/0022-2836(81)90087-5

Studier, J. A., & Keppler, K. J. (1988). A note on the neighbor-joining algorithm of Saitou and Nei. Molecular Biology and Evolution, 5(6), 729–731. https://doi.org/10.1093/oxfordjournals.molbev.a040527


### AI Usage
Katelyn: 
I used Claude Opus 5.5 for this project. I used it to help me understand the background material, and to work through small examples by hand. It guided me through writing the pseudocode and Python for `build_distance_matrix` and `neighbor_joining` by asking questions and giving me hints if I needed them. It helped me debug syntax, and suggested the fix for extracting the short sequence IDs in `read_fasta`. It also helped me with Git and VS Code setup as well as improving my line comments. For the README, I also used it to help me write the Documentation and Usage section, as well as the Design Decisions, and the Formulas, which I reviewed and edited. It gave me feedback on my description and reflection, which I wrote myself.

Aamna:
ChatGPT (GPT-6), 2026.
I used ChatGPT throughout this project to help me understand the Neighbor-Joining algorithm and work through parts of my implementation. I especially needed help understanding how the corrected distance matrix was calculated, how branch lengths were determined, and how the matrix was updated after joining two nodes. I also used it to help troubleshoot parts of my code, to test whether my implementation was producing the expected results, and to
help me with syntax and optimize commenting.

I also used ChatGPT when comparing different approaches to the project, particularly recursion versus iteration, stopping at two versus three nodes, and different ways of normalizing similarity scores and converting them into distances. It helped me understand the reasoning behind these decisions so I could better explain my choices during our group discussions.

Mara:
ChatGPT (GPT-6), 2026.
I also used ChatGPT to brainstorm pseudocode and to subsequently troubleshoot bugs. 
Prompt: "I am coding a neighbor joining phylogeny based on a normalized smith-waterman algorithm and am writing pseudocode. Do not write any code for me."