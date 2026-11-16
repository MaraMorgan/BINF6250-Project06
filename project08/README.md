# BINF6250
# Project 08: The Viterbi Algorithm

## Description
This project implements the Viterbi algorithm for hidden Markov models (HMMs). Given an HMM (initial, transition, and emission probabilities) and a sequence of observed symbols, the Viterbi algorithm finds the single most probable sequence of hidden states that could have produced the observations.

The algorithm uses dynamic programming: it fills a matrix with the best log score for each state at each position, stores a traceback pointer for each cell, and then follows the pointers back from the best final state to recover the path. All calculations are done in log space to prevent numerical underflow on long sequences.

The code is organized as a `HiddenMarkovModel` class so that later projects can add the forward, backward, and Baum-Welch algorithms as new methods on the same model.


## Documentation/Usage
### Setup
Requirements:
- Python 3
- NumPy

Clone the repository and switch to the project branch:

    git clone https://github.com/dsouzaka/BINF6250.git
    cd BINF6250
    git checkout project08_PR

Install NumPy if needed:

    pip install numpy

Open the project notebook (`project08MM.ipynb`) in Jupyter or VS Code and run the cells in order. The class definition must be run before the example cells.

### Running the project
Create a model by passing the states, the observation symbols, and the three probability tables as NumPy arrays. The order of `states` must match the row order of every array, and the order of `observation_symbols` must match the column order of the emission array.

    import numpy as np

    states = ["Island", "Genome"]            # row 0 = Island, row 1 = Genome
    observations = ["A", "C", "G", "T"]      # emission columns in this order

    init_probs = np.array([0.2, 0.8])        # P(start in Island), P(start in Genome)

    trans_probs = np.array([
        [0.7, 0.3],                          # from Island: to Island, to Genome
        [0.1, 0.9]                           # from Genome: to Island, to Genome
    ])

    emit_probs = np.array([
        [0.1, 0.4, 0.4, 0.1],                # Island: A, C, G, T
        [0.3, 0.2, 0.2, 0.3]                 # Genome: A, C, G, T
    ])

    hmm = HiddenMarkovModel(
        states=states,
        observation_symbols=observations,
        initial_probability=init_probs,
        transition_probability=trans_probs,
        emission_probability=emit_probs
    )

    hmm.observation_sequence = "GGCACTGAA"
    path = hmm.viterbi()
    print(path)
    print(hmm.best_final_score)

Output:

    ['Genome', 'Genome', 'Genome', 'Genome', 'Genome', 'Genome', 'Genome', 'Genome', 'Genome']
    -13.929108456051065

`viterbi()` returns the most probable path as a list of state names. The natural-log probability of that path is stored in `hmm.best_final_score`.

### Functions
- **`__init__(self, states, observation)_symbol,initial_probability,transition_probability, emission_probability)`**: Initialize objects used in HMM

- **`convert_to_log(self)`**: Helper function to convert probabilities to log to avoid underflow

- **`initialization_viterbi(self)`**:  Initialize stuff for probabilities and traceback matrix. Set up a probability matrix using the initial probabilities and the first observation. Initialize a traceback matrix for path reconstruction. No transition probs here

- **`calculate_viterbi_score(self)`**: Function to calculate the maximum probabilities of each using the Viterbi algorithm. Max is for probability matrix. Argmax is for traceback matrix

- **`termination_viterbi(self)`**: Function to determine the maximums

- **`traceback_viterbi(self)`**: Set up traceback path for Viterbi algorithm for most-probable paths.

- **`viterbi(self)`**: Function in which we call calculate_viterbi_score and traceback_viterbi to find the optimal path.





### Output (Testing)
# Example initial probabilities (probability of starting in each state)
states = ["Genome", "Island"]
observations = ["A", "C", "G", "T"]
init_probs = [0.2, 0.8]
# "Island": 0.2, "Genome": 0.8


# Example transition probabilities (probability of moving from one state to another)
trans_probs = np.array([
    [0.7, 0.3],
    [0.1, 0.9]])

# "Island": {"I": 0.7, "G": 0.3}, "Genome": {"I": 0.1, "G": 0.9}
# Example emission probabilities (probability of observing a symbol in a given state)
emit_probs = np.array([
    [0.1, 0.4, 0.4, 0.1],
    [0.3, 0.2, 0.2, 0.3]
])

# "Island": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
# "Genome": {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}


# Example observation sequence
obs_seq = "GGCACTGAA"


hmm = HiddenMarkovModel(
    states=states,
    observation_symbols=observations,
    initial_probability=init_probs,
    transition_probability=trans_probs,
    emission_probability=emit_probs
)
hmm.observation_sequence = obs_seq
result = hmm.viterbi()
print(result)


['Genome', 'Genome', 'Genome', 'Island', 'Island', 'Island', 'Island', 'Island', 'Island']



## Pseudocode
### convert_to_log
'''
Convert each probability array to log probability array; store as new attribute

'''

### initialization_viterbi
'''
Get the number of hidden states and the length of the observation sequence
    # Rows = number of hidden states
    # Columns = length of observation sequence
Initialize dimensions of both probability and traceback matrices
    # viterbi_matrix = np.zeros((rows, columns))
    # viterbi_traceback_matrix = np.zeros((rows, columns), dtype=int)
Initialize the first column of the traceback with a character because there aren't previous states to the first state
Getting the first observation and its emission probability
    # Establish the first observation at index 0 for observation sequence
    # Create a symbol index to get the probability for the symbol at observation sequence index 0
    # Get the emission probabilities for all rows in column 1 at symbol index for that position in sequence
Calculation for first column of probabilities
    # initial probability from symbol, index + the first emission probability (remember we are in log, so its addition)
Store in first column of the probability matrix
'''

### calculate_viterbi_score
'''
Iterate through each observation
iterate through each possible hidden state at that position
calculate probabilities for each hidden state from previous state
find the max
Add that to emission probability and store 
store the previous state that produced it for the traceback
'''

### termination_viterbi
'''
Determine the last index of the observation sequence
Determine the best final state
Determine the best final probability
'''

### traceback_viterbi
''' 
Initialize with the best final state
Set current state equal to that
iterate backwards through the observations
look at previous state based on current state and position
update current state to previous state
add previous state into path
iterate until start
reverse the final order
'''


### Formulas
### The Viterbi Algorithm
Dynamic programming

$\pi_i$ = initial probabilities

* the chance that the first position is in state $i$

$a_{ij}$ = transition probabilities

* the chance that the next position is in state $j$ when this position is in state $i$

$b_j(k)$ = emission probabilities

* the chance that a position in state $j$ shows symbol $k$

Each row of the transition table and each row of the emission table adds up to 1, because each row lists every possible outcome from one starting point.

$\lambda = (A, B, \pi)$

* the three sets together

Assumptions (Markov property):

* The next state depends only on the current state, not on any earlier position
* Each observation depends only on the state at its own position

### Formal Description by Durbin et al.

Durbin writes emissions as $e_\ell(x_i)$ instead of $b_j(k)$; they mean the same thing. Durbin also adds a silent begin/end state, numbered 0.

#### Initialization ($i = 0$):

$$v_0(0) = 1, \quad v_k(0) = 0 \text{ for } k > 0$$

#### Recursion ($i = 1 \ldots L$):

$$v_\ell(i) = e_\ell(x_i) \cdot \max_k \big( v_k(i-1) \, a_{k\ell} \big)$$

$$\text{ptr}_i(\ell) = \arg\max_k \big( v_k(i-1) \, a_{k\ell} \big)$$

#### Termination:

$$P(x, \pi^*) = \max_k \big( v_k(L) \, a_{k0} \big)$$

$$\pi_L^* = \arg\max_k \big( v_k(L) \, a_{k0} \big)$$

#### Traceback ($i = L \ldots 1$):

$$\pi_{i-1}^* = \text{ptr}_i(\pi_i^*)$$

#### Glossary

* $\max$: the highest value
* $\arg\max$: which state gave the highest value (stored as the traceback pointer)
* $e_\ell(x_i)$: emission probability of symbol $x_i$ from state $\ell$
* $a_{k\ell}$: transition probability from state $k$ to state $\ell$
* $k$: the previous state (the one maximized over)
* $\ell$: the current state (the cell being filled)
* $L$: the length of the observation sequence
* $i$: the position in the sequence
* $x$: the observation sequence; $x_i$ is the symbol at position $i$
* $\pi^*$: the most probable path; $\pi_i^*$ is its state at position $i$ (not the same as the initial probabilities $\pi_i$)

### In our implementation

* Our model has no end state, so termination drops $a_{k0}$:

$$P(x, \pi^*) = \max_k v_k(L), \qquad \pi_L^* = \arg\max_k v_k(L)$$

* We are given initial probabilities directly, so column 1 is $v_k(1) = \pi_k \, e_k(x_1)$ and the recursion starts at $i = 2$.
* To avoid underflow, we work in log space (max stays max; products become sums):

$$\log v_\ell(i) = \log e_\ell(x_i) + \max_k \big( \log v_k(i-1) + \log a_{k\ell} \big)$$

### Design Decisions

## Sucesses

## Struggles


## Reflections
### Group Leader: Katelyn
Starting this project, I had somewhat of an understanding of the Viterbi algorithm. Mara created an outline including the formulas and a basic structure for us to follow, and having that to reference and to follow was extremely helpful since I definitely struggled with not having any skeleton code this project. We all worked on creating our own implementations from this structure and then had a meeting to discuss any thoughts we had, like how to address different data types and what structures would work best. Our group worked well together to help fill gaps in eachother's understanding.
I used Claude to help walk me through examples of the Viterbi algorithm, which I felt really helped me understand what we were doing with this project. I also had it help me troubleshoot areas where I was stuck with the pseudocode, and I definitely found this project more difficult than past ones. Once we had the pseudocode, the Python was not too bad. 
Since this project is the start of many others, we had to make sure that we prepared for that. This was definitely different than other projects because we didn't previously need to change things so they could be expanded upon, but this process seems to better represent what this process would look like in the real world. I am interested to see how this project will be expanded upon in the future modules and how we can add more algorithms to our HMM class.

### Other Members: 

### Other Members: 



## Appendix (AI usage, citations, ect.)
### References
Feng, D.-F., & Doolittle, R. F. (1987). Progressive sequence alignment as a prerequisite to correct phylogenetic trees. Journal of Molecular Evolution, 25(4), 351–360. https://doi.org/10.1007/BF02603120

Gascuel, O., & Steel, M. (2006). Neighbor-joining revealed. Molecular Biology and Evolution, 23(11), 1997–2000. https://doi.org/10.1093/molbev/msl072

Kosakovsky Pond, S. L., Weaver, S., Leigh Brown, A. J., & Wertheim, J. O. (2018). HIV-TRACE (TRAnsmission Cluster Engine): A tool for large scale molecular epidemiology of HIV-1 and other rapidly evolving pathogens. Molecular Biology and Evolution, 35(7), 1812–1819. https://doi.org/10.1093/molbev/msy016

Rhee, S.-Y., Gonzales, M. J., Kantor, R., Betts, B. J., Ravela, J., & Shafer, R. W. (2003). Human immunodeficiency virus reverse transcriptase and protease sequence database. Nucleic Acids Research, 31(1), 298–303. https://doi.org/10.1093/nar/gkg100

Saitou, N., & Nei, M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution, 4(4), 406–425. https://doi.org/10.1093/oxfordjournals.molbev.a040454

Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. Journal of Molecular Biology, 147(1), 195–197. https://doi.org/10.1016/0022-2836(81)90087-5

Studier, J. A., & Keppler, K. J. (1988). A note on the neighbor-joining algorithm of Saitou and Nei. Molecular Biology and Evolution, 5(6), 729–731. https://doi.org/10.1093/oxfordjournals.molbev.a040527 


### AI Usage
Katelyn: I used Claude (Anthropic) throughout this project, mainly to understand the algorithm and check work. The group's main implementation was written by Mara. I used Claude to walk me through the Viterbi algorithm step by step, explain related concepts such as conditional independence and object-oriented programming, review and correct pseudocode and the formatting of the formulas for the notebook, guide me through writing a dictionary-based implementation by pointing out errors and giving me hints while I wrote most of the methods myself, test implementations against hand-checked examples to identify bugs and edge cases we missed, to draft the Setup and Running the Project sections of the README, and revise my Description section. 
