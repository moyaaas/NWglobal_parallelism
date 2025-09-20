import argparse
import numpy as np
from multiprocessing import Pool, cpu_count

DIAG = 0
UP = 1
LEFT = 2


def initialize_matrix(seq1: str, seq2: str, gap_penalty: int) -> np.ndarray: # la classe che rappresenta gli array in numpy è np.ndarray, np.array è il costruttore
    ''' 
    This function initializes the score and directions matrices for the global alignment.
    
    Parameters
    ----------
    seq1: str
        The first input sequence.
    seq2: str
        The second input sequence.
    gap_penalty: int
        The penalty to apply for introducing a gap.
    
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        A tuple containing:
            - the scoring matrix (filled with only gap penalties for the first row and column)
            - the directions matrix (with one or more directions taken to reach each cell)
    '''
    nrow, ncol = len(seq1), len(seq2)
    matrix = np.zeros((nrow + 1, ncol + 1), dtype = int)
    directions_matrix = np.empty((nrow + 1, ncol + 1), dtype = object)
    for m in range(nrow +1):
        for n in range(ncol + 1):
            directions_matrix[m,n] = []
    for n in range(1, ncol + 1): 
        matrix[0][n] = matrix[0][n - 1] + gap_penalty
        directions_matrix[0][n] = [LEFT]
    for m in range(1, nrow + 1):
        matrix[m][0] = matrix[m -1][0] + gap_penalty
        directions_matrix[m][0] = [UP]
    return matrix, directions_matrix



def calculate_cell_score(up, diag, left, char1, char2, match_score:int, mismatch_score:int, gap_penalty:int):
    '''
    Calculate the optimal score and directions for a specific cell.

    Parameters
    ----------
    matrix: np.ndarray
        The score matrix.
    up : int 
        score from the cell above
    diag : int
        score from the diagonal cell
    left: int
        score from the cell to the left
    char1: str
        character from the seq1, current row.
    char2: str
        character from the seq2, current column
    match_score: int
        The score to assign for a match.
    mismatch_score: int
        the score to assign for a mismatch.
    gap_penalty: int
        The penalty to apply for introducing a gap.
    
    Returns
    -------
    tuple((int, list[int]))
       - the score for the current cell 
       - a list of possible directions to follow to reach the current cell
    '''
    global DIAG, UP, LEFT
    up_score = up + gap_penalty
    diag_score = diag + (match_score if char1 == char2 else mismatch_score) # ternary operator
    left_score = left + gap_penalty
    max_score = max(up_score, diag_score, left_score) 
    candidates = [(UP, up_score), (DIAG, diag_score), (LEFT, left_score)]
    directions = [dir for dir, score in candidates if score == max_score]
    
    return max_score, directions

def compute_antidiagonal(nrow : int, ncol : int): 
    '''
    Compute the antidiagonals of a matrix with dimensions nrow, ncol.
    
    Parameters
    ----------
    nrow : int 
        length of seq1 that is equals to the number of rows.
    ncol : int 
        length of seq2 that is equal to the number of columns.  

    Returns
    -------
    list[list[tuple(int, int)]]
        list of antidiagonals of the matrix, each list is a list of positions.
   '''
    antidiagonal = [[] for i in range(nrow + ncol +1)]
    for m in range(nrow + 1):
        for n in range(ncol +1): 
            antidiagonal[m+n].append((m,n)) 
    return antidiagonal
            

def fill_matrix(matrix: np.ndarray, seq1: str, seq2:str, directions_matrix: np.ndarray, antidiagonal:list, match_score:int, mismatch_score:int, gap_penalty:int):
    '''
    For each antidiagonal, for each cell it's computed the score and stored the direction in a separate matrix.

    Parameters
    ----------
        Parameters
    ----------
    matrix: np.ndarray
        The score matrix.
    seq1: str
        The first input sequence.
    seq2: str
        The second input sequence.
    directions_matrix:
        The matrix storing directions taken to reach each cell.
    antidiagonal: list
        List of matrix coordinates grouped by antidiagonals.
    match_score: int
        The score to assign for a match.
    mismatch_score: int
        the score to assign for a mismatch.
    gap_penalty: int
        The penalty to apply for introducing a gap.
    
    Returns
    -------
    None
    '''
    with Pool() as pool:
        for antid in antidiagonal:
            args_list = [] # to store arguments for each cell
            coords = [] # to store coordinates of the cell
            for m,n in antid:
                if m == 0 or n == 0: # first row and column have been already initialized, so they can be skipped
                    continue
                # extract arguments needed to calculate the score
                up = matrix[m-1, n]
                diag = matrix[m-1,n-1]
                left = matrix[m, n-1]
                char1 = seq1[m-1]
                char2= seq2[n-1]
                # store arguments
                args_list.append((up, diag, left, char1, char2, match_score, mismatch_score, gap_penalty))
                coords.append((m,n))
            if not args_list: 
                continue
            # compute score and direction for all cells in the antidiagonal in parallel
            results = pool.starmap(calculate_cell_score, args_list)
            #updating the matrices
            for i, (score, dirs) in enumerate(results):
                m,n = coords[i]
                matrix[m,n] = score
                directions_matrix[m,n] = dirs


def traceback (directions_matrix: np.ndarray, seq1:str, seq2:str):
    global DIAG, UP, LEFT
    '''
    Performs the traceback step to reconstruct all optimal alignments.
    Start from the last cell (bottom-right) and follows the directions back to the first cell (upper-left).

    Parameters:
    ----------
     directions_matrix:
        The matrix storing directions taken to reach each cell.
    seq1: str
        The first input sequence.
    seq2: str
        The second input sequence.
   
    Returns:
    --------
    set[tuple[str,str]]
        A set of tuple containing two aligned sequences.
    '''
    nrow, ncol = len(seq1), len(seq2)
    stack = [(nrow, ncol, "", "")]
    aligned_seq = set()
    while stack:
        m, n, a1, a2 = stack.pop()
        if m == 0 and n == 0:
            aligned_seq.add((a1,a2))
            continue
        for dir in directions_matrix[m][n]:
            if dir == DIAG:
                new_a1 = seq1[m-1] +a1
                new_a2 = seq2[n-1] +a2
                stack.append((m-1,n-1, new_a1, new_a2))
            elif dir == UP:
                new_a1 = seq1[m-1] +a1
                new_a2 = '-' + a2
                stack.append((m-1,n, new_a1, new_a2))
            elif dir == LEFT:
                new_a1 = '-' + a1
                new_a2 = seq2[n-1] +a2
                stack.append((m,n-1, new_a1, new_a2))
    return aligned_seq




def main():
    '''
    Executes the Needleman - Wunsch global alignment.
        - Parses input sequences and scoring parameters from the command line
        - Initializes the scoring and direction matrices
        - Fills the matrices using antidiagonal 
        - Performs traceback to generate all optimal alignments
        - Prints and saves the alignments and score to a file
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument("seq1", type = str, help = 'first input sequence')
    parser.add_argument("seq2", type = str, help = 'second input sequence')
    parser.add_argument("match_score", type = int, help = 'a positive integer')
    parser.add_argument("mismatch_score", type = int, help = 'a negative integer')
    parser.add_argument("gap_penalty", type = int, help = 'a negative integer')

    args = parser.parse_args()

# Input sequences
    seq1 = args.seq1
    seq2 = args.seq2

# Dimensions of the alignment matrix
# m rows
# n columns
    nrow, ncol = len(seq1), len(seq2)

# Alignment Parameters
    match_score = args.match_score
    mismatch_score = args.mismatch_score 
    gap_penalty = args.gap_penalty

# Initializing matrices
    matrix, directions_matrix = initialize_matrix(seq1, seq2, gap_penalty)

# Computing antidiagonals   
    antidiagonal = compute_antidiagonal(len(seq1), len(seq2))

# Filling the scoring matrix using antidiagonals
    fill_matrix(matrix, seq1, seq2, directions_matrix, antidiagonal, match_score, mismatch_score, gap_penalty)

# Performing traceback
    aligned_seq = traceback(directions_matrix, seq1, seq2)
    final_score = matrix[nrow][ncol]

# Writing output file
    with open("globalNW_results.txt", "w") as f:
        f.write(f'Number of optimal alignments found: {len(aligned_seq)}\n\n')
        f.write(f"Alignment score: {final_score:}\n\n")
        for a1, a2 in aligned_seq:
            f.write(f"{a1}\n")
            f.write(f"{a2}\n\n")
            


if __name__ == "__main__":
    main()

