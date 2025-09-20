# Global Alignment with Needleman-Wunsch 
This project implements the **Needleman-Wunsch algorithm** for global sequence alignment with parallel computation of antidiagonals, using Python's **multiprocessing** module. It calculetes the optimal algnments between two sequences, using a score and a direction matrix and traceback.

## Features
- Global Aignment
- Customizable scoring system
- Parallel computation of antidiagonal scores
- Traceback to find the optimal alignments
- Output all optimal alignments
- Results saved to a file

## File Structure
globalNW.py - Main script 
test_globalNW.py - Test script
README.md - Documentation, this file
globalNW_results.txt - Output file with alignments results

## Requirements
- Python 3.x
- Numpy
- Multiprocessing
- unittest (for testing)
- argparse (for command line arguments)

## Running the main script
To run the script, execute the following command in your terminal:
'''bash
python globalNW.py <sequence1> <sequence2> <match_score> <mismatch_score> <gap_penalty>
'''
example:
'''bash
python globalNW.py ACGC GACTAC 0 -1 -1
'''
This will align the sequences ACGC and GACTAC with a match score of 0, a mismatch score of -1, and a gap penalty of -1.

You can also run the script with the '-h' or '--help' option to see the help message:
'''bash
python globalNW.py -h
'''
This will display the following help message:
usage: globalNW.py [-h] seq1 seq2 match_score mismatch_score gap_penalty

positional arguments:
  seq1            first input sequence
  seq2            second input sequence
  match_score     a positive integer
  mismatch_score  a negative integer
  gap_penalty     a negative integer

options:
  -h, --help      show this help message and exit

-h --help: Show help message and exit

## Output
The output will be saved in the file 'globalNW_results.txt'.
The output file will contain:
Number of optimal alignments found: 2

Alignment score: -3

-AC-GC
GACTAC

-ACG-C
GACTAC

## Testing
To run the tests, execute the following command in your terminal:
'''bash
python -m unittest test_globalNW.py
'''
This will run the unit tests in the 'test_globalNW.py' file.

## Functions
- `initialize_matrix(seq1, seq2, gap_penalty)`: Initializes the scoring and direction matrices.
- `calculate_cell_score(up, diag, left, char1, char2, match_score, mismatch_score, gap_penalty)`: Calculates the score for a cell in the scoring matrix based on the scores of neighboring cells and the sequence characters.
- `compute_antidiagonal(nrow, ncol)`: Computes the antidiagonal indices for the scoring matrix.
- `fill_matrix(matrix, seq1, seq2, directions_matrix, antidiagonal, match_score, mismatch_score, gap_penalty)`: Fills the scoring matrix using the Needleman-Wunsch algorithm and multiprocessing. The function fills the scoring matrix and the direction matrix by antidiagonals, allowing for parallel computation. Along an antidiagonal, all cells are indipendent, and can be computed in parallel.
It gathers all the input arguments for each cell, due to the fact that the previous antidiagonal has been already computed, and then computes the scores in parallel.
- `traceback(directions_matrix, seq1, seq2)`: Performs traceback to find the optimal alignments.

## Author
Martina Montonati
Master's Degree in Bioinformatics for Computational Genomics at the Università degli Studi di Milan and Politecnico di Milano
Scientific Programming course, 2025