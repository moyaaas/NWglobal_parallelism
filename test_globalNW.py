from globalNW import initialize_matrix, calculate_cell_score, fill_matrix, compute_antidiagonal, traceback, DIAG, UP, LEFT
import numpy as np
import unittest
from unittest.mock import patch, MagicMock

class TestInitialization(unittest.TestCase):
    def test_initialize_matrix(self):
        seq1 = "ACGC"
        seq2 = "GACTAC"
        gap_penalty = -1
        matrix, directions = initialize_matrix(seq1,seq2, gap_penalty)
        self.assertEqual(matrix.shape, (len(seq1)+1,len(seq2)+1))
        self.assertEqual(matrix[0][0], 0)
        self.assertEqual(matrix[0][1], -1)
        self.assertEqual(matrix[1][0], -1)
        self.assertEqual(matrix[4][0], -4)
        self.assertEqual(matrix[0][6], -6)
        self.assertEqual(directions[0][1], [LEFT])
        self.assertEqual(directions[1][0], [UP])

class TestCellScore(unittest.TestCase):
    def setUp(self):
        self.params = {
            "match_score" : 1,
            "mismatch_score" : -1, 
            "gap_penalty": -1
        }
    def test_match(self):
        score, dirs = calculate_cell_score(
            up = 0, diag = 0, left = 0,
            char1 = "A", char2 = "A",
            ** self.params
        )
        self.assertEqual(score, 1)
        self.assertEqual(dirs, [DIAG])
    
    def test_mismatch(self):
        score, dirs = calculate_cell_score(
            up = 0, diag = 0, left = 0,
            char1 = "A", char2 = "T",
            ** self.params
        )
        self.assertEqual(score, -1)
        self.assertEqual(sorted(dirs), [DIAG,UP,LEFT])

    def test_breaking(self):
        score, dirs = calculate_cell_score(
            up = 3, diag = 3, left = 3,
            char1 = "A", char2 = "T",
            ** self.params 
        )
        # 3 + (-1) = 2? diag=3 + mismatch(-1)=2, up=3 + gap(-1)=2, left=3 + gap(-1)=2
        self.assertEqual(score,2)
        self.assertEqual(sorted(dirs), [DIAG, UP, LEFT])

class TestAntidiagonal(unittest.TestCase):
    def test_mini_matrix(self):
        antidiagonals = compute_antidiagonal(2,2)
        expected = [[(0,0)], 
                    [(0,1),(1,0)], 
                    [(0,2),(1,1),(2,0)],
                    [(1,2),(2,1)],
                    [(2,2)]]
        self.assertEqual(antidiagonals, expected)

class TestFillingMatrix(unittest.TestCase):
    def test_alignment(self):
        matrix = np.zeros((2,2), dtype = int)
        dirs = np.empty((2,2), dtype = object)
        for m in range (2):
            for n in range(2):
                dirs[m,n] = []
        matrix[0,1] = -1
        matrix[1,0] = -1
        antidiagonals = compute_antidiagonal(1,1)
        fill_matrix(matrix, "A","A", dirs, antidiagonals, 1, -1, -1)
        self.assertEqual(matrix[1,1],1)
        self.assertEqual(dirs[1,1], [DIAG])

class TestTraceback(unittest.TestCase):
    def test_one_alignment(self):
        dirs = np.empty((3,3), dtype = object)
        for m in range(3):
            for n in range(3):
                dirs[m,n] = []
        dirs[2,2] = [DIAG]
        dirs[1,1] = [DIAG]
        alignments = traceback(dirs, "AT", "AT")
        self.assertEqual(alignments, {("AT","AT")})



class TestParallelism(unittest.TestCase):
    @patch("globalNW.Pool")
    def test_parallelism(self, mock_pool_class):
        mock_pool = MagicMock()
        mock_pool.__enter__.return_value.starmap.return_value = [(1, [DIAG])]
        mock_pool_class.return_value = mock_pool

        matrix = np.zeros((2,2), dtype = int)
        dirs = np.empty((2,2), dtype = object)
        for m in range(2):
            for n in range(2):
                dirs[m,n] = []     
        matrix[0,0] = 0 #diag
        matrix[0,1] = -1 #up
        matrix[1,0] = -1 #left   
        antidiagonals = [[(1,1)]]
        fill_matrix(matrix, "A", "A", dirs, antidiagonals, 1, -1, -1)

        #to verify that parallel process is called, Pool() and starmap
        mock_pool_class.assert_called_once()
        mock_pool.__enter__.return_value.starmap.assert_called_once()
        
        # arguments passed to starmap --> call_args --> a tuple (args, kwargs) 
        call_args = mock_pool.__enter__.return_value.starmap.call_args # a tuple
        args, kwargs = call_args
        # args is a tuple: the function 'calculate_cell_score' and a list of tuples for each cell in the antidiagonal
        self.assertEqual(len(args),2)
        args_list = args[1] # list of tuples
        self.assertEqual(len(args_list), 1) # we have only one cell in the antidiagonal
        # check expected arguments
        # up = matrix[0,1] = -1
        # diag = matrix[0,0] = 0
        # left = matrix[1,0] = -1
        expected_args = (-1,0,-1, "A","A", 1,-1,-1)
        self.assertEqual(args_list[0], expected_args)
        
        # verify results
        self.assertEqual(matrix[1,1], 1)
        self.assertEqual(dirs[1,1], [DIAG])
        
        
if __name__ == "__main__":
    print("Test...")
    unittest.main()


