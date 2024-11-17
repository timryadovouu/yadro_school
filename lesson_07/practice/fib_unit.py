import unittest
from fib import Matrix, fibonacci

class TestMatrixUnit(unittest.TestCase):
    def test_matrix_multiplication(self):
        matrix_a = Matrix([[1, 1], [1, 0]])
        matrix_b = Matrix([[1, 1], [1, 0]])
        result = matrix_a * matrix_b
        expected = Matrix([[2, 1], [1, 1]])
        self.assertEqual(result, expected)

    def test_identity_matrix(self):
        identity = Matrix.identity(2)
        expected = Matrix([[1, 0], [0, 1]])
        self.assertEqual(identity, expected)

    def test_matrix_power(self):
        matrix = Matrix([[1, 1], [1, 0]])
        result = matrix.power(2)
        expected = Matrix([[2, 1], [1, 1]])
        self.assertEqual(result, expected)

    def test_get_element(self):
        matrix = Matrix([[1, 2], [3, 4]])
        self.assertEqual(matrix.get(0, 1), 2)
        self.assertEqual(matrix.get(1, 0), 3)

if __name__ == "__main__":
    unittest.main()