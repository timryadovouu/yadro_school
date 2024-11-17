import unittest
from fib import Matrix, fibonacci

class TestFibonacciIntegration(unittest.TestCase):
    def test_fibonacci_with_matrix_power(self):
        # Проверка первых чисел Фибоначчи
        self.assertEqual(fibonacci(0), 0)
        self.assertEqual(fibonacci(1), 1)
        self.assertEqual(fibonacci(2), 1)
        self.assertEqual(fibonacci(3), 2)
        self.assertEqual(fibonacci(4), 3)
        self.assertEqual(fibonacci(5), 5)
        self.assertEqual(fibonacci(6), 8)
        self.assertEqual(fibonacci(10), 55)

    def test_large_fibonacci(self):
        # Проверка числа Фибоначчи для больших значений n
        self.assertEqual(fibonacci(20), 6765)
        self.assertEqual(fibonacci(30), 832040)
        self.assertEqual(fibonacci(40), 102334155)
        self.assertEqual(fibonacci(50), 12586269025)

    def test_matrix_power_large(self):
        # Проверка возведения матрицы в большую степень
        matrix = Matrix([[1, 1], [1, 0]])
        result = matrix.power(10)
        expected = Matrix([[89, 55], [55, 34]])  # Matrix for F(11) and F(10)
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()