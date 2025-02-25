class Matrix:
    def __init__(self, data):
        self.data = data
        self.rows = len(data)
        self.cols = len(data[0]) if self.rows > 0 else 0

    def __mul__(self, other):
        """Умножение двух матриц"""
        if self.cols != other.rows:
            raise ValueError("Матрицы нельзя перемножить: несоответствующие размеры.")
        
        result = [[sum(self.data[i][k] * other.data[k][j] for k in range(self.cols)) 
                   for j in range(other.cols)] for i in range(self.rows)]
        return Matrix(result)

    def __eq__(self, other):
        """Сравнение двух матриц на равенство"""
        return self.data == other.data

    @staticmethod
    def identity(size):
        """Создать единичную матрицу заданного размера"""
        return Matrix([[1 if i == j else 0 for j in range(size)] for i in range(size)])

    def power(self, n):
        """Возведение матрицы в степень n"""
        result = Matrix.identity(self.rows)
        base = self

        while n > 0:
            if n % 2 == 1:
                result = result * base
            base = base * base
            n //= 2

        return result

    def get(self, row, col):
        """Получить элемент матрицы по индексу"""
        return self.data[row][col]

def fibonacci(n):
    """Вычислить n-е число Фибоначчи с использованием матричного возведения в степень"""
    if n == 0:
        return 0
    if n == 1:
        return 1

    fib_matrix = Matrix([[1, 1], [1, 0]])
    result_matrix = fib_matrix.power(n - 1)
    return result_matrix.get(0, 0)