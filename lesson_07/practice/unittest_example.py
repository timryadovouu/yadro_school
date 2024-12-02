import unittest
from unittest.mock import MagicMock


class MathTestCase(unittest.TestCase):
    def test_addition(self):
        self.assertEqual(1 + 1, 2)

    def test_subtraction(self):
        self.assertEqual(5 - 3, 2)
        
        
class DatabaseTestCase(unittest.TestCase):
    def setUp(self):
        # Этот метод выполняется перед каждым тестом
        self.db = "Подключение к базе данных"

    def tearDown(self):
        # Этот метод выполняется после каждого теста
        self.db = None

    def test_connection(self):
        self.assertIsNotNone(self.db)
        
        
    def test_connection(self):
        self.assertIsNotNone(self.db)

class AssertionTestCase(unittest.TestCase):
    def test_assert_methods(self):
        self.assertEqual(3 + 2, 5)
        self.assertTrue(5 > 3)
        self.assertIn(3, [1, 2, 3])

    def test_assert_raises(self):
        with self.assertRaises(ZeroDivisionError):
            1 / 0

class SkipTestCase(unittest.TestCase):
    @unittest.skip("Пропускаем этот тест")
    def test_skipped(self):
        self.fail("Этот тест пропущен")

    @unittest.skipIf(1 > 0, "Пропускаем, так как 1 > 0")
    def test_skipped_conditionally(self):
        self.fail("Этот тест также пропущен")
        
        
class ClassSetupTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ### Один раз перед и после выполнения всех тестов в классе
        cls.resource = "Общая настройка ресурса"

    @classmethod
    def tearDownClass(cls):
        cls.resource = None

    def test_example(self):
        self.assertEqual(self.resource, "Общая настройка ресурса")


class APITestCase(unittest.TestCase):
    def test_api_call(self):
        mock_api = MagicMock()
        mock_api.get_data.return_value = {"data": "example"}
        
        response = mock_api.get_data()
        self.assertEqual(response, {"data": "example"})
        mock_api.get_data.assert_called_once()
        
if __name__ == '__main__':
    unittest.main()