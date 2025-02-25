import tempfile
import os
import shutil
import unittest

class FileProcessor:
    def __init__(self, filename):
        self.filename = filename

    def read_file(self):
        """Чтение содержимого файла"""
        with open(self.filename, 'r') as file:
            return file.read()

    def count_words(self):
        """Подсчёт количества слов в файле"""
        content = self.read_file()
        return len(content.split())

    def count_lines(self):
        """Подсчёт количества строк в файле"""
        with open(self.filename, 'r') as file:
            return sum(1 for _ in file)

class TestFileProcessor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Создание временной директории перед выполнением всех тестов"""
        cls.temp_dir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        """Удаление временной директории после завершения всех тестов"""
        cls.temp_dir.cleanup()

    def setUp(self):
        """Создание временного файла перед каждым тестом"""
        self.temp_file = tempfile.NamedTemporaryFile(delete=False, dir=self.temp_dir.name, mode='w+')
        self.temp_file.write("Hello World\nThis is a test file.\nIt contains multiple lines.\n")
        self.temp_file.seek(0)  # Возвращаем курсор в начало файла для последующего чтения
        self.processor = FileProcessor(self.temp_file.name)

    def tearDown(self):
        """Удаление временного файла после каждого теста"""
        self.temp_file.close()
        os.remove(self.temp_file.name)

    def test_read_file(self):
        """Тестирование метода read_file"""
        content = self.processor.read_file()
        expected_content = "Hello World\nThis is a test file.\nIt contains multiple lines.\n"
        self.assertEqual(content, expected_content)

    def test_count_words(self):
        """Тестирование метода count_words"""
        word_count = self.processor.count_words()
        self.assertEqual(word_count, 11)  # 11 слов в тестовом файле

    def test_count_lines(self):
        """Тестирование метода count_lines"""
        line_count = self.processor.count_lines()
        self.assertEqual(line_count, 3)  # 3 строки в тестовом файле

if __name__ == "__main__":
    unittest.main()