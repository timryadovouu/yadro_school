import os
import pytest
from tempfile import NamedTemporaryFile

def test_addition():
    assert 1 + 1 == 2

def test_subtraction():
    assert 5 - 3 == 2
    
    
@pytest.fixture
def sample_list():
    # Создаём тестовые данные
    return [1, 2, 3, 4, 5]

def test_sample_list_length(sample_list):
    assert len(sample_list) == 5

def test_sample_list_sum(sample_list):
    assert sum(sample_list) == 15
    

@pytest.mark.parametrize("a, b, expected", [
    (1, 1, 2),
    (2, 3, 5),
    (3, 5, 8),
    (10, 20, 30)
])
def test_addition(a, b, expected):
    assert a + b == expected
    

@pytest.mark.slow
def test_slow_function():
    import time
    time.sleep(2)
    assert True

@pytest.mark.fast
def test_fast_function():
    assert True
    
    
@pytest.fixture()
def db():
    pass

@pytest.fixture()
def tempfile():
    with NamedTemporaryFile(dir=os.getcwd()) as f:
        yield f

def test_with_tempfile_and_db(tempfile, db):
    pass

def test_with_db(db):
    pass

def test_with_tempfile(tempfile):
    pass


from weather import WeatherService

def mock_get_weather(self, city):
    return {
        "city": city,
        "temperature": 20,
        "description": "Sunny"
    }

def test_get_weather_with_monkeypatch(monkeypatch):
    # Заменяем метод get_weather на фиктивную функцию
    monkeypatch.setattr(WeatherService, "get_weather", mock_get_weather)

    service = WeatherService(api_url="http://fakeapi.com")
    result = service.get_weather("Paris")

    assert result == {
        "city": "Paris",
        "temperature": 20,
        "description": "Sunny"
    }