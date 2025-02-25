import unittest
import requests
from weather import WeatherService
from unittest.mock import patch, MagicMock


class TestWeatherService(unittest.TestCase):

    @patch('weather.requests.get')
    def test_get_weather_success(self, mock_get):
        # Настройка мока для успешного ответа от API
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "city": "London",
            "temperature": 15,
            "description": "Partly cloudy"
        }
        mock_get.return_value = mock_response

        # Создаем экземпляр WeatherService и вызываем метод get_weather
        service = WeatherService(api_url="http://fakeapi.com")
        result = service.get_weather("London")

        # Проверяем, что результат соответствует ожидаемому
        self.assertEqual(result, {
            "city": "London",
            "temperature": 15,
            "description": "Partly cloudy"
        })
        
        # Проверяем, что requests.get был вызван с нужными аргументами
        mock_get.assert_called_once_with("http://fakeapi.com/weather?city=London")

    @patch('weather.requests.get')
    def test_get_weather_failure(self, mock_get):
        # Настройка мока для ответа с ошибкой от API
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError("404 Not Found")
        mock_get.return_value = mock_response

        # Создаем экземпляр WeatherService
        service = WeatherService(api_url="http://fakeapi.com")

        # Проверяем, что вызывается исключение HTTPError при ошибке 404
        with self.assertRaises(requests.exceptions.HTTPError):
            service.get_weather("UnknownCity")

        # Проверяем, что requests.get был вызван с нужными аргументами
        mock_get.assert_called_once_with("http://fakeapi.com/weather?city=UnknownCity")

if __name__ == '__main__':
    unittest.main()