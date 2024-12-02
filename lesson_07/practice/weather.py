import requests

class WeatherService:
    def __init__(self, api_url):
        self.api_url = api_url

    def get_weather(self, city):
        """Запрашивает погоду для указанного города через API"""
        response = requests.get(f"{self.api_url}/weather?city={city}")
        if response.status_code == 200:
            return response.json()
        else:
            response.raise_for_status()