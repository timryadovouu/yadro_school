class Person:
    def __init__(self, name, age, hobbies):
        self.name = name
        self.age = age
        self.hobbies = hobbies

    def greet(self):
        print(f"Hello, my name is {self.name} and I am {self.age} years old.")


def calculate_area(radius, pi=3.14159):
    return pi * (radius ** 2)


def example_function(a, b, c, d=None, e=None):
    if a == b:
        print("a is equal to b")
        print("That's interesting!")
    elif a > b:
        print("a is greater than b")
        if c:
            print("c is True")
            if d:
                print("d is also True")
    else:
        print("a is less than b")


def complex_calculation(x, y, z):
    result = x + y * z - (x / y) * z
    return result


# Multiline list and dictionary
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
data = {
    "name": "Alice",
    "age": 30,
    "city": "Wonderland",
    "profession": "Adventurer",
    "skills": ["jumping", "running", "hiding"],
    "strengths": {"intelligence": 10, "agility": 8, "strength": 7},
}


def unformatted_function_with_long_parameters(
    parameter_one,
    parameter_two,
    parameter_three,
    parameter_four,
    parameter_five,
    parameter_six,
    parameter_seven,
    parameter_eight,
):
    print(
        "This function has a lot of parameters, and should ideally be reformatted for readability."
    )


class UnformattedClassWithLongMethods:
    def __init__(self, data):
        self.data = data

    def long_method_with_lots_of_logic(self, val):
        if val > 100:
            print("That's a big value!")
        elif val < 0:
            print("Negative values are not allowed!")
            if isinstance(self.data, dict):
                for k, v in self.data.items():
                    print(f"Key {k}, Value {v}")
            else:
                print("Data is not a dictionary")
        return val


def very_long_list():
    return [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
    ]


def unaligned_code_block():
    for i in range(5):
        print(f"Number {i}")
    print("Loop done")
    print("More code")
    if True:
        print("Condition met")
    else:
        print("Condition not met")


# Large dictionary with various values
large_dict = {
    "item1": "A very long string that probably needs to be wrapped across multiple lines for better readability",
    "item2": "Another long string",
    "item3": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    "item4": {"key1": "value1", "key2": "value2", "key3": "value3"},
}


def unnecessary_inline_lambda():
    return (lambda x: x ** 2)(5)
