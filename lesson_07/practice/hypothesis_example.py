from hypothesis import given, settings, Verbosity
from hypothesis import strategies as st

@given(st.lists(st.integers()))
def test_sort_list_properties(lst):
    sorted_lst = list(sorted(lst))

    # Свойство 1: Длина отсортированного списка равна длине исходного списка
    assert len(sorted_lst) == len(lst)

    # Свойство 2: Отсортированный список содержит те же элементы, что и исходный
    assert sorted(sorted_lst) == sorted(lst)

    # Свойство 3: Отсортированный список отсортирован в неубывающем порядке
    assert all(sorted_lst[i] <= sorted_lst[i + 1] for i in range(len(sorted_lst) - 1))