import pandas as pd
from marbel.data_generations import calc_zero_ratio, maximize, minimize
import pytest


@pytest.mark.parametrize("a, b, expected", [
    (3, 4, True),
    (4, 4, False),
    (123, 1, False),
])
def test_minimize(a, b, expected):
    assert minimize(a, b) == expected


@pytest.mark.parametrize("a, b, expected", [
    (3, 4, False),
    (4, 4, False),
    (123, 1, True),
])
def test_maximize(a, b, expected):
    assert maximize(a, b) == expected


def test_calc_zero_ratio():
    lst = [[1, 2, 3, 0, 5, 1],
           [2, 4, 5, 2, 3, 1],
           [1, 0, 5, 0, 1, 0],
           [1, 0, 5, 2, 8, 9],
           [1, 5, 6, 0, 1, 7]]
    df = pd.DataFrame(lst)
    assert calc_zero_ratio(df) == 0.2
