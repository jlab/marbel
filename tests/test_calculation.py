import pandas as pd
from marbel.data_generations import calc_zero_ratio, maximize, minimize


def test_minimize():
    assert minimize(3, 4)
    assert not minimize(4, 4)
    assert not minimize(123, 1)


def test_maximize():
    assert not maximize(3, 4)
    assert not maximize(4, 4)
    assert maximize(123, 1)


def test_calc_zero_ratio():
    lst = [[1, 2, 3, 0, 5, 1],
           [2, 4, 5, 2, 3, 1],
           [1, 0, 5, 0, 1, 0],
           [1, 0, 5, 2, 8, 9],
           [1, 5, 6, 0, 1, 7]]
    df = pd.DataFrame(lst)
    assert calc_zero_ratio(df) == 0.2

    df = pd.DataFrame({'a': [0, 2, 3], 'b': [4, 0, 6], 'c': [7, 8, 9], 'd': [10, 11, 0]})
    ratio = calc_zero_ratio(df)
    assert ratio == 0.25
