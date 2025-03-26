import pandas as pd
from marbel.data_generations import calc_zero_ratio, maximize, minimize


def test_minimize():
    assert minimize(3, 4) == True
    assert minimize(4, 4) == False
    assert minimize(123, 1) == False

def test_maximize():
    assert maximize(3, 4) == False
    assert maximize(4, 4) == False
    assert maximize(123, 1) == True

def test_calc_zero_ratio():
    lst = [[1,2,3,0,5,1],
           [2,4,5,2,3,1],
           [1,0,5,0,1,0],
           [1,0,5,2,8,9],
           [1,5,6,0,1,7]]
    df = pd.DataFrame(lst)

    assert calc_zero_ratio(df) == 0.2

