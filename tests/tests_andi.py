from marbel.meta_tran_sim import checknegative


def test_checknegative():
    input_value = 1.0
    result = checknegative(input_value)
    assert result == input_value
