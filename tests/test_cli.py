import pandas as pd
from marbel.cli import limitthreads, compression_level_callback
import pytest
import typer


def test_limitthreads():
    limit = limitthreads(200)
    assert limit == 128


@pytest.mark.parametrize("valid", [
    0, 1 ,2, 3, 4, 5, 6, 7, 8, 9, 11
])
def test_compression_level_callback_valid(valid):
    callback = compression_level_callback(valid)
    assert callback == valid


@pytest.mark.parametrize("invalid", [
    -1, 10, 12, 13, 14
])
def test_compression_level_callback_invalid(invalid):
    with pytest.raises(typer.BadParameter):
        compression_level_callback(invalid)
