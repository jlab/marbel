import pandas as pd
from marbel.cli import limitthreads
import pytest


def test_limitthreads():
    limit = limitthreads(200)
    assert limit == 128