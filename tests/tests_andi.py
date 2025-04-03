import pytest
import pandas as pd
import numpy as np
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, draw_orthogroups_by_rate, draw_dge_factors
from marbel.presets import AVAILABLE_SPECIES, pg_overview, DGE_LOG_2_CUTOFF_VALUE
from marbel.meta_tran_sim import checknegative


def test_checknegative():
    input_value = 1.0
    result = checknegative(input_value)
    assert result == input_value