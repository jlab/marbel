import pytest
import pandas as pd
import numpy as np
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, draw_orthogroups_by_rate
from marbel.presets import AVAILABLE_SPECIES, pg_overview

#Tests for draw_random_species

def test_correct_length():
    number_of_species = 5
    result = draw_random_species(number_of_species)
    assert len(result) == number_of_species


def test_unique_species():
    number_of_species = 5
    result = draw_random_species(number_of_species)
    assert len(result) == len(set(result))


def test_too_many_species():
    number_of_species = len(AVAILABLE_SPECIES) + 1
    with pytest.raises(ValueError):
        draw_random_species(number_of_species)


def test_too_few_species():
    number_of_species = 0
    with pytest.raises(ValueError):
        draw_random_species(number_of_species)

#Tests for create_ortholgous_group_rates


def test_correct_number_of_orthogroups():
    number_of_orthogroups = 10
    max_species = 5
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species)
    
    assert len(result) == number_of_orthogroups, f"Expected {number_of_orthogroups}, but got {len(result)}"

def test_max_species_per_group_respected():
    number_of_orthogroups = 20
    max_species = 7
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species)
    
    assert np.all(result <= max_species), f"Some orthogroups exceed the max species of {max_species}"

def test_non_negative_groups():
    number_of_orthogroups = 15
    max_species = 8
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species)
    
    assert np.all(result >= 0), "There are negative values in the orthogroup sizes"

def test_random_seed_reproducibility():
    number_of_orthogroups = 5
    max_species = 10
    seed = 42
    
    result_1 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=seed)
    result_2 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=seed)
    
    np.testing.assert_array_equal(result_1, result_2, err_msg="Results are not identical for the same seed")

def test_different_seed_variation():
    number_of_orthogroups = 5
    max_species = 10
    
    result_1 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=1)
    result_2 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=2)
    
    assert not np.array_equal(result_1, result_2), "Results should differ for different seeds"

def test_scaled_orthogroups_not_exceed_max():
    number_of_orthogroups = 10
    max_species = 12
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species)
    
    assert np.max(result) <= max_species, f"Max species per group exceeded: {np.max(result)}"

def test_empty_orthogroups():
    number_of_orthogroups = 0
    max_species = 5
    with pytest.raises(ValueError):
        create_ortholgous_group_rates(number_of_orthogroups, max_species)


#Tests for filter_by_seq_id_and_phylo_dist

def test_no_filters():
    result = filter_by_seq_id_and_phylo_dist()
    pd.testing.assert_frame_equal(result, pg_overview)


def test_max_phylo_distance_filter():
    max_phylo_distance = 0.9
    expected = pg_overview[pg_overview["tip2tip_distance"] <= max_phylo_distance]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance)
    pd.testing.assert_frame_equal(result, expected)


def test_min_identity_filter():
    min_identity = 90
    expected = pg_overview[pg_overview["medium_identity"] >= min_identity]
    result = filter_by_seq_id_and_phylo_dist(min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)


def test_both_filters():
    max_phylo_distance = 0.9
    min_identity = 90
    expected = pg_overview[(pg_overview["tip2tip_distance"] <= max_phylo_distance) & (pg_overview["medium_identity"] >= min_identity)]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance, min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)


def test_no_matching_records():
    max_phylo_distance = 0.1
    min_identity = 100
    expected = pg_overview[(pg_overview["tip2tip_distance"] <= max_phylo_distance) & (pg_overview["medium_identity"] >= min_identity)]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance, min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)

