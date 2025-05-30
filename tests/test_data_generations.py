import pytest
import pandas as pd
import numpy as np
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, draw_dge_factors
from marbel.data_generations import calc_zero_ratio, get_all_zero_genes, add_counts_to_large_orthogroups
from marbel import data_generations as dg
from marbel.presets import MAX_SPECIES, ErrorModel
from marbel.preload import get_pg_overview

# Tests for draw_random_species


def test_correct_length():
    number_of_species = 5
    result = draw_random_species(number_of_species)
    assert len(result) == number_of_species


def test_unique_species():
    number_of_species = 5
    result = draw_random_species(number_of_species)
    assert len(result) == len(set(result))


def test_too_many_species():
    number_of_species = MAX_SPECIES + 1
    with pytest.raises(ValueError):
        draw_random_species(number_of_species)


def test_too_few_species():
    number_of_species = 0
    with pytest.raises(ValueError):
        draw_random_species(number_of_species)


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


# Tests for filter_by_seq_id_and_phylo_dist
def test_no_filters():
    result = filter_by_seq_id_and_phylo_dist()
    pg_overview = get_pg_overview()
    pd.testing.assert_frame_equal(result, pg_overview)


def test_max_phylo_distance_filter():
    max_phylo_distance = 0.9
    pg_overview = get_pg_overview()
    expected = pg_overview[pg_overview["tip2tip_distance"] <= max_phylo_distance]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance)
    pd.testing.assert_frame_equal(result, expected)


def test_min_identity_filter():
    min_identity = 90
    pg_overview = get_pg_overview()
    expected = pg_overview[pg_overview["mean_identity"] >= min_identity]
    result = filter_by_seq_id_and_phylo_dist(min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)


def test_both_filters():
    max_phylo_distance = 0.9
    min_identity = 90
    pg_overview = get_pg_overview()
    expected = pg_overview[(pg_overview["tip2tip_distance"] <= max_phylo_distance) & (pg_overview["mean_identity"] >= min_identity)]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance, min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)


def test_no_matching_records():
    max_phylo_distance = 0.1
    min_identity = 100
    pg_overview = get_pg_overview()
    expected = pg_overview[(pg_overview["tip2tip_distance"] <= max_phylo_distance) & (pg_overview["mean_identity"] >= min_identity)]
    result = filter_by_seq_id_and_phylo_dist(max_phylo_distance=max_phylo_distance, min_identity=min_identity)
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dge_ratio, num_genes", [
    (0.25, 200),
    (0.4, 150),
])
def test_draw_dge_factors_output_shape(dge_ratio, num_genes):
    result = draw_dge_factors(dge_ratio, num_genes)
    assert isinstance(result, np.ndarray)
    assert result.shape[0] == num_genes


@pytest.mark.parametrize("dge_ratio", [-1, -0.01, 1.0, 1.5, 2.0])
def test_draw_dge_factors_ratio_outdie_threshold_raises_value_error(dge_ratio):
    with pytest.raises(ValueError):
        draw_dge_factors(dge_ratio, 100)


def test_draw_dge_factors_valid_ratio_range():
    with pytest.raises(ValueError):
        draw_dge_factors(-0.1, 100)
    with pytest.raises(ValueError):
        draw_dge_factors(1.1, 100)


@pytest.mark.parametrize("dge_ratio", [0.1, 0.25, 0.49])
@pytest.mark.parametrize("num_genes", [10, 100, 1000])
def test_draw_dge_factors_distribution_properties(dge_ratio, num_genes):
    result = draw_dge_factors(dge_ratio, num_genes)
    assert np.all(result > 0)
    assert not np.any(np.isnan(result))


def test_draw_dge_factors_extreme_case():
    result = draw_dge_factors(0.01, 10)
    assert isinstance(result, np.ndarray)
    assert len(result) == 10


def test_draw_dge_factors_empty_output():
    result = draw_dge_factors(0.49, 0)
    assert result.size == 0


def test_calc_zero_ratio_no_zero():
    df = pd.DataFrame({'a': [2, 2, 3], 'b': [4, 3, 6], 'c': [7, 8, 9], 'd': [10, 11, 5]})
    ratio = calc_zero_ratio(df)
    assert ratio == 0


def test_get_all_zero_genes():
    data = {
        "gene_name": ["geneA", "geneB", "geneC"],
        "sample1": [0, 1, 0],
        "sample2": [0, 0, 0],
        "sample3": [0, 2, 0],
    }
    df = pd.DataFrame(data)
    result = get_all_zero_genes(df)
    assert result == {"geneA", "geneC"}, f"Expected {{'geneA', 'geneC'}}, got {result}"


def test_add_counts_to_large_orthogroups():
    data = {
        "gene_name": ["geneA", "geneB", "geneC"],
        "orthogroup": ["og1", "og1", "og2"],
        "sample1": [0, 1, 0],
        "sample2": [0, 0, 0],
        "sample3": [0, 2, 0],
    }
    df = pd.DataFrame(data)

    add_counts_to_large_orthogroups(df, 2)
    all_zero_genes = get_all_zero_genes(df)

    assert all_zero_genes == {"geneC"}, f"Expected to filter {{'geneC'}}, not {all_zero_genes}"


def test_determine_mode_and_model():
    assert dg.determine_mode_and_model(ErrorModel.basic, 150) == (ErrorModel.basic, None)
    assert dg.determine_mode_and_model(ErrorModel.perfect, 100) == (ErrorModel.perfect, None)
    assert dg.determine_mode_and_model(ErrorModel.HiSeq, 150) == ("kde", ErrorModel.HiSeq)
    assert dg.determine_mode_and_model(ErrorModel.NextSeq, None) == ("kde", ErrorModel.NextSeq)
