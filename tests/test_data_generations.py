import pytest
import pandas as pd
import numpy as np
import random
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, draw_dge_factors
from marbel.data_generations import calc_zero_ratio, get_all_zero_genes, add_counts_to_large_orthogroups, select_species_with_criterion
from marbel import data_generations as dg
from marbel.presets import MAX_SPECIES, ErrorModel, SelectionCriterion
from marbel.preload import get_pg_overview

# Tests for draw_random_species
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
RNG = np.random.default_rng(SEED)


def test_correct_length():
    number_of_species = 5
    result = draw_random_species(number_of_species, rng=RNG)
    assert len(result) == number_of_species


def test_unique_species():
    number_of_species = 5
    result = draw_random_species(number_of_species, rng=RNG)
    assert len(result) == len(set(result))


def test_too_many_species():
    number_of_species = MAX_SPECIES + 1
    with pytest.raises(ValueError):
        draw_random_species(number_of_species, rng=RNG)


def test_too_few_species():
    number_of_species = 0
    with pytest.raises(ValueError):
        draw_random_species(number_of_species, rng=RNG)


def test_correct_number_of_orthogroups():
    number_of_orthogroups = 10
    max_species = 5
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)
    assert len(result) == number_of_orthogroups, f"Expected {number_of_orthogroups}, but got {len(result)}"


def test_max_species_per_group_respected():
    number_of_orthogroups = 20
    max_species = 7
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)
    assert np.all(result <= max_species), f"Some orthogroups exceed the max species of {max_species}"


def test_non_negative_groups():
    number_of_orthogroups = 15
    max_species = 8
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)
    assert np.all(result >= 0), "There are negative values in the orthogroup sizes"


def test_random_seed_reproducibility():
    number_of_orthogroups = 5
    max_species = 10
    result_1 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)
    result_2 = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)
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
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)

    assert np.max(result) <= max_species, f"Max species per group exceeded: {np.max(result)}"


def test_empty_orthogroups():
    number_of_orthogroups = 0
    max_species = 5
    with pytest.raises(ValueError):
        create_ortholgous_group_rates(number_of_orthogroups, max_species, seed=SEED)


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
    result = draw_dge_factors(dge_ratio, num_genes, seed=SEED)
    assert isinstance(result, np.ndarray)
    assert result.shape[0] == num_genes


@pytest.mark.parametrize("dge_ratio", [-1, -0.01, 1.0, 1.5, 2.0])
def test_draw_dge_factors_ratio_outdie_threshold_raises_value_error(dge_ratio):
    with pytest.raises(ValueError):
        draw_dge_factors(dge_ratio, 100, seed=SEED)


def test_draw_dge_factors_valid_ratio_range():
    with pytest.raises(ValueError):
        draw_dge_factors(-0.1, 100, seed=SEED)
    with pytest.raises(ValueError):
        draw_dge_factors(1.1, 100, seed=SEED)


@pytest.mark.parametrize("dge_ratio", [0.1, 0.25, 0.49])
@pytest.mark.parametrize("num_genes", [10, 100, 1000])
def test_draw_dge_factors_distribution_properties(dge_ratio, num_genes):
    result = draw_dge_factors(dge_ratio, num_genes, seed=SEED)
    assert np.all(result > 0)
    assert not np.any(np.isnan(result))


def test_draw_dge_factors_extreme_case():
    result = draw_dge_factors(0.01, 10, seed=SEED)
    assert isinstance(result, np.ndarray)
    assert len(result) == 10


def test_draw_dge_factors_empty_output():
    result = draw_dge_factors(0.49, 0, seed=SEED)
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

    add_counts_to_large_orthogroups(df, 2, rng=RNG)
    all_zero_genes = get_all_zero_genes(df)

    assert all_zero_genes == {"geneC"}, f"Expected to filter {{'geneC'}}, not {all_zero_genes}"


def test_determine_mode_and_model():
    assert dg.determine_mode_and_model(ErrorModel.basic, 150) == (ErrorModel.basic, None)
    assert dg.determine_mode_and_model(ErrorModel.perfect, 100) == (ErrorModel.perfect, None)
    assert dg.determine_mode_and_model(ErrorModel.HiSeq, 150) == ("kde", ErrorModel.HiSeq)
    assert dg.determine_mode_and_model(ErrorModel.NextSeq, None) == ("kde", ErrorModel.NextSeq)


def test_speedup_select_species_with_criterion():
    exp_2449615081 = [
        'ALL_Companilactobacillus_futsaii_strain_Y97_NZ_CP040736',
        'Latilactobacillus_curvatus_strain_WiKim52_NZ_CP016602',
        'Ligilactobacillus_agilis_strain_AM_LB8_NZ_CP104397',
        'ALL_Ligilactobacillus_salivarius_strain_LPM01_chromosome_I_NZ_LT604074',
        'ALL_Ligilactobacillus_faecis_strain_WILCCON_0062_NZ_CP123639',
        'ALL_Liquorilactobacillus_mali_strain_LM596_NZ_CP045035',
        'ALL_Pediococcus_claussenii_strain_TMW_2_53_NZ_CP014933',
        'Pediococcus_acidilactici_strain_PMC65_NZ_CP053421',
        'Apilactobacillus_kunkeei_strain_DSMZ_12361_NZ_CP080568',
        'Limosilactobacillus_reuteri_strain_2010_NZ_PUXG01000029',
        'ALL_Limosilactobacillus_fermentum_strain_EFEL6800_NZ_CP124737',
        'ALL_Lactiplantibacillus_paraplantarum_strain_FL_8_NZ_CP118745',
        'Lactiplantibacillus_plantarum_strain_SRCM100442_NZ_CP028221',
        'ALL_Lactiplantibacillus_pentosus_strain_DSM_20314_NZ_CP032757',
        'ALL_Levilactobacillus_brevis_strain_LMT1_73_NZ_CP033885',
        'ALL_Lentilactobacillus_kefiri_strain_DH5_NZ_CP029971',
        'Lactobacillus_acetotolerans_strain_LA749_NZ_CP044496',
        'Lactobacillus_johnsonii_strain_NCK2677_NZ_CP059055',
        'Lactobacillus_acidophilus_strain_HN017_NZ_CP101537',
        'ALL_Lactobacillus_amylolyticus_strain_L5_NZ_CP031835',
        'Lactobacillus_helveticus_strain_DS3_8_NZ_QAZU01000009',
        'Lacticaseibacillus_rhamnosus_strain_1_0320_NZ_CP040780',
        'Streptococcus_equinus_strain_pGA_7_NZ_FNJX01000001',
        'ALL_Streptococcus_mutans_strain_NCH105_NZ_CP044221',
        'ALL_Streptococcus_agalactiae_strain_NGBS128_NZ_CP012480',
        'Streptococcus_pyogenes_strain_NCTC12064_chromosome_1_NZ_LS483338',
        'Streptococcus_suis_strain_NCTC10234_chromosome_1_NZ_LS483418',
        'Streptococcus_anginosus_strain_J4211_NZ_CP012805',
        'Streptococcus_mitis_strain_OT25_NZ_JYGP01000003',
        'Streptococcus_oralis_strain_FDAARGOS_1075_NZ_CP066041']
    obs = select_species_with_criterion(number_of_species=30, seed=2449615081, selection_criterion=SelectionCriterion.maximize)
    assert exp_2449615081 == obs

    exp_2802522352 = [
        'ALL_Staphylococcus_haemolyticus_strain_S167_NZ_CP013911',
        'ALL_Staphylococcus_hominis_strain_FDAARGOS_575_NZ_CP033732',
        'ALL_Staphylococcus_capitis_strain_FDAARGOS_378_NZ_CP023966',
        'ALL_Staphylococcus_epidermidis_strain_ATCC_14990_NZ_CP035288',
        'Staphylococcus_simiae_strain_NCTC13838_chromosome_1_NZ_LT906460',
        'ALL_Staphylococcus_saprophyticus_strain_1A_NZ_CP031196',
        'Staphylococcus_xylosus_strain_SMQ_121_NZ_CP008724',
        'ALL_Staphylococcus_equorum_strain_FDAARGOS_1068_NZ_CP066013',
        'ALL_Staphylococcus_saccharolyticus_strain_13T0028_NZ_CP068029',
        'Staphylococcus_carnosus_strain_NCTC13825_NZ_UHCT01000001',
        'ALL_Staphylococcus_pseudintermedius_strain_SP_11304_3A_NZ_CP065921',
        'Mammaliicoccus_lentus_strain_Colony453_NZ_CP075503',
        'Staphylococcus_aureus_strain_Sa_Newman_UoM_chromosome_1_NZ_LT598688',
        'Candidatus_Karelsulcia_muelleri_strain_PUNC_NZ_CP013212',
        'Candidatus_Baumannia_cicadellinicola_strain_BGSS_NZ_CP008985',
        'Apilactobacillus_kunkeei_strain_DSMZ_12361_NZ_CP080568',
        'Streptococcus_equinus_strain_pGA_7_NZ_FNJX01000001',
        'ALL_Streptococcus_thermophilus_strain_STH_CIRM_65_chromosome_STHERMO_NZ_LR822015',
        'Streptococcus_pyogenes_strain_NCTC12064_chromosome_1_NZ_LS483338',
        'ALL_Streptococcus_agalactiae_strain_NGBS128_NZ_CP012480',
        'Streptococcus_suis_strain_NCTC10234_chromosome_1_NZ_LS483418',
        'Streptococcus_anginosus_strain_J4211_NZ_CP012805',
        'Streptococcus_oralis_strain_FDAARGOS_1075_NZ_CP066041',
        'Streptococcus_mitis_strain_OT25_NZ_JYGP01000003',
        'Streptococcus_downii_strain_CECT_9732_NZ_MZNX01000009',
        'Streptococcus_gordonii_strain_NCTC10231_chromosome_1_NZ_LR594049',
        'Streptococcus_sinensis_strain_HKU4_NZ_JPEN01000059',
        'ALL_Streptococcus_pneumoniae_strain_Hu17_NZ_CP020549',
        'Streptococcus_sanguinis_strain_NCTC11085_chromosome_1_NZ_LS483346',
        'ALL_Streptococcus_mutans_strain_NCH105_NZ_CP044221']
    obs = select_species_with_criterion(number_of_species=30, seed=2802522352, selection_criterion=SelectionCriterion.maximize)
    assert exp_2802522352 == obs

    exp_3813447723 = [
        'Francisella_tularensis_strain_Schu_S4_249_NZ_CP073129',
        'Candidatus_Baumannia_cicadellinicola_strain_BGSS_NZ_CP008985',
        'Haemophilus_ducreyi_strain_VAN2_NZ_CP015425',
        'Haemophilus_influenzae_strain_477_NZ_CP007470',
        'Haemophilus_parainfluenzae_strain_EL1_NZ_CP125087',
        'Aggregatibacter_actinomycetemcomitans_strain_4S_NZ_CP085096',
        'ALL_Pasteurella_multocida_strain_FDAARGOS_218_NZ_CP020405',
        'Actinobacillus_succinogenes_strain_GXAS137_NZ_NHRD01000002',
        'ALL_Glaesserella_parasuis_strain_YHP1818_NZ_CP071487',
        'ALL_Mannheimia_haemolytica_strain_USDA_ARS_USMARC_191_NZ_CP023044',
        'Avibacterium_paragallinarum_strain_ESV_135_NZ_CP050316',
        'Candidatus_Karelsulcia_muelleri_strain_PUNC_NZ_CP013212',
        'Ehrlichia_japonica_strain_HF_NZ_CP007474',
        'Wolbachia_pipientis_strain_wStv_NZ_CP050531',
        'Neisseria_meningitidis_strain_11_7_NZ_CP021520',
        'ALL_Neisseria_gonorrhoeae_strain_TUM19854_NZ_AP023069',
        'ALL_Conchiformibius_steedae_strain_DSM_2580_NZ_CP059563',
        'ALL_Moraxella_catarrhalis_strain_CCRI_195ME_NZ_CP018059',
        'Candidatus_Thioglobus_autotrophicus_strain_EF1_NZ_CP010552',
        'Candidatus_Liberibacter_asiaticus_strain_JXGC_NZ_CP019958',
        'Liberibacter_crescens_strain_BT_0_NZ_CP010522',
        'ALL_Edwardsiella_tarda_strain_KC_Pc_HB1_NZ_CP023706',
        'ALL_Edwardsiella_hoshinae_strain_FDAARGOS_940_NZ_CP065626',
        'ALL_Edwardsiella_piscicida_strain_18EpOKYJ_NZ_CP090968',
        'ALL_Yersinia_ruckeri_strain_KMM821_NZ_CP071802',
        'Yersinia_entomophaga_strain_MH96_NZ_CP010029',
        'ALL_Yersinia_enterocolitica_strain_NW56_NZ_CP107102',
        'ALL_Hafnia_alvei_strain_A23BA_NZ_CP050150',
        'Enterobacter_sp_S_33_NZ_CP117838',
        'ALL_Enterobacter_asburiae_strain_17Nkhm_UP2_NZ_AP019630']
    obs = select_species_with_criterion(number_of_species=30, seed=3813447723, selection_criterion=SelectionCriterion.maximize)
    assert exp_3813447723 == obs
