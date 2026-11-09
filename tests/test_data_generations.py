import pytest
import pandas as pd
import numpy as np
import os
import random
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, draw_dge_factors
from marbel.data_generations import calc_zero_ratio, get_all_zero_genes, add_counts_to_large_orthogroups, select_species_with_criterion, select_orthogroups
from marbel import data_generations as dg
from marbel.presets import MAX_SPECIES, ErrorModel, SelectionCriterion, DESEQ2_FITTED_A0, DESEQ2_FITTED_A1
from marbel.preload import get_pg_overview
import copy

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
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)
    assert len(result) == number_of_orthogroups, f"Expected {number_of_orthogroups}, but got {len(result)}"


def test_max_species_per_group_respected():
    number_of_orthogroups = 20
    max_species = 7
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)
    assert np.all(result <= max_species), f"Some orthogroups exceed the max species of {max_species}"


def test_non_negative_groups():
    number_of_orthogroups = 15
    max_species = 8
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)
    assert np.all(result >= 0), "There are negative values in the orthogroup sizes"


def test_random_seed_reproducibility():
    number_of_orthogroups = 5
    max_species = 10
    rng1 = np.random.default_rng(SEED)
    rng2 = np.random.default_rng(SEED)
    result_1 = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=rng1)
    result_2 = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=rng2)
    np.testing.assert_array_equal(result_1, result_2, err_msg="Results are not identical for the same seed")


def test_different_seed_variation():
    number_of_orthogroups = 5
    max_species = 10
    result_1 = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)
    result_2 = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)
    assert not np.array_equal(result_1, result_2), "Results should differ for different rng draws"


def test_scaled_orthogroups_not_exceed_max():
    number_of_orthogroups = 10
    max_species = 12
    result = create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)

    assert np.max(result) <= max_species, f"Max species per group exceeded: {np.max(result)}"


def test_empty_orthogroups():
    number_of_orthogroups = 0
    max_species = 5
    with pytest.raises(ValueError):
        create_ortholgous_group_rates(number_of_orthogroups, max_species, rng=RNG)


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
    result = draw_dge_factors(dge_ratio, num_genes, rng=RNG)
    assert isinstance(result, np.ndarray)
    assert result.shape[0] == num_genes


@pytest.mark.parametrize("dge_ratio", [-1, -0.01, 1.0, 1.5, 2.0])
def test_draw_dge_factors_ratio_outdie_threshold_raises_value_error(dge_ratio):
    with pytest.raises(ValueError):
        draw_dge_factors(dge_ratio, 100, rng=RNG)


def test_draw_dge_factors_valid_ratio_range():
    with pytest.raises(ValueError):
        draw_dge_factors(-0.1, 100, rng=RNG)
    with pytest.raises(ValueError):
        draw_dge_factors(1.1, 100, rng=RNG)


@pytest.mark.parametrize("dge_ratio", [0.1, 0.25, 0.49])
@pytest.mark.parametrize("num_genes", [10, 100, 1000])
def test_draw_dge_factors_distribution_properties(dge_ratio, num_genes):
    result = draw_dge_factors(dge_ratio, num_genes, rng=RNG)
    assert np.all(result > 0)
    assert not np.any(np.isnan(result))


def test_draw_dge_factors_extreme_case():
    result = draw_dge_factors(0.01, 10, rng=RNG)
    assert isinstance(result, np.ndarray)
    assert len(result) == 10


def test_draw_dge_factors_empty_output():
    result = draw_dge_factors(0.49, 0, rng=RNG)
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


EXP_SPECIES = {
    2449615081: [
        'ALL_Bradyrhizobium_guangxiense_strain_CCBAU_53363_NZ_CP022219',
        'Bradyrhizobium_diazoefficiens_strain_110spc4_NZ_CP032617',
        'Bradyrhizobium_japonicum_strain_5038_NZ_CP058354',
        'Bradyrhizobium_yuanmingense_strain_ALSPC2061_NZ_JARFMK010000009',
        'Bradyrhizobium_arachidis_strain_CCBAU_051107_NZ_CP030050',
        'Bradyrhizobium_zhanjiangense_strain_CCBAU_51778_NZ_CP022221',
        'Bradyrhizobium_sp_ORS_278_strain_ORS278_NC_009445',
        'ALL_Rhodopseudomonas_palustris_strain_RCB100_NZ_CP066699',
        'Candidatus_Karelsulcia_muelleri_strain_PUNC_NZ_CP013212',
        'Liberibacter_crescens_strain_BT_0_NZ_CP010522',
        'Candidatus_Liberibacter_asiaticus_strain_JXGC_NZ_CP019958',
        'Candidatus_Baumannia_cicadellinicola_strain_BGSS_NZ_CP008985',
        'Brucella_canis_strain_2009004498_chromosome_1_NZ_CP016974',
        'Brucella_suis_strain_Human_AR_US_1981_chromosome_I_NZ_CP010851',
        'Brucella_abortus_strain_BFY_chromosome_1_NZ_CP007737',
        'Brucella_melitensis_strain_CIT31_chromosome_1_NZ_CP025822',
        'Brucella_intermedia_strain_NCTC12171_NZ_UGSH01000001',
        'Brucella_haematophila_strain_CCUG_38531_NZ_VCPE01000009',
        'ALL_Agrobacterium_tumefaciens_strain_12D1_chromosome_circular_NZ_CP033032',
        'ALL_Agrobacterium_fabrum_strain_1D132_chromosome_circular_NZ_CP033023',
        'ALL_Agrobacterium_sp_CGMCC_11546_chromosome_1_NZ_CP047388',
        'ALL_Agrobacterium_rosae_strain_NCPPB_1650_NZ_NXEJ01000009',
        'Rhizobium_rhizophilum_strain_7209_2_NZ_STGT01000009',
        'ALL_Agrobacterium_vitis_strain_NCPPB_3554_NZ_LMVL02000008',
        'ALL_Rhizobium_rhizogenes_strain_LBA9402_chromosome_1_NZ_CP044123',
        'ALL_Rhizobium_tropici_strain_CB3060_NZ_CP104158',
        'ALL_Rhizobium_etli_strain_NXC12_NZ_CP020906',
        'ALL_Rhizobium_phaseoli_strain_R650_NZ_CP013532',
        'ALL_Rhizobium_leguminosarum_strain_SM52_NZ_SILH01000009',
        'Sinorhizobium_medicae_strain_USDA1037_NZ_VITA01000068'],
    2802522352: [
        'ALL_Rhodococcus_opacus_strain_S8_NZ_CP093380',
        'Rhodococcus_jostii_strain_DSM_44719_NZ_FNTL01000001',
        'ALL_Rhodococcus_oxybenzonivorans_strain_S2_17_NZ_CP021354',
        'ALL_Rhodococcus_aetherivorans_strain_N1_NZ_CP106982',
        'ALL_Rhodococcus_ruber_strain_C1_NZ_CP044211',
        'ALL_Rhodococcus_rhodochrous_strain_LH_B3_NZ_CP120356',
        'ALL_Rhodococcus_sp_p52_NZ_CP016819',
        'ALL_Prescottella_equi_strain_DSSKP_R_001_NZ_CP027793',
        'ALL_Rhodococcus_qingshengii_strain_CS98_NZ_AP023172',
        'Gordonia_terrae_strain_NCTC10669_NZ_CABEIC010000002',
        'ALL_Gordonia_rubripertincta_strain_SD5_NZ_CP059694',
        'Dietzia_natronolimnaea_strain_S_XJ_1_NZ_NTGA01000099',
        'Rhodoluna_lacicola_strain_MWH_Ta8_NZ_CP007490',
        'Candidatus_Karelsulcia_muelleri_strain_PUNC_NZ_CP013212',
        'Candidatus_Baumannia_cicadellinicola_strain_BGSS_NZ_CP008985',
        'Corynebacterium_jeikeium_strain_DSM_7171_NZ_CP063195',
        'Corynebacterium_diphtheriae_strain_BQ11_NZ_CP029644',
        'Corynebacterium_glutamicum_strain_ATCC_13032_NZ_CP115148',
        'ALL_Micrococcus_luteus_strain_CW_Ay_NZ_CP082331',
        'ALL_Glutamicibacter_nicotianae_strain_OTC_16_NZ_CP033081',
        'ALL_Microbacterium_sediminis_strain_YLB_01_NZ_CP038256',
        'Dermacoccus_abyssi_strain_MT1_1_NZ_JAVYJM010000004',
        'ALL_Mycolicibacterium_smegmatis_strain_FDAARGOS_679_NZ_CP054795',
        'Mycolicibacterium_fortuitum_strain_MF_GZ001_NZ_CP107719',
        'Mycolicibacterium_iranicum_strain_DSM_45541_NZ_LQPC01000082',
        'ALL_Mycobacteroides_abscessus_strain_GZ002_NZ_CP034181',
        'Mycobacteroides_chelonae_strain_MCHE08_NZ_CP058976',
        'Mycobacterium_canetti_strain_Percy1130_chromosome_McanPercy1130_NZ_OU917819',
        'Mycobacterium_tuberculosis_strain_LE13_NZ_CP023602',
        'Mycobacterium_marinum_strain_050012_NZ_CP118910'],
    3813447723: [
        'Leptospira_yasudae_strain_F1_NZ_QHCU01000009',
        'Leptospira_tipperaryensis_strain_GWTS1_chromosome_1_NZ_CP015218',
        'Leptospira_adleri_strain_201601302_NZ_RQHG01000008',
        'Leptospira_gomenensis_strain_201800299_NZ_RQFA01000093',
        'ALL_Leptospira_weilii_strain_CUDO6_chromosome_I_NZ_CP040841',
        'Candidatus_Karelsulcia_muelleri_strain_PUNC_NZ_CP013212',
        'Candidatus_Baumannia_cicadellinicola_strain_BGSS_NZ_CP008985',
        'Ehrlichia_japonica_strain_HF_NZ_CP007474',
        'Wolbachia_pipientis_strain_wStv_NZ_CP050531',
        'Candidatus_Liberibacter_asiaticus_strain_JXGC_NZ_CP019958',
        'Liberibacter_crescens_strain_BT_0_NZ_CP010522',
        'Haemophilus_ducreyi_strain_VAN2_NZ_CP015425',
        'Haemophilus_influenzae_strain_477_NZ_CP007470',
        'Haemophilus_parainfluenzae_strain_EL1_NZ_CP125087',
        'Aggregatibacter_actinomycetemcomitans_strain_4S_NZ_CP085096',
        'ALL_Pasteurella_multocida_strain_FDAARGOS_218_NZ_CP020405',
        'Actinobacillus_succinogenes_strain_GXAS137_NZ_NHRD01000002',
        'ALL_Glaesserella_parasuis_strain_YHP1818_NZ_CP071487',
        'ALL_Mannheimia_haemolytica_strain_USDA_ARS_USMARC_191_NZ_CP023044',
        'Avibacterium_paragallinarum_strain_ESV_135_NZ_CP050316',
        'Neisseria_meningitidis_strain_11_7_NZ_CP021520',
        'ALL_Neisseria_gonorrhoeae_strain_TUM19854_NZ_AP023069',
        'ALL_Conchiformibius_steedae_strain_DSM_2580_NZ_CP059563',
        'ALL_Moraxella_catarrhalis_strain_CCRI_195ME_NZ_CP018059',
        'Candidatus_Thioglobus_autotrophicus_strain_EF1_NZ_CP010552',
        'Mesoplasma_florum_strain_MouA_2_NZ_CP022508',
        'Mycoplasma_mycoides_strain_Caprivi_NZ_AP024602',
        'ALL_Edwardsiella_tarda_strain_KC_Pc_HB1_NZ_CP023706',
        'ALL_Edwardsiella_hoshinae_strain_FDAARGOS_940_NZ_CP065626',
        'ALL_Edwardsiella_piscicida_strain_18EpOKYJ_NZ_CP090968']
 }


@pytest.mark.parametrize("seed", list(EXP_SPECIES.keys()))
def test_speedup_select_species_with_criterion(seed):
    rng = np.random.default_rng(seed)
    obs = select_species_with_criterion(number_of_species=30, rng=rng, selection_criterion=SelectionCriterion.maximize)
    assert EXP_SPECIES[seed] == obs


@pytest.mark.parametrize("seed", list(EXP_SPECIES.keys()))
def test_select_species_rng(seed):
    rng1 = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed)
    initial_state = copy.deepcopy(rng1.bit_generator.state)

    obs1 = select_species_with_criterion(
        number_of_species=30,
        rng=rng1,
        selection_criterion=SelectionCriterion.maximize,
    )

    obs2 = select_species_with_criterion(
        number_of_species=30,
        rng=rng2,
        selection_criterion=SelectionCriterion.maximize,
    )

    assert obs1 == obs2

    # RNG state was consumed
    assert rng1.bit_generator.state != initial_state


@pytest.mark.parametrize("seed", list(EXP_SPECIES.keys()))
def test_select_species_rng(seed):
    rng1 = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed)
    initial_state = copy.deepcopy(rng1.bit_generator.state)

    obs1 = select_species_with_criterion(
        number_of_species=30,
        rng=rng1,
        selection_criterion=SelectionCriterion.minimize,
    )

    obs2 = select_species_with_criterion(
        number_of_species=30,
        rng=rng2,
        selection_criterion=SelectionCriterion.minimize,
    )

    assert obs1 == obs2

    # RNG state was consumed
    assert rng1.bit_generator.state != initial_state


@pytest.mark.parametrize("seed", list(EXP_SPECIES.keys()))
def test_select_orthogroups(seed):
    exp = pd.read_csv(os.path.dirname(__file__) + '/data/orthogroups_high.tsv', sep="\t").set_index(['Seed', 'axis'])
    rng = np.random.default_rng(seed)

    filtered_orthog_groups = dg.filter_by_seq_id_and_phylo_dist(None, None)
    obs = select_orthogroups(filtered_orthog_groups, EXP_SPECIES[seed], 10000, rng, minimize=False, force=True)

    delta = 0.001
    print(obs)
    assert np.isclose(exp.loc[seed, 'species'], (obs.iloc[:, :-1] != '-').sum().describe(), atol=delta).all()
    assert np.isclose(exp.loc[seed, 'genes'], (obs.iloc[:, :-1] != '-').sum(axis=1).describe(), atol=delta).all()


def test_odd_DGE_distribution():
    seed = 2449615081
    species = EXP_SPECIES[seed]
    filtered_orthog_groups = dg.filter_by_seq_id_and_phylo_dist(None, None)
    rng = np.random.default_rng(seed)
    selected_ortho_groups = select_orthogroups(filtered_orthog_groups, EXP_SPECIES[seed], 10000, rng, minimize=False, force=True)
    species_abundances = dg.generate_species_abundance(len(EXP_SPECIES[seed]), rng)
    number_of_selected_genes = selected_ortho_groups["group_size"].sum()
    read_mean_counts = dg.generate_read_mean_counts(number_of_selected_genes, rng)
    gene_summary_df = dg.aggregate_gene_data(species, species_abundances, selected_ortho_groups, read_mean_counts)

    # issue is in this random number generation
    dge_factors = draw_dge_factors(0.1, number_of_selected_genes, rng)

    gene_summary_df["simulation_fold_change"] = dge_factors
    sample_group_1 = dg.create_sample_values(gene_summary_df, 10, True, DESEQ2_FITTED_A0, DESEQ2_FITTED_A1, rng)
    gene_summary_df = pd.merge(gene_summary_df, sample_group_1, on="gene_name")
    sample_group_2 = dg.create_sample_values(gene_summary_df, 10, False, DESEQ2_FITTED_A0, DESEQ2_FITTED_A1, rng)
    gene_summary_df = pd.merge(gene_summary_df, sample_group_2, on="gene_name")

    # for diagnosis
    gene_summary_df['counts_group1'] = gene_summary_df[[c for c in gene_summary_df.columns if 'group_1_sample_' in c]].sum(axis=1)
    gene_summary_df['counts_group2'] = gene_summary_df[[c for c in gene_summary_df.columns if 'group_2_sample_' in c]].sum(axis=1)
    gene_summary_df['counts_all'] = gene_summary_df[['counts_group1', 'counts_group2']].sum(axis=1)
    nonzero_genes = gene_summary_df[gene_summary_df['counts_all'] > 0]
    simulated_up_regulated_genes = sum(nonzero_genes["simulation_fold_change"] >= 2.0)
    simulated_down_regulated_genes = sum(nonzero_genes["simulation_fold_change"] <= 0.5)

    assert abs((simulated_up_regulated_genes / (simulated_up_regulated_genes + simulated_down_regulated_genes)) - 0.5) < 0.1, 'odd distribution of non-negative DGE genes'


def test_reproducibility_draw_dge_factors():
    seed = 2449615081
    exp = [1.033994979662841, 1.6805846532613182, 1.147713164635003, 
           0.6923367274605421, 2.5190819560518323, 0.9591528484220874]
    rng = np.random.default_rng(seed)
    obs = draw_dge_factors(0.1, 10000, rng)
    assert np.allclose(list(obs)[:3], exp[:3])
    assert np.allclose(list(obs)[-3:], exp[-3:])