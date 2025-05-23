import numpy as np
import pandas as pd
from scipy import stats
from Bio import SeqIO
import os
import random
import argparse
from insilicoseq_marbeldep.app import generate_reads
from joblib import Parallel, delayed
import polars as pl
import pymc as pm
import sys

from marbel.presets import MAX_SPECIES, PATH_TO_GROUND_GENES_INDEX, DGE_LOG_2_CUTOFF_VALUE, PANGENOME_OVERVIEW
from marbel.presets import ErrorModel, LibrarySizeDistribution, __version__, SelectionCriterion

from marbel.preload import get_pg_overview, get_species_tree, get_species_stats_dict, get_pymc_model


def draw_random_species(number_of_species):
    """
    Draws a random sample of species from the list of available species set in the presets.

    Parameters:
    number_of_species (int): The number of species to draw.

    Returns:
    list: A list of randomly drawn species from the list of available species.
    """
    if number_of_species < 1 or number_of_species > MAX_SPECIES:
        raise ValueError(f"Number of species must be between 1 and {MAX_SPECIES}.")
    available_species = pl.read_parquet(PANGENOME_OVERVIEW).columns[:MAX_SPECIES]
    return random.sample(available_species, number_of_species)


def create_ortholgous_group_rates(number_of_orthogous_groups, max_species_per_group, seed=None):
    """
    Creates a list of group sizes for orthogroups, such that the maximum group size is less than or equal to
    the specified maximum species per group and the total number of orthogroups matches the specified number
    of orthogroups.

    Parameters:
        number_of_orthogous_groups (int): The total number of orthogroups.
        max_species_per_group (int): The maximum number of species per orthogroup.
        seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
        numpy.ndarray: A numpy array containing the group size for each orthogroup.
    """
    if number_of_orthogous_groups < 1:
        raise ValueError("Number of orthogroups must be greater than 0.")
    model = get_pymc_model()
    with model:
        orthologues_samples = pm.sample_prior_predictive(number_of_orthogous_groups, var_names=['ortho'], random_seed=seed)
    orthogroups = orthologues_samples.to_dataframe()["ortho"].to_numpy()
    scaled_orthogroups = orthogroups * (max_species_per_group / np.max(orthogroups))
    scaled_orthogroups = np.minimum(np.ceil(scaled_orthogroups), max_species_per_group)
    return scaled_orthogroups


def filter_by_seq_id_and_phylo_dist(max_phylo_distance=None, min_identity=None):
    """
    Filters the orthologous group overview dataframe by both the mean phylogenetic distance to the root of the tree
    and the mean sequence identity if they are specified.

    Parameters:
        max_phylo_distance (float, optional): The maximum median phylogenetic distance to the root of the tree. Defaults to None.
        min_identity (float, optional): The minimum median sequence identity. Defaults to None.

    Returns:
        pandas.DataFrame: The filtered orthologous group overview dataframe.
    """
    pg_overview = get_pg_overview()
    if max_phylo_distance is not None:
        filtered_slice = pg_overview[pg_overview["tip2tip_distance"] <= max_phylo_distance]
    else:
        filtered_slice = pg_overview
    if min_identity is not None:
        filtered_slice = filtered_slice[filtered_slice["mean_identity"] >= min_identity]
    if max_phylo_distance is None and min_identity is None:
        return pg_overview
    return filtered_slice


# randomization based on rates calculated from the pdf
def draw_orthogroups_by_rate(orthogroup_slice, orthogroup_rates, species):
    """
    Draws orthologous groups based on their rates. Given a dataframe slice of orthologous groups,
    a list of rates for each orthologous group, and a list of species, randomly samples orthologous groups
    based on their rates. If there are not enough orthologous groups with a given rate, it returns None.

    Parameters:
        orthogroup_slice (pandas.DataFrame): The orthologous group dataframe.
        orthogroup_rates (list): A list of rates for each orthologous group.
        species (list): A list of species.

    Returns:
        pandas.DataFrame or None: A randomly sampled dataframe of orthologous groups based on their rates and filtered by species,
        or None if there are not enough orthologous groups with a given rate.
    """
    orthogroup_rates = pd.Series(orthogroup_rates)
    value_counts = orthogroup_rates.value_counts()
    orthogroups = orthogroup_slice[species].copy()
    number_of_species = len(species)
    orthogroups["group_size"] = orthogroups.apply(lambda x: number_of_species - len(x[x == "-"]), axis=1)
    orthogroups = orthogroups[orthogroups["group_size"] > 0]
    sampled_groups = pd.DataFrame(columns=orthogroups.columns)
    for index in value_counts.index:
        subset_orthogroups = orthogroups[orthogroups["group_size"] == index]
        if len(subset_orthogroups) < value_counts[index]:
            print(f"Warning: Not enough orthogroups with size {index} to satisfy the rate, switching to sampling the orthogroups without the gamma function rates.")
            return None
        subset_orthogroups = subset_orthogroups.sample(value_counts[index])
        sampled_groups = pd.concat([sampled_groups, subset_orthogroups])
    return sampled_groups


def draw_orthogroups(orthogroup_slice, number_of_orthogous_groups, species, force):
    """
    Draws orthologous groups based on actual occurences in the dataset, instead of rates based on the pdf.
    Given a dataframe slice of orthologous groups, a number of orthologous groups to be drawn, and a list of species,
    randomly samples orthologous groups based on their actual occurences. If there are not enough orthologous groups
    to satisfy the parameters, it raises an error.

    Parameters:
        orthogroup_slice (pandas.DataFrame): The orthologous group dataframe.
        number_of_orthogous_groups (int): The number of orthologous groups to be drawn.
        species (list): A list of species.
        force (bool): If True, returns all available orthogroups when there aren't enough to satisfy the request.
                      If False, exits with an error message when there aren't enough orthogroups.

    Returns:
        pandas.DataFrame: A randomly sampled dataframe of orthologous groups based on their actual occurences and filtered by species.
    """
    orthogroups = orthogroup_slice[species].copy()
    orthogroups["group_size"] = orthogroups.apply(lambda x: len(species) - len(x[x == "-"].index.to_list()), axis=1)
    orthogroups = orthogroups[orthogroups["group_size"] > 0]

    check_force_result = check_enough_orthogroups(orthogroups, number_of_orthogous_groups, force)
    if check_force_result is not None:
        return check_force_result

    orthogroups_sample = orthogroups.sample(n=number_of_orthogous_groups)
    return orthogroups_sample


def check_enough_orthogroups(orthogroups, required_count, force):
    if orthogroups.shape[0] < required_count:
        print("Not enough orthogroups to satisfy the parameters.")
        print("Try increasing the number of species, relaxing filters, or lowering the count.")
        print(f"Available: {orthogroups.shape[0]}, Required: {required_count}")

        if force:
            print("Returning all available orthogroups due to --force flag.")
            return orthogroups
        else:
            print("Exiting. Use --force to override.")
            sys.exit(1)
    return None


def generate_species_abundance(number_of_species, seed=None):
    """
    Generates a list of species abundances based on the distribution in the preset model.

    Parameters:
    number_of_species (int): The number of species for which to generate abundances.
    seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
    list: A list of species abundances.
    """
    model = get_pymc_model()
    with model:
        species_samples = pm.sample_prior_predictive(number_of_species, var_names=['species'], random_seed=seed)
    return species_samples.to_dataframe()['species'].to_list()


def generate_read_mean_counts(number_of_reads, seed=None):
    """
    Generates a list of read mean counts based on the distribution in the preset model.

    Parameters:
    number_of_reads (int): The number of reads for which to generate counts.
    seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
    list: A list of read mean counts.
    """
    model = get_pymc_model()
    with model:
        reads = pm.sample_prior_predictive(number_of_reads, var_names=['read_counts'], random_seed=seed)
    return reads.to_dataframe()["read_counts"].to_list()


def filter_genes_from_ground(gene_list, output_fasta, gtf_path):
    """
    Writes genes from a list to a FASTA file. Duplicates will also be written, the order is preserved.
    """
    ground_genes = SeqIO.index_db(PATH_TO_GROUND_GENES_INDEX)
    records = []
    gtf_entries = []

    for gene in gene_list:
        try:
            record = ground_genes[gene]
            records.append(record)
            gtf_entries.append({
                "cds": gene,
                "block_start": 1,
                "end": len(record.seq),
            })
        except KeyError:
            print(f"Critical Warning: Gene {gene} not found in the ground genes.")
    SeqIO.write(records, output_fasta, "fasta")

    blocks_df = pl.DataFrame(gtf_entries).with_columns(
        pl.lit("marbel").alias("source"),
        pl.lit("gene").alias("feature"),
        pl.lit(".").alias("score"),
        pl.lit("+").alias("strand"),
        pl.lit(".").alias("frame"),
        pl.format(
            "gene_id {}; transcript_id {};",
            pl.col("cds"), pl.col("cds")
        ).alias("attributes"),
        (pl.col("block_start") + 1).alias("gtf_start"),
    )

    blocks_df.select([
        "cds", "source", "feature", "gtf_start", "end", "score", "strand", "frame", "attributes"
    ]).write_csv(gtf_path, separator="\t", include_header=False)


def aggregate_gene_data(species, species_abundances, selected_ortho_groups, read_mean_counts):
    """
    Aggregates gene names, species, read mean counts, base read mean counts, species abundances and orthologous groups from given data.

    Parameters:
    species (list): A list of species.
    species_abundances (list): A list of species abundances.
    selected_ortho_groups (pandas.DataFrame): A dataframe of orthologous groups.
    read_mean_counts (list): A list of read mean counts.

    Returns:
    pandas.DataFrame: A dataframe with gene names, species, read mean counts, base read mean counts, species abundances and orthologous groups.
    """
    gene_summary_df = pd.DataFrame()
    scaled_read_mean_counts, all_species_genes = [], []
    current_read_index, i = 0, 0
    species_weights = species_abundances / np.sum(species_abundances)
    species_weight_col = []
    origin_species = []
    origin_orthogroup = []

    for sp in species:
        species_genes_list = selected_ortho_groups[selected_ortho_groups[sp] != "-"][sp].to_list()
        origin_orthogroup += [f"og{ortho_group}" for ortho_group in selected_ortho_groups[selected_ortho_groups[sp] != "-"].index.to_list()]
        scaled_read_mean_counts += [species_weights[i] * c for c in read_mean_counts[current_read_index:(current_read_index + len(species_genes_list))]]
        current_read_index += len(species_genes_list)
        all_species_genes += species_genes_list
        origin_species += [sp] * len(species_genes_list)
        species_weight_col += [species_weights[i]] * len(species_genes_list)
        i += 1

    gene_summary_df["gene_name"] = all_species_genes
    gene_summary_df["origin_species"] = origin_species
    gene_summary_df["read_mean_count"] = scaled_read_mean_counts
    gene_summary_df["base_read_mean"] = read_mean_counts
    gene_summary_df["species_abundance"] = species_weight_col
    gene_summary_df["orthogroup"] = origin_orthogroup

    return gene_summary_df


def write_parameter_summary(number_of_orthogous_groups, number_of_species, number_of_sample, outdir, max_phylo_distance,
                            min_identity, deg_ratio, seed, output_format, error_model, read_length, library_size, library_distribution, library_sizes, min_sparsity,
                            force, actual_orthogroups, min_overlap, summary_dir):
    """
    Writes the simulation parameters to the result_file.

    Args:
        number_of_orthogous_groups (int): The number of orthologous groups.
        number_of_species (int): The number of species.
        number_of_sample (tuple): The number of samples (group 1, group 2).
        outdir (str): The output directory.
        max_phylo_distance (float): The maximum phylogenetic distance.
        min_identity (float): The minimum sequence identity.
        deg_ratio (float): The ratio of up and down regulated genes.
        seed (int): The seed for the simulation.
        compressed (bool): Compression of files.
        read_length (int): The read length.
        result_file (file): The file to write the summary to.
    """
    with open(f"{summary_dir}/marbel_params.txt", "w") as result_file:
        result_file.write(f"Marbel version: {__version__}\n")
        result_file.write(f"Number of orthogroups: {number_of_orthogous_groups}\n")
        result_file.write(f"Number of species: {number_of_species}\n")
        result_file.write(f"Number of samples: {number_of_sample}\n")
        result_file.write(f"Output directory: {outdir}\n")
        result_file.write(f"Max phylogenetic distance: {max_phylo_distance}\n")
        result_file.write(f"Min identity: {min_identity}\n")
        result_file.write(f"Ratio of up and down regulated genes: {deg_ratio}\n")
        result_file.write(f"Seed: {seed}\n")
        result_file.write(f"File compression: {output_format}\n")
        result_file.write(f"Model used: {error_model}\n")
        result_file.write(f"Read length: {read_length}\n")
        result_file.write(f"Library size: {library_size}\n")
        result_file.write(f"Library size distribution: {library_distribution}\n")
        result_file.write(f"Library sizes for samples: {library_sizes}\n")
        result_file.write(f"Minimum sparsity: {min_sparsity}\n")
        result_file.write(f"Forced creation: {force}\n")
        result_file.write(f"Actual orthogroups (if force was used): {actual_orthogroups}\n")
        result_file.write(f"Minimum overlap for Overlap Blocks: {min_overlap}\n")


def generate_report(summary_dir, gene_summary, number_of_dropped_genes, specified_orthogroups):
    """
    Generates a report of the simulation parameters.

    Parameters:
        summary_dir (str): The output directory for the summary
        gene_summary (pandas.DataFrame): The summary of genes.
    """
    species_tree = get_species_tree()
    species_stats_dict = get_species_stats_dict()
    with open(f"{summary_dir}/species_tree.newick", "w") as f:
        species_subtree = species_tree.copy()
        species_subtree.prune(gene_summary["origin_species"].unique().tolist())
        f.write(species_subtree.write())
    species = gene_summary["origin_species"].unique().tolist()
    species_subset = {k: v for k, v in species_stats_dict.items() if k in species}
    species_info_df = pd.DataFrame.from_dict(species_subset, orient="index")
    num_sampled_genes = gene_summary.groupby("origin_species").size()
    species_info_df["num_sampled_genes"] = num_sampled_genes
    species_info_df.to_csv(f"{summary_dir}/species_stats.csv")

    number_of_genes = gene_summary.shape[0]
    simulated_orthogroups = gene_summary["orthogroup"].unique().shape[0]
    simulated_up_regulated_genes = sum(gene_summary["simulation_fold_change"] >= 2.0)
    simulated_down_regulated_genes = sum(gene_summary["simulation_fold_change"] <= 0.5)
    up_regulated_genes = sum(gene_summary["actual_log2fc"] >= 1.0)
    down_regulated_genes = sum(gene_summary["actual_log2fc"] <= -1.0)

    gene_summary = move_column(gene_summary, "actual_log2fc", 6)

    gene_summary = gene_summary.rename(columns={"simulation_fold_change": "Fold Change used for simulating counts, differs from the actual fold change",
                                                "actual_log2fc": "Actual LOG2 Fold Change (includes +-inf for zero counts in one group)"})

    gene_summary.to_csv(f"{summary_dir}/gene_summary.csv", index=False)

    with open(f"{summary_dir}/simulation_stats.txt", "w") as f:
        f.write(f"Number of simulated genes: {number_of_genes}\n")
        f.write(f"Number of simulated orthogroups:  {simulated_orthogroups}\n")
        f.write(f"Dropped genes due to all zero assignment by distribution: {number_of_dropped_genes}\n")
        f.write(f"Dropped orthogroups due to all zero assignment by distribution: {specified_orthogroups - simulated_orthogroups}\n")
        f.write(f"Number of simulated up regulated genes: {simulated_up_regulated_genes} (percentage = {simulated_up_regulated_genes / number_of_genes})\n")
        f.write(f"Number of simulated down regulated genes: {simulated_down_regulated_genes} (percentage = {simulated_down_regulated_genes / number_of_genes})\n")
        f.write(f"Number of up regulated genes (accoring to actual fold change): {up_regulated_genes} (percentage = {up_regulated_genes / number_of_genes})\n")
        f.write(f"Number of down regulated genes (according to actual fold change): {down_regulated_genes} (percentage = {down_regulated_genes / number_of_genes})\n")


def move_column(df, col_name, new_pos):
    cols = list(df.columns)
    cols.insert(new_pos, cols.pop(cols.index(col_name)))
    return df[cols]


def create_sample_values(gene_summary_df, number_of_samples, first_group, a0, a1):
    """
    Generates a sparse matrix of sample values based on DESeq2 dispersion assumptions.

    Parameters:
        gene_summarary_df (pandas.DataFrame): The summary of genes.
        number_of_samples (int): The number of samples.

    Returns:
        pandas.DataFrame: The summary df including the count matrix for the samples.
    """
    if first_group:
        group = "group_1"
    else:
        group = "group_2"
        gene_summary_df = gene_summary_df.copy()
        gene_summary_df["read_mean_count"] = gene_summary_df["read_mean_count"] * gene_summary_df["simulation_fold_change"]

    dispersion_df = pd.DataFrame({
        "gene_name": gene_summary_df["gene_name"],
        f"estimated_dispersion_{group}" : [(a0 / mu) + a1 for mu in list(gene_summary_df["read_mean_count"])]
    })

    means = list(gene_summary_df["read_mean_count"])
    dispersions = 1 / dispersion_df[f"estimated_dispersion_{group}"].values

    with pm.Model() as _:
        _ = pm.NegativeBinomial(f"{group}_counts", mu=means, alpha=dispersions, shape=len(means))
        prior_predictive = pm.sample_prior_predictive(draws=number_of_samples)

    simulated_counts = prior_predictive.prior[f"{group}_counts"].values[0]

    sample_columns = [f"{group}_sample_{i + 1}" for i in range(number_of_samples)]
    simulated_data_matrix = pd.DataFrame(simulated_counts.T, columns=sample_columns)
    simulated_data_matrix.insert(0, "gene_name", dispersion_df["gene_name"])

    sample_disp_df = pd.merge(simulated_data_matrix, dispersion_df, on="gene_name")
    return sample_disp_df


def create_fastq_file(sample_df, sample_name, output_dir, gzip, mode, model, seed, read_length, threads):
    """
    Creates a fastq file for the sample using the InSilicoSeq (iss) package.

    Parameters:
        sample_df (pandas.DataFrame): The dataframe with the gene names and according absolute counts.
        sample_name (str): The name of the sample.
        output_dir (str): The output directory for the fastq files.
        gzip (bool): Whether the fastq files should be gzipped.
        mode (str): The mode for the simulation. Can be 'kde', 'perfect' or 'basic'.
        model (str): The error model for the reads Illumina (InSilicoSeq Simulation). Can be None.
        seed (int): The seed for the simulation. Can be None.
        read_length (int): The read length. Will only be used if the model is 'basic' or 'perfect'.
        threads (int): The number of threads to use.
    """
    read_count_file = f"{output_dir}/{sample_name}.tsv"
    number_of_pairs = sample_df[["gene_name", "absolute_numbers"]].copy()
    number_of_pairs["absolute_numbers"] = number_of_pairs["absolute_numbers"] * 2
    number_of_pairs.to_csv(read_count_file, sep="\t", index=False, header=False)

    args = argparse.Namespace(
        mode=mode,
        seed=seed,
        model=model,
        fragment_length=None,     # 300
        fragment_length_sd=None,  # 20
        read_length=read_length,
        store_mutations=True,
        genomes=[f"{output_dir}/summary/metatranscriptome_reference.fasta"],
        draft=None,
        ncbi=False,
        n_genomes_ncbi=0,
        output=f"{output_dir}/{sample_name}",
        n_genomes=None,
        readcount_file=read_count_file,
        abundance_file=None,
        coverage_file=None,
        coverage=None,
        abundance=None,
        n_reads=None,
        cpus=threads,
        sequence_type="metagenomics",
        gc_bias=False,
        compress=gzip,
        debug=True,
        quiet=False
    )
    generate_reads(args)
    if os.path.exists(read_count_file):
        os.remove(read_count_file)
    else:
        print(f"Warning: Could not remove {read_count_file}.")


def draw_library_sizes(library_size, library_size_distribution, number_of_samples):
    """
    Draws the library sizes for each sample according to the specified distribution.

    Parameters:
        library_size (int): The base library size.
        library_size_distribution (LibrarySizeDistribution): The distribution of the library sizes.
        number_of_samples (int): The number of samples, this is how many library sizes will be drawn.

    Returns:
        list: The library sizes for each sample.
    """
    match library_size_distribution.distribution_name:
        case LibrarySizeDistribution.poisson:
            poisson_lambda = library_size_distribution.poisson
            sample_library_sizes = library_size * (np.random.poisson(poisson_lambda, number_of_samples) / poisson_lambda)
        case LibrarySizeDistribution.uniform:
            sample_library_sizes = [library_size] * number_of_samples
        case LibrarySizeDistribution.negative_binomial:
            nbin_n, nbin_p = library_size_distribution.nbin_n, library_size_distribution.nbin_p
            expected_mean = nbin_n * (1 - nbin_p) / nbin_p
            sample_library_sizes = (library_size * (np.random.negative_binomial(nbin_n, nbin_p, number_of_samples) / expected_mean)).round()
    return sample_library_sizes


def scale_fastq_samples(gene_summary_df, sample_library_sizes):
    """
    Scales each sample in the gene_summary_df to scaled library size.

    Parameters:
        gene_summary_df (pandas.DataFrame): Dataframe containing simuzlation information
        sample_library_sizes (list): The library sizes for each sample.
        read_length (int): The read length. Will only be used if the model is 'basic' or 'perfect'.
        threads (int): The number of threads to use.
    Returns:
        pandas.DataFrame: The filtered dataframe with all-zero columns removed.
        """
    # scale to library size
    gene_summary_df["gene_mean_scaled_to_library_size"] = (gene_summary_df["read_mean_count"] / gene_summary_df["read_mean_count"].sum()) * sample_library_sizes[0]
    # multiply each sample with scaled library size, scaling and ceiling results to avoid float values
    for sample, sample_library_size in zip([col for col in gene_summary_df.columns if "sample" in col], sample_library_sizes):
        gene_summary_df[sample] = (gene_summary_df[sample] / gene_summary_df[sample].sum()) * sample_library_size
        gene_summary_df[sample] = np.ceil(gene_summary_df[sample]).astype(int)

    return gene_summary_df


def determine_mode_and_model(model, read_length):
    mode = "kde"
    if model == ErrorModel.basic or model == ErrorModel.perfect:
        mode = model
        model = None
    elif read_length:
        print("Warning: Read length is ignored if model is not 'basic' or 'perfect'.")
        # TODO write the read length of the selected model
    return mode, model


def create_fastq_samples(gene_summary_df, outdir, compression, mode, model, seed, read_length, threads, bar):
    """
    Calls the create_fastq_file function for each sample in the gene_summary_df.
    Parameters:
        outdir (str): The output directory for the fastq files.
        compression (bool): Whether the fastq files should be gzipped.
        model (ErrorModel): The error model for the reads (Illumina).
        seed (int): The seed for the simulation. Can be None.
    """
    for sample in [col for col in gene_summary_df.columns if "sample" in col]:
        sample_copy = gene_summary_df[["gene_name", sample]].copy()
        sample_copy.rename(columns={sample: "absolute_numbers"}, inplace=True)
        create_fastq_file(sample_copy, sample, outdir, compression, mode, model, seed, read_length, threads)
        bar.next()
    bar.finish()


def get_all_zero_genes(gene_summary_df):
    """
    Get all genes of the gene_summary_df with all-zero columns.
    Parameters:
        gene_summary_df (pandas.DataFrame): The gene summary dataframe.
    Returns:
        set: A set of all zero genes.
    """
    pl_df = pl.DataFrame(gene_summary_df)
    sample_cols = [col for col in pl_df.columns if "sample" in col]
    all_zero_genes = (
        pl_df.filter(pl.all_horizontal([pl.col(col) == 0 for col in sample_cols]))
        .select("gene_name")
        .to_series()
        .to_list()
    )
    return set(all_zero_genes)


def draw_dge_factors(dge_ratio, number_of_selected_genes):
    """"
    Draws the log2 DGE factors from a normal distribution adjusted to the specified ratio of up and downregulated genes.
    We use the normal and log2 because it is easier to calculate and then transform with exp2 to get the actual fold changes

    Parameters:
        dge_ratio (float): The ratio of up and downregulated genes
        number_of_selected_genes (int): The number of selected genes

    Returns:
        numpy.ndarray: The differentialy expressed factors
    """
    dge_ratio = dge_ratio / 2

    if number_of_selected_genes == 0:
        return np.array([])

    if not (0 < dge_ratio < 0.5):
        raise ValueError("dge_ratio must be between 0 and 0.5 (exclusive)")

    z_score = stats.norm.ppf(1 - dge_ratio)
    sigma = DGE_LOG_2_CUTOFF_VALUE / z_score

    with pm.Model() as _:
        _ = pm.Normal("dge_ratios", mu=0, sigma=sigma)
        prior_predictive = pm.sample_prior_predictive(draws=number_of_selected_genes)

    simulated_ratios = prior_predictive.prior['dge_ratios'].values[0]
    simulated_ratios = np.exp2(simulated_ratios)

    return simulated_ratios


def calculate_species_identity(species_id, chosen_species, id_sets):
    combined_species = chosen_species + [species_id]
    filtered_id_sets = id_sets[id_sets.apply(lambda x: any(s_id in x for s_id in combined_species))]
    average_group_identity = filtered_id_sets.apply(lambda x: sum(s_id in x for s_id in combined_species)).mean()
    return species_id, average_group_identity


def maximize(x, y):
    return x > y


def minimize(x, y):
    return x < y


def select_species_with_criterion(number_of_species, number_of_threads, selection_criterion: SelectionCriterion = SelectionCriterion.maximize):
    if selection_criterion not in SelectionCriterion:
        raise ValueError(f"Invalid selection criterion: {selection_criterion}. Must be one of {list(SelectionCriterion)}")
    match selection_criterion:
        case SelectionCriterion.maximize:
            comparator = maximize
            initial_best_value = 0
        case SelectionCriterion.minimize:
            comparator = minimize
            initial_best_value = 10000

    # TODO: question: should I make it random or should I start with the highest pair? -> ask Stefan
    random_species = random.randint(0, MAX_SPECIES)
    chosen_species = [random_species]
    # this takes about 2 minutes 38sec, so it could be precomputed, would increase load on LFS and increase precompution steps
    pg_overview = get_pg_overview()
    id_sets = pg_overview.apply(lambda x: frozenset([i for i in range(0, MAX_SPECIES) if x.iloc[i] != "-"]), axis=1)

    for _ in range(number_of_species - 1):
        species_left = [i for i in range(0, MAX_SPECIES) if i not in chosen_species]
        # i hope id_sets is thread safe, as i only view it
        results = Parallel(n_jobs=number_of_threads)(
            delayed(calculate_species_identity)(species_id, chosen_species, id_sets)
            for species_id in species_left
        )
        best_value = initial_best_value
        best_species = -1
        for species_id, average_group_identity in results:
            if comparator(average_group_identity, best_value):
                best_value = average_group_identity
                best_species = species_id
        chosen_species.append(best_species)

    index_species_dict = dict(zip(range(0, MAX_SPECIES), pg_overview.columns[:MAX_SPECIES]))
    return [index_species_dict[species] for species in chosen_species]


def select_orthogroups(orthogroup_slice, species, number_of_groups, minimize=True, force=False):
    orthogroups = orthogroup_slice[species].copy()
    number_of_species = len(species)
    orthogroups["group_size"] = orthogroups.apply(lambda x: number_of_species - (x == "-").sum(), axis=1)
    orthogroups = orthogroups[orthogroups["group_size"] > 0]
    orthogroups = orthogroups.sample(frac=1).reset_index(drop=True)
    orthogroups = orthogroups.sort_values(by="group_size", ascending=minimize)
    check_force_result = check_enough_orthogroups(orthogroups, number_of_groups, force)
    if check_force_result is not None:
        return check_force_result

    return orthogroups.head(number_of_groups)


def calc_zero_ratio(df):
    df = df[~(df == 0).all(axis=1)]
    return (df == 0).sum().sum() / df.size


def add_extra_sparsity(gene_summary_df, sparsity_target, seed):
    """"
    Adds more zeros to dataframe to reach sparsity target. Will avoid all-zero rows. This means target sparsity may not be reached.

    Parameters:
        gene_summary_df (pd.Dataframe): The ratio of up and downregulated genes
        sparsity_target (float): Sparsity target, between 0 and 1
        seed (int): Random seed for reproducibility.

    Returns:
        pd.Dataframe: The dataframe with added zeros, provided sparsity target is not already reached.
    """
    gene_summary_df = pl.DataFrame(gene_summary_df)
    sample_cols = [col for col in gene_summary_df.columns if "sample" in col]
    df_filtered = gene_summary_df.select(sample_cols)
    df_np = df_filtered.to_numpy().copy()
    total_cells = df_np.size
    current_zeros = np.count_nonzero(df_np == 0)
    target_zeros = int(np.ceil(total_cells * sparsity_target))

    n_to_zero = target_zeros - current_zeros

    if n_to_zero <= 0:
        return gene_summary_df.to_pandas()

    rng = np.random.default_rng(seed)
    safe_indices = []

    for i, row in enumerate(df_np):
        nonzero_cols = np.flatnonzero(row != 0)
        if len(nonzero_cols) > 1:
            reserved = rng.choice(nonzero_cols, 1)
            available = [c for c in nonzero_cols if c != reserved]
            safe_indices.extend((i, c) for c in available)

    if not safe_indices:
        return gene_summary_df.to_pandas()

    safe_indices = np.array(safe_indices)

    n_to_zero = min(n_to_zero, len(safe_indices))

    if n_to_zero < target_zeros - current_zeros:
        print(f"Can only zero {n_to_zero} cells without causing all-zero rows.")

    selected = rng.choice(len(safe_indices), n_to_zero, replace=False)
    rows, cols = safe_indices[selected].T
    df_np[rows, cols] = 0

    updated_sample_cols = [
        pl.Series(name, df_np[:, i]) for i, name in enumerate(sample_cols)
    ]
    return gene_summary_df.with_columns(updated_sample_cols).to_pandas()


def add_actual_log2fc(gene_summary_df):
    gene_summary_df = pl.DataFrame(gene_summary_df)
    sample_cols_group_1 = [col for col in gene_summary_df.columns if "group_1" in col]
    sample_cols_group_2 = [col for col in gene_summary_df.columns if "group_2" in col]

    gene_summary_df = gene_summary_df.with_columns([
        pl.sum_horizontal(sample_cols_group_1).alias("sum_group1"),
        pl.sum_horizontal(sample_cols_group_2).alias("sum_group2"),
    ]).with_columns([
        (pl.col("sum_group2") / (pl.col("sum_group1"))).alias("fold_change")
    ]).with_columns([
        pl.col("fold_change").log(2).alias("actual_log2fc"),
    ]).drop(["sum_group1", "sum_group2", "fold_change"])

    return gene_summary_df.to_pandas()


def add_counts_to_large_orthogroups(gene_summary, species_count):
    og_sizes = gene_summary["orthogroup"].value_counts()
    large_orthogroups = og_sizes[og_sizes == species_count].index.to_list()
    og_to_modify = random.choice(large_orthogroups)

    filtered_gene_summary = gene_summary[gene_summary["orthogroup"] == og_to_modify]

    count_cols = [col for col in gene_summary.columns if "sample" in col]
    row_indices = []

    for row in filtered_gene_summary[count_cols].iterrows():
        if sum(row[1]) == 0:
            row_indices.append(row[0])

    for row_index in row_indices:
        sample_to_one = random.choice(count_cols)
        gene_summary.at[row_index, sample_to_one] = 1
