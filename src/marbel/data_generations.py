import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, IntVector
from Bio import SeqIO, bgzf
from pathlib import Path
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

from marbel.presets import AVAILABLE_SPECIES, model, pm, pg_overview, species_tree, PATH_TO_GROUND_GENES_INDEX, DEFAULT_PHRED_QUALITY


def draw_random_species(number_of_species):
    """
    Draws a random sample of species from the list of available species set in the presets.

    Parameters:
    number_of_species (int): The number of species to draw.

    Returns:
    list: A list of randomly drawn species from the list of available species.
    """
    return random.sample(AVAILABLE_SPECIES, number_of_species)


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
    with model:
        orthologues_samples = pm.sample_prior_predictive(number_of_orthogous_groups, var_names=['ortho'], random_seed=seed)
    orthogroups = orthologues_samples.to_dataframe()["ortho"].to_numpy()
    scaled_orthogroups = orthogroups * (max_species_per_group / np.max(orthogroups))
    scaled_orthogroups = np.ceil(scaled_orthogroups)
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
    if max_phylo_distance is not None:
        filtered_slice = pg_overview[pg_overview["tip2tip_distance"] <= max_phylo_distance]
    if min_identity is not None:
        filtered_slice = filtered_slice[filtered_slice["medium_identity"] >= min_identity]
    if max_phylo_distance is None and min_identity is None:
        return pg_overview
    return filtered_slice


#  randomization based on rates calculated from the pdf
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


def draw_orthogroups(orthogroup_slice, number_of_orthogous_groups, species):
    """
    Draws orthologous groups based on actual occurences in the dataset, instead of rates based on the pdf.
    Given a dataframe slice of orthologous groups, a number of orthologous groups to be drawn, and a list of species,
    randomly samples orthologous groups based on their actual occurences. If there are not enough orthologous groups
    to satisfy the parameters, it raises an error.

    Parameters:
        orthogroup_slice (pandas.DataFrame): The orthologous group dataframe.
        number_of_orthogous_groups (int): The number of orthologous groups to be drawn.
        species (list): A list of species.

    Returns:
        pandas.DataFrame: A randomly sampled dataframe of orthologous groups based on their actual occurences and filtered by species.
    """
    orthogroups = orthogroup_slice[species].copy()
    orthogroups["group_size"] = orthogroups.apply(lambda x: len(species) - len(x[x == "-"].index.to_list()), axis=1)
    orthogroups = orthogroups[orthogroups["group_size"] > 0]
    if orthogroups.shape[0] < number_of_orthogous_groups:
        print("Error: Not enough orthogroups to satisfy the parameters, specify different parameters, i.e. lower orthogroups and less stringent sequence similarity and allow more phygenetic distance.")
        quit()
    orthogroups_sample = orthogroups.sample(n=number_of_orthogous_groups)
    return orthogroups_sample


def generate_species_abundance(number_of_species, seed=None):
    """
    Generates a list of species abundances based on the distribution in the preset model.

    Parameters:
    number_of_species (int): The number of species for which to generate abundances.
    seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
    list: A list of species abundances.
    """
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
    with model:
        reads = pm.sample_prior_predictive(number_of_reads, var_names=['reads'], random_seed=seed)
    return reads.to_dataframe()["reads"].to_list()


def generate_fold_changes(number_of_transcripts, dge_ratio):
    """
    Generates a list of fold changes for each transcript based on the given ratio of up and down regulated genes.

    Parameters:
    number_of_transcripts (int): The total number of transcripts.
    dge_ratio (tuple): A tuple representing the ratio of up regulated genes and down regulated genes. The first value
                        represents the ratio of up regulated genes, the second represents the ratio of down regulated
                        genes.

    Returns:
    list: A list of fold changes for each transcript. Each fold change is represented as a list of two values
          representing the fold change for up and down regulation respectively.
    """
    dge_genes = random.sample(range(number_of_transcripts), int(sum(dge_ratio) * number_of_transcripts))
    up_regulated = dge_genes[:int((dge_ratio[0] / sum(dge_ratio)) * len(dge_genes))]
    fold_changes = [[1, 1] if i not in dge_genes else [2, 1] if i in up_regulated else [1, 2] for i in range(number_of_transcripts)]
    return fold_changes


def generate_reads(gene_summarary_df, replicates_per_sample, filtered_genes_file, outdir, dge_ratio, seed, read_length):
    """
    Generates reads for a given dataset using the polyester package.

    Parameters:
    gene_summarary_df (pandas.DataFrame): A DataFrame containing information about the genes.
    replicates_per_sample (list): A list of integers representing the number of replicates per sample. Should contain two values.
    filtered_genes_file (str): The path to the filtered genes file. Based on this file, the reads will be generated.
    outdir (str): The output directory for the generated reads.
    dge_ratio (tuple): A tuple representing the ratio of up regulated genes and down regulated genes. The first value
                        represents the ratio of up regulated genes, the second represents the ratio of down regulated
                        genes.
    seed (int, optional): Random seed for reproducibility. Defaults to None.
    read_length (int): The length of the reads.
    """
    polyester = importr('polyester')
    base_expression_values = gene_summarary_df["read_mean_count"].to_list()

    reads_per_transcript = FloatVector(base_expression_values)
    num_reps = IntVector(replicates_per_sample)
    number_of_transcripts = len(base_expression_values)
    fold_changes = generate_fold_changes(number_of_transcripts, dge_ratio)
    fold_changes_r = robjects.r.matrix(FloatVector([elem for sublist in fold_changes for elem in sublist]), nrow=number_of_transcripts, byrow=True)
    gene_summarary_df["fold_change_ratio"] = [float(i[0]) / float(i[1]) for i in fold_changes]

    if seed:
        polyester.simulate_experiment(
            filtered_genes_file,
            reads_per_transcript=reads_per_transcript,
            num_reps=num_reps,
            fold_changes=fold_changes_r,
            outdir=outdir,
            seed=seed,
            readlen=read_length
        )
    else:
        polyester.simulate_experiment(
            filtered_genes_file,
            reads_per_transcript=reads_per_transcript,
            num_reps=num_reps,
            fold_changes=fold_changes_r,
            outdir=outdir
        )


def filter_genes_from_ground(gene_list, output_fasta):
    """
    Writes genes from a list to a FASTA file. Duplicates will also be written, the order is perserved.

    Parameters:
    gene_list (list): A list of gene names.
    output_fasta (str): Path to the output FASTA file where the sequences will be saved.
    """
    ground_genes = SeqIO.index_db(PATH_TO_GROUND_GENES_INDEX)
    with open(output_fasta, "w") as outfile:
        for gene in gene_list:
            if gene in ground_genes:
                SeqIO.write(ground_genes[gene], outfile, "fasta")
            else:
                print(f"Critical Warning: Gene {gene} not found in the ground genes.")


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


def convert_fasta_dir_to_fastq_dir(fasta_dir, gzipped=True):
    """
    Converts a directory containing .fasta files to a directory containing .fastq files. If gzipped is True, the output files will be gzipped.

    Parameters:
    - fasta_dir (str): The path to the directory containing the .fasta files.
    - gzipped (bool): Whether the output files should be gzipped.

    Note that the input .fasta files will be removed after the conversion is done.
    """
    fasta_dir = Path(fasta_dir)
    futures = []
    max_workers = multiprocessing.cpu_count()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for fa_path in fasta_dir.glob("*.fasta"):
            fq_path = fasta_dir / fa_path.with_suffix(".fastq.gz" if gzipped else ".fastq").name
            if gzipped:
                futures.append(executor.submit(write_as_fastq_gz, fa_path, fq_path))
            else:
                futures.append(executor.submit(write_as_fastq, fa_path, fq_path))
        for future in as_completed(futures):
            future.result()
    for fa_path in fasta_dir.glob("*.fasta"):
        os.remove(fa_path)


def write_as_fastq_gz(fa_path, fq_path):
    """
    Converts a .fasta file to a .fastq.gz file. The function reads the .fasta file,
    adds phred quality scores to each sequence and writes the output to a .fastq.gz file.

    Args:
        fa_path (str): Path to the input .fasta file.
        fq_path (str): Path to the output .fastq.gz file.
    """
    with open(fa_path, "r") as fasta, bgzf.BgzfWriter(fq_path, "wb") as fastq_gz:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [DEFAULT_PHRED_QUALITY] * len(record)
            SeqIO.write(record, fastq_gz, "fastq")


def write_as_fastq(fa_path, fq_path):
    """
    Converts a .fasta file to a .fastq file. The function reads the .fasta file,
    adds phred quality scores to each sequence and writes the output to a .fastq file.

    Args:
        fa_path (str): Path to the input .fasta file.
        fq_path (str): Path to the output .fastq file.
    """
    with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [DEFAULT_PHRED_QUALITY] * len(record)
            SeqIO.write(record, fastq, "fastq")


def summarize_parameters(number_of_orthogous_groups, number_of_species, number_of_sample,
                         outdir, max_phylo_distance, min_identity, deg_ratio, seed, output_format, read_length, result_file):
    """
    Writes the simulation parameters to the result_file.

    Args:
        number_of_orthogous_groups (int): The number of orthologous groups.
        number_of_species (int): The number of species.
        number_of_sample (tuple): The number of samples (group 1, group 2).
        outdir (str): The output directory.
        max_phylo_distance (float): The maximum phylogenetic distance.
        min_identity (float): The minimum sequence identity.
        deg_ratio (tuple): The ratio of up and down regulated genes (up, down).
        seed (int): The seed for the simulation.
        output_format (str): The output format.
        read_length (int): The read length.
        result_file (file): The file to write the summary to.
    """
    result_file.write(f"Number of orthogroups: {number_of_orthogous_groups}\n")
    result_file.write(f"Number of species: {number_of_species}\n")
    result_file.write(f"Number of samples: {number_of_sample}\n")
    result_file.write(f"Output directory: {outdir}\n")
    result_file.write(f"Max phylogenetic distance: {max_phylo_distance}\n")
    result_file.write(f"Min identity: {min_identity}\n")
    result_file.write(f"Up and down regulated genes: {deg_ratio}\n")
    result_file.write(f"Seed: {seed}\n")
    result_file.write(f"Output format: {output_format}\n")
    result_file.write(f"Read length: {read_length}\n")


def generate_report(number_of_orthogous_groups, number_of_species, number_of_sample,
                    outdir, max_phylo_distance, min_identity, deg_ratio, seed, output_format, gene_summary, read_length):
    """
    Generates a report of the simulation results.

    Parameters:
        number_of_orthogous_groups (int): The number of orthologous groups.
        number_of_species (int): The number of species.
        number_of_sample (tuple): The number of samples (group 1, group 2).
        outdir (str): The output directory.
        max_phylo_distance (float): The maximum phylogenetic distance.
        min_identity (float): The minimum sequence identity.
        deg_ratio (tuple): The ratio of up and down regulated genes (up, down).
        seed (int): The seed for the simulation.
        output_format (str): The output format.
        gene_summary (pandas.DataFrame): The summary of genes.
        read_length (int): The read length.
    """
    summary_dir = f"{outdir}/summary"
    with open(f"{summary_dir}/marbel_params.txt", "w") as f:
        summarize_parameters(number_of_orthogous_groups, number_of_species, number_of_sample, outdir,
                             max_phylo_distance, min_identity, deg_ratio, seed, output_format, read_length, f)
    gene_summary.to_csv(f"{summary_dir}/gene_summary.csv", index=False)
    with open(f"{summary_dir}/species_tree.newick", "w") as f:
        species_subtree = species_tree.copy()
        species_subtree.prune(gene_summary["origin_species"].unique().tolist())
        f.write(species_subtree.write())
