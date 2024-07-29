import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, IntVector
from Bio import SeqIO, bgzf
from pathlib import Path
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

from presets import *


def draw_random_species(number_of_species):
    return random.sample(AVAILABLE_SPECIES, number_of_species)


def create_ortholgous_group_rates(number_of_orthogous_groups, max_species_per_group, seed=None):
    with model:
        orthologues_samples = pm.sample_prior_predictive(number_of_orthogous_groups, var_names=['ortho'], random_seed=seed)
    orthogroups = orthologues_samples.to_dataframe()["ortho"].to_numpy()
    scaled_orthogroups = orthogroups * (max_species_per_group / np.max(orthogroups))
    scaled_orthogroups = np.ceil(scaled_orthogroups)
    return scaled_orthogroups


def filter_by_seq_id_and_phylo_dist(max_phylo_distance, min_identity):
    if not max_phylo_distance is None:
        filtered_slice = pg_overview[pg_overview["tip2tip_distance"] <= max_phylo_distance]
    if not min_identity is None:
        filtered_slice = filtered_slice[filtered_slice["medium_identity"] >= min_identity]
    if max_phylo_distance is None and min_identity is None:
        return pg_overview
    return filtered_slice


#randomization based on rates calculated from the pdf
def draw_orthogroups_by_rate(orthogroup_slice, orthogroup_scales, species):
    orthogroup_scales = pd.Series(orthogroup_scales)
    value_counts = orthogroup_scales.value_counts()
    orthogroups = orthogroup_slice[species].copy()
    number_of_species = len(species)
    orthogroups["group_size"] = orthogroups.apply(lambda x: number_of_species - len(x[x=="-"]), axis=1)
    orthogroups = orthogroups[orthogroups["group_size"] > 0]
    sampled_groups = pd.DataFrame(columns=orthogroups.columns)
    for index in value_counts.index:
        subset_orthogroups = orthogroups[orthogroups["group_size"] == index]
        subset_orthogroups = subset_orthogroups.sample(value_counts[index])
        sampled_groups = pd.concat([sampled_groups, subset_orthogroups])
    return sampled_groups


#randomization based on actual occurences, instead based on the pdf
def draw_orthogroups(number_of_orthogous_groups, species):
    pg_overview_filter = pg_overview[species]
    max_species_per_group = len(species)
    orthogroups = pg_overview_filter.sample(n=number_of_orthogous_groups)
    species_subset = draw_random_species(max_species_per_group)
    orthogroups = orthogroups[species_subset]
    orthogroups = orthogroups[orthogroups["group_size"] > 0]
    return orthogroups.sample(n=number_of_orthogous_groups)


def generate_species_abundance(number_of_species, seed=None):
    with model:
        species_samples = pm.sample_prior_predictive(number_of_species, var_names=['species'], random_seed=seed)
    return species_samples.to_dataframe()['species'].to_list()


def generate_read_mean_counts(number_of_reads, seed=None):
    with model:
        reads = pm.sample_prior_predictive(number_of_reads, var_names=['reads'], random_seed=seed)
    return reads.to_dataframe()["reads"].to_list()


def generate_reads(base_expression_values, replicates_per_sample, filtered_genes_file, outdir, dge_ratio, seed):
    polyester = importr('polyester')
    reads_per_transcript = FloatVector(base_expression_values)
    num_reps = IntVector(replicates_per_sample)
    number_of_transcripts = len(base_expression_values)
    dge_genes = random.sample(range(number_of_transcripts), int(sum(dge_ratio) * number_of_transcripts))
    up_regulated = dge_genes[:int((dge_ratio[0]/sum(dge_ratio)) * len(dge_genes))]
    fold_changes = [[1,1] if i not in dge_genes else [2, 1] if i in up_regulated else [1, 2] for i in range(number_of_transcripts)]
    fold_changes_r = robjects.r.matrix(FloatVector([elem for sublist in fold_changes for elem in sublist]), nrow=number_of_transcripts, byrow=True)
    if seed:
        polyester.simulate_experiment(
            filtered_genes_file,
            reads_per_transcript=reads_per_transcript,
            num_reps=num_reps,
            fold_changes=fold_changes_r,
            outdir=outdir,
            seed=seed
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
    gene_list (list of str): List of gene names.
    output_fasta (str): Path to the output FASTA file where the sequences will be saved.
    """
    ground_genes = SeqIO.index_db(PATH_TO_GROUND_GENES_INDEX)
    with open(output_fasta, "w") as outfile:
        for gene in gene_list:
            if gene in ground_genes:
                SeqIO.write(ground_genes[gene], outfile, "fasta")
            else:
                print(f"Critical Warning: Gene {gene} not found in the ground genes.")


def extract_combined_gene_names_and_weigths(species, species_abundances, selected_ortho_groups, read_mean_counts):
    scaled_read_mean_counts, all_species_genes = [], []
    current_read_index, i = 0, 0
    species_weights = species_abundances / np.sum(species_abundances)

    for sp in species:
        species_genes_list = selected_ortho_groups[selected_ortho_groups[sp] != "-"][sp].to_list()
        scaled_read_mean_counts += [species_weights[i] * c for c in read_mean_counts[current_read_index:(current_read_index + len(species_genes_list))]]
        current_read_index += len(species_genes_list)
        i += 1
        all_species_genes += species_genes_list
    return scaled_read_mean_counts, all_species_genes


def convert_fasta_dir_to_fastq_dir(fasta_dir, gzipped=True):
    fasta_dir = Path(fasta_dir)
    futures = []
    max_workers = multiprocessing.cpu_count()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for fa_path in fasta_dir.glob("*.fasta"):
            fq_path = fasta_dir / fa_path.with_suffix(".fq.gz" if gzipped else ".fq").name
            if gzipped:
                futures.append(executor.submit(write_as_fastq_gz, fa_path, fq_path))
            else:
                futures.append(executor.submit(write_as_fastq, fa_path, fq_path))
        for future in as_completed(futures):
            future.result()


def write_as_fastq_gz(fa_path, fq_path):
    with open(fa_path, "r") as fasta, bgzf.BgzfWriter(fq_path, "wb") as fastq_gz:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [DEFAULT_PHRED_QUALITY] * len(record)
            SeqIO.write(record, fastq_gz, "fastq")


def write_as_fastq(fa_path, fq_path):
    with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [DEFAULT_PHRED_QUALITY] * len(record)
            SeqIO.write(record, fastq, "fastq")
