import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, IntVector
from Bio import SeqIO

from presets import *


def draw_random_species(number_of_species):
    return random.sample(AVAILABLE_SPECIES, number_of_species)


def create_ortholgous_group_rates(number_of_orthogous_groups, max_species_per_group):
    with model:
        orthologues_samples = pm.sample_prior_predictive(number_of_orthogous_groups, var_names=['ortho'])
    orthogroups = orthologues_samples.to_dataframe()["ortho"].to_numpy()
    scaled_orthogroups = orthogroups * (max_species_per_group / np.max(orthogroups))
    scaled_orthogroups = np.ceil(scaled_orthogroups)
    return scaled_orthogroups


#randomization based on rates calculated from the pdf
def draw_orthogroups_by_rate(orthogroup_scales, species):
    orthogroup_scales = pd.Series(orthogroup_scales)
    value_counts = orthogroup_scales.value_counts()
    orthogroups = pg_overview[species]
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


def generate_species_abundance(number_of_species):
    with model:
        species_samples = pm.sample_prior_predictive(number_of_species, var_names=['species'])
    return species_samples.to_dataframe()['species'].to_list()


def generate_read_mean_counts(number_of_reads):
    with model:
        reads = pm.sample_prior_predictive(number_of_reads, var_names=['reads'])
    return reads.to_dataframe()["reads"].to_list()


def generate_reads(base_expression_values, replicates_per_sample, filtered_genes_file):
    polyester = importr('polyester')
    reads_per_transcript = FloatVector(base_expression_values)
    num_reps = IntVector(replicates_per_sample)
    number_of_transcripts = len(base_expression_values)
    fold_changes = [[1,1] for _ in range(number_of_transcripts)]
    fold_changes_r = robjects.r.matrix(FloatVector([elem for sublist in fold_changes for elem in sublist]), nrow=number_of_transcripts, byrow=True)
    polyester.simulate_experiment(
        filtered_genes_file,
        reads_per_transcript=reads_per_transcript,
        num_reps=num_reps,
        fold_changes=fold_changes_r,
        outdir='simulated_reads'
    )


def filter_genes_from_ground(gene_list, output_fasta):
    """
    Filters sequences from a FASTA file based on a list of gene names and writes them to a new file.
    
    Parameters:
    gene_list (list of str): List of gene names to be filtered.
    input_fasta (str): Path to the input FASTA file.
    output_fasta (str): Path to the output FASTA file where filtered sequences will be saved.
    """
    with open(PATH_TO_GROUND_GENES, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if any(gene in record.id for gene in gene_list):
                SeqIO.write(record, outfile, "fasta")