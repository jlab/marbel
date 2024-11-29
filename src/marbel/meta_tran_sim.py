from typing import Optional, Tuple, List
from typing_extensions import Annotated

import typer
import random
import numpy as np
import os
import pandas as pd
import uuid
import re
from progress.bar import Bar

from marbel.presets import __version__, MAX_SPECIES, MAX_ORTHO_GROUPS, rank_distance, LibrarySizeDistribution, Rank, ErrorModel, DESEQ2_FITTED_A0, DESEQ2_FITTED_A1
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, create_sample_values, create_fastq_samples, draw_library_sizes
from marbel.data_generations import draw_orthogroups_by_rate, draw_orthogroups, generate_species_abundance, generate_read_mean_counts, aggregate_gene_data, filter_genes_from_ground, generate_report
from marbel.data_generations import draw_dge_factors, write_parameter_summary

app = typer.Typer()


def version_callback(value: bool):
    if value:
        print(f"marbel Version: {__version__}")
        raise typer.Exit()


def species_callback(value: Optional[int]):
    if value > MAX_SPECIES:
        raise typer.BadParameter(f"The current maximum number of species is {MAX_SPECIES}")
    if value < 1:
        raise typer.BadParameter("The number of species is 1")
    return value


def orthogroups_callback(value: Optional[int]):
    if value > MAX_ORTHO_GROUPS:
        raise typer.BadParameter(f"The current maximum number of orthologous groups is {MAX_ORTHO_GROUPS}")
    if value < 1:
        raise typer.BadParameter("The number of orthologous groups is 1")
    return value


def checknegative(value: float):
    if value < 0:
        raise typer.BadParameter("Deseq2 dispersion values cannot be negative.")
    return value


def sample_callback(value: Optional[Tuple[int, int]]):
    if value[0] < 1 or value[1] < 1:
        raise typer.BadParameter("The minimum number of samples has to be 1")
    return value


def dge_ratio_callback(value: float):
    if value < 0:
        raise typer.BadParameter("Ratio cannot be negative")
    if value >= 1:
        raise typer.BadParameter("DGE ratio must be smaller than 1")
    return value


def library_size_distribution_callback(value):
    value = value.split(",")
    if value[0] not in LibrarySizeDistribution.possible_distributions:
        raise typer.BadParameter(f"Library size distribution {value[0]} is not a valid distribution. Choose from {LibrarySizeDistribution.possible_distributions}")
    if value[0] == "negative_binomial":
        if len(value[1:]) > 2 or len(value[1:]) == 1:
            raise typer.BadParameter("Negative binomial distribution requires two parameters")
        if value[1:] == 2:
            try:
                value[1] = int(value[1])
                value[2] = float(value[2])
            except ValueError:
                raise typer.BadParameter(f"Negative binomial distribution requires n to be an integer and p to be a float. Given: n: {value[1]} and p: {value[2]}")
            if value[1] < 1:
                raise typer.BadParameter(f"Negative binomial distribution requires n to be greater than 1. Given: {value[1]}")
            if value[2] < 0 or value[2] > 1:
                raise typer.BadParameter(f"Negative binomial distribution requires p to be between 0 and 1. Given: {value[2]}")
            return LibrarySizeDistribution(value[0], nbin_n=value[1], nbin_p=value[2])
        return LibrarySizeDistribution(value[0])
    if value[0] == "poisson":
        try:
            value[1] = int(value[1])
        except ValueError:
            raise typer.BadParameter(f"Poisson distribution requires an integer parameter. Given: {value[1]}")
        if len(value[1:]) > 1:
            raise typer.BadParameter("Poisson distribution requires one parameter")
        if len(value[1:]) == 1:
            return LibrarySizeDistribution(value[0], poisson=value[1])
        return LibrarySizeDistribution(value[0])
    return LibrarySizeDistribution(value[0])


def rank_species_callback(value: Optional[str]):
    if value is None:
        return None
    try:
        return rank_distance[value]
    except KeyError:
        try:
            return float(value)
        except ValueError:
            raise typer.BadParameter(f"Rank {value} is not a valid rank or a valid float. Choose from {list(rank_distance.keys())} or specify a float.")


def bar_next(bar):
    bar.next()
    print()


@app.command()
def main(n_species: Annotated[int, typer.Option(callback=species_callback,
                                                help="Number of species to be drawn for the metatranscriptomic in silico dataset")] = 20,
         n_orthogroups: Annotated[int,
                                  typer.Option(callback=orthogroups_callback,
                                               help="Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset")] = 1000,
         n_samples: Annotated[Tuple[int, int],
                              typer.Option(callback=sample_callback,
                                           help="Number of samples to be created for the metatranscriptomic in silico dataset"
                                           + "the first number is the number of samples for group 1 and"
                                           + "the second number is the number of samples for group 2"
                                           )] = [10, 10],
         outdir: Annotated[str, typer.Option(help="Output directory for the metatranscriptomic in silico dataset")] = "simulated_reads",
         max_phylo_distance: Annotated[Rank, typer.Option(callback=rank_species_callback, help="Maximimum mean phylogenetic distance for orthologous groups."
                                                          + "specify stricter limit, if you want to avoid orthologous groups"
                                                          + "with a more diverse phylogenetic distance.")] = None,
         min_identity: Annotated[float, typer.Option(help="Minimum mean sequence identity score for an orthologous groups."
                                                          + "Specify for more ")] = None,
         dge_ratio: Annotated[float, typer.Option(callback=dge_ratio_callback, help="Ratio of up and down regulated genes. Must be between 0 and 1."
                                                  "This is a random drawing process from normal distribution, so the actual ratio might vary.")] = 0.1,
         seed: Annotated[int, typer.Option(help="Seed for the sampling. Set for reproducibility")] = None,
         error_model: Annotated[ErrorModel, typer.Option(help="Sequencer model for the reads, use basic or perfect (no errors) for custom read length")] = ErrorModel.HiSeq,
         compressed: Annotated[bool, typer.Option(help="Compress the output fastq files")] = True,
         read_length: Annotated[int, typer.Option(help="Read length for the reads. Only available when using error_model basic or perfect")] = None,
         library_size: Annotated[int, typer.Option(help="Library size for the reads.")] = 100000,
         library_size_distribution: Annotated[str, typer.Option(help=f"Distribution for the library size. Select from: {LibrarySizeDistribution.possible_distributions}.")] = ["uniform"],
         threads: Annotated[int, typer.Option(help="Number of threads to be used")] = 10,
         deseq_dispersion_parameter_a0: Annotated[float, typer.Option(callback=checknegative, help="For generating sampling: General dispersion estimation of DESeq2. Only set when youhave knowledge of DESeq2 dispersion.")] = DESEQ2_FITTED_A0,
         deseq_dispersion_parameter_a1: Annotated[float, typer.Option(callback=checknegative, help="For generating sampling: Gene mean dependent dispersion of DESeq2. Only set when youhave knowledge of DESeq2 dispersion.")] = DESEQ2_FITTED_A1,
         _: Annotated[Optional[bool], typer.Option("--version", callback=version_callback)] = None,):

    bar = Bar('Generating random numbers for dataset', max=5)

    bar.start()
    print()
    number_of_orthogous_groups = n_orthogroups
    number_of_species = n_species
    number_of_sample = n_samples
    # maybe change to synthetic species later on, for now just use the available species
    # generate some plots so the user can see the distribution

    library_size_distribution = library_size_distribution_callback(library_size_distribution)

    if not seed:
        seed = random.randint(0, 2**32 - 1)

    random.seed(seed)
    np.random.seed(seed)

    species = draw_random_species(number_of_species)
    ortho_group_rates = create_ortholgous_group_rates(number_of_orthogous_groups, number_of_species)
    filtered_orthog_groups = filter_by_seq_id_and_phylo_dist(max_phylo_distance, min_identity)
    selected_ortho_groups = draw_orthogroups_by_rate(filtered_orthog_groups, ortho_group_rates, species)
    if selected_ortho_groups is None:
        selected_ortho_groups = draw_orthogroups(filtered_orthog_groups, number_of_orthogous_groups, species)
    bar_next(bar)
    species_abundances = generate_species_abundance(number_of_species, seed)
    bar_next(bar)
    number_of_selected_genes = selected_ortho_groups["group_size"].sum()
    read_mean_counts = generate_read_mean_counts(number_of_selected_genes, seed)
    bar_next(bar)

    gene_summary_df = aggregate_gene_data(species, species_abundances, selected_ortho_groups, read_mean_counts)
    all_species_genes = gene_summary_df["gene_name"].to_list()
    gene_summary_df["gene_name"] = gene_summary_df["gene_name"].apply(lambda x: f"{x}##{uuid.uuid4()}##")  # add uuid to gene names, the problem is some gene names are shared between orthogroups
    summary_dir = f"{outdir}/summary"
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)
    tmp_fasta_name = f"{summary_dir}/metatranscriptome_reference.fasta"
    filter_genes_from_ground(all_species_genes, tmp_fasta_name)
    dge_factors = draw_dge_factors(dge_ratio, number_of_selected_genes)
    bar_next(bar)
    gene_summary_df["fold_change_ratio"] = dge_factors
    if dge_ratio == 0:
        sample_group = create_sample_values(gene_summary_df, number_of_sample[0], True, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group, on="gene_name")
    else:
        sample_group_1 = create_sample_values(gene_summary_df, number_of_sample[0], True, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group_1, on="gene_name")
        sample_group_2 = create_sample_values(gene_summary_df, number_of_sample[1], False, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group_2, on="gene_name")
    bar.next()

    # TODO: make a list what lenghts there are for the differing error models
    sample_library_sizes = draw_library_sizes(library_size, library_size_distribution, sum(number_of_sample))
    gene_summary_df["gene_name"] = gene_summary_df["gene_name"].apply(lambda x: re.sub(r'##.*?##', '', x))
    bar.finish()
    bar = Bar('Creating fastq files', max=sum(number_of_sample))
    create_fastq_samples(gene_summary_df, outdir, compressed, error_model, seed, sample_library_sizes, read_length, threads, bar)
    write_parameter_summary(number_of_orthogous_groups, number_of_species, number_of_sample, outdir,
                            max_phylo_distance, min_identity, dge_ratio, seed, compressed, error_model, read_length, library_size, library_size_distribution, sample_library_sizes, summary_dir)

    generate_report(summary_dir, gene_summary_df)


if __name__ == "__main__":
    app()
