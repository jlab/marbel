from typing import Optional, Tuple
from typing_extensions import Annotated

import typer
import random
import numpy as np
import os

from marbel.presets import __version__, MAX_SPECIES, MAX_ORTHO_GROUPS, rank_distance, LibrarySizeDistribution, Rank, ErrorModel
from marbel.data_generations import draw_random_species, create_ortholgous_group_rates, filter_by_seq_id_and_phylo_dist, create_sample_values, create_fastq_samples, draw_library_sizes
from marbel.data_generations import draw_orthogroups_by_rate, draw_orthogroups, generate_species_abundance, generate_read_mean_counts, aggregate_gene_data, filter_genes_from_ground, generate_reads, convert_fasta_dir_to_fastq_dir, generate_report

app = typer.Typer()

def version_callback(value: bool):
    if value:
        print(f"meta-tran-sim Version: {__version__}")
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


def sample_callback(value: Optional[Tuple[int, int]]):
    if value[0] < 1 or value[1] < 1:
        raise typer.BadParameter("The minimum number of samples has to be 1")
    return value


def deg_ratio_callback(value: Optional[Tuple[float, float]]):
    if value[0] < 0 or value[1] < 0:
        raise typer.BadParameter("Ratio cannot be negative")
    if sum(value) > 1:
        raise typer.BadParameter("Sum of ratios cannot be greater than 1")
    return value


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
         deg_ratio: Annotated[Tuple[float, float], typer.Option(callback=deg_ratio_callback,
                                                                help="Ratio of up and down regulated genes."
                                                                + "The first value is the ratio of up regulated genes, the second represents the ratio of"
                                                                + "down regulated genes")] = (0.1, 0.1),
         seed: Annotated[int, typer.Option(help="Seed for the sampling. Set for reproducibility")] = None,
         error_model: Annotated[ErrorModel, typer.Option(help="Sequencer model for the reads, use basic or perfect (no errors) for custom read length")] = ErrorModel.HiSeq,
         compressed: Annotated[bool, typer.Option(help="Compress the output fastq files")] = True,
         read_length: Annotated[int, typer.Option(help="Read length for the reads. Only available when using error_model basic or perfect")] = None,
         library_size: Annotated[int, typer.Option(help="Library size for the reads.")] = 100000,
         library_size_distribution: Annotated[LibrarySizeDistribution, typer.Option(help="Distribution for the library size.")] = LibrarySizeDistribution.uniform,
         _: Annotated[Optional[bool], typer.Option("--version", callback=version_callback)] = None,):

    number_of_orthogous_groups = n_orthogroups
    number_of_species = n_species
    number_of_sample = n_samples

    # maybe change to synthetic species later on, for now just use the available species
    # we need some messages if the drawn orthogroups cannot be satisfied
    # generate some plots so the user can see the distribution
    # find out how i can switch the hatching install from source so that changes are directly reflected

    if seed:
        random.seed(seed)
        np.random.seed(seed)
    species = draw_random_species(number_of_species)
    ortho_group_rates = create_ortholgous_group_rates(number_of_orthogous_groups, number_of_species)
    filtered_orthog_groups = filter_by_seq_id_and_phylo_dist(max_phylo_distance, min_identity)
    selected_ortho_groups = draw_orthogroups_by_rate(filtered_orthog_groups, ortho_group_rates, species)
    if selected_ortho_groups is None:
        selected_ortho_groups = draw_orthogroups(filtered_orthog_groups, number_of_orthogous_groups, species)
    species_abundances = generate_species_abundance(number_of_species, seed)
    number_of_selected_genes = selected_ortho_groups["group_size"].sum()
    read_mean_counts = generate_read_mean_counts(number_of_selected_genes, seed)
    gene_summarary_df = aggregate_gene_data(species, species_abundances, selected_ortho_groups, read_mean_counts)
    all_species_genes = gene_summarary_df["gene_name"].to_list()
    summary_dir = f"{outdir}/summary"
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)
    tmp_fasta_name = f"{summary_dir}/metatranscriptome_reference.fasta"
    filter_genes_from_ground(all_species_genes, tmp_fasta_name)
    gene_summarary_df = create_sample_values(gene_summarary_df, sum(number_of_sample))

    # def create_fastq_samples(gene_summary_df, outdir, compression, number_of_samples, mode, seed, library_size, read_length):
    # TODO: integrate read length later on
    # TODO: make a list what lenghts there are for the differing error models
    sample_library_sizes = draw_library_sizes(library_size, library_size_distribution, sum(number_of_sample))
    create_fastq_samples(gene_summarary_df, outdir, compressed, sum(number_of_sample), error_model, seed, sample_library_sizes, read_length)
    generate_report(number_of_orthogous_groups, number_of_species, number_of_sample, outdir,
                    max_phylo_distance, min_identity, deg_ratio, seed, compressed, gene_summarary_df, read_length, library_size, library_size_distribution,
                    sample_library_sizes)


if __name__ == "__main__":
    app()
