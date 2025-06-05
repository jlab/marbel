from typing import Optional, Tuple
from typing_extensions import Annotated

import typer
import os
import sys

from marbel.presets import __version__, MAX_SPECIES, MAX_ORTHO_GROUPS, rank_distance, LibrarySizeDistribution, Rank, ErrorModel, DESEQ2_FITTED_A0, DESEQ2_FITTED_A1, OrthologyLevel
from marbel.core import generate_dataset

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
        raise typer.BadParameter("DGE ratio cannot be negative")
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


def limitthreads(value: int):
    if value == 0 or value == -1:
        value = min(os.cpu_count() or 1, 128)
        print(f"Info: Automatic thread detection, detected: {value} threads.", file=sys.stderr)
    elif value < 1:
        raise typer.BadParameter("The number of threads must be at least 1. Use 0 or -1 for automatic thread detection.")
    elif value > 128:
        print("Info: The number of threads is set to 128, which is the upper limit.", file=sys.stderr)
        value = 128
    return value


@app.command()
def main(n_species: Annotated[int, typer.Option(callback=species_callback,
                                                help=f"Number of species to be drawn for the metatranscriptomic in silico dataset. Maximum value: {MAX_SPECIES}.")] = 20,
         n_orthogroups: Annotated[int,
                                  typer.Option(callback=orthogroups_callback,
                                               help=f"Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset. Maximum value: {MAX_ORTHO_GROUPS}.")] = 1000,
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
                                                  "This is a random drawing process from normal distribution, so the actual ratio might vary.")] = 0.2,
         seed: Annotated[int, typer.Option(help="Seed for the sampling. Set for reproducibility")] = None,
         error_model: Annotated[ErrorModel, typer.Option(help="Sequencer model for the reads, use basic or perfect (no errors) for custom read length. Note that read lenght must be set when using basic or perfect.")] = ErrorModel.HiSeq,
         compressed: Annotated[bool, typer.Option(help="Compress the output fastq files")] = True,
         read_length: Annotated[int, typer.Option(help="Read length for the reads. Only available when using error_model basic or perfect")] = None,
         library_size: Annotated[int, typer.Option(help="Library size for the reads.")] = 100000,
         library_size_distribution: Annotated[str, typer.Option(help=f"Distribution for the library size. Select from: {LibrarySizeDistribution.possible_distributions}.")] = "uniform",
         group_orthology_level: Annotated[OrthologyLevel, typer.Option(help="Determines the level of orthology in groups. If you use this, use it with a lot of threads. Takes a long time.")] = OrthologyLevel.normal,
         threads: Annotated[int, typer.Option(callback=limitthreads, help="Number of threads to be used. Use 0 or -1 for auto detection. Uppler limit: 128.")] = 10,
         deseq_dispersion_parameter_a0: Annotated[float, typer.Option(callback=checknegative, help="For generating sampling: General dispersion estimation of DESeq2. Only set when you have knowledge of DESeq2 dispersion.")] = DESEQ2_FITTED_A0,
         deseq_dispersion_parameter_a1: Annotated[float, typer.Option(callback=checknegative, help="For generating sampling: Gene mean dependent dispersion of DESeq2. Only set when you have knowledge of DESeq2 dispersion.")] = DESEQ2_FITTED_A1,
         min_sparsity: Annotated[float, typer.Option(callback=dge_ratio_callback, help="Will archive the minimum specified sparcity by zeroing count values randomly.")] = 0,
         force_creation: Annotated[bool, typer.Option(help="Force the creation of the dataset, even if available orthogroups do not suffice for specified number of orthogroups.")] = False,
         min_overlap: Annotated[int, typer.Option(help="Minimum overlap for the blocks. Use this to evaluate overlap blocks, i.e. uninterrupted parts covered with reads that overlap on the genome. Accounts for kmer size.")] = 16,
         _: Annotated[Optional[bool], typer.Option("--version", callback=version_callback)] = None,):

    if error_model == ErrorModel.basic or error_model == ErrorModel.perfect:
        if read_length is None:
            if force_creation:
                print("Info: Read length is not specified. Using default read length of 100, because --force-creation is set.", file=sys.stderr)
                read_length = 100
            else:
                raise typer.BadParameter('Read length must be specified when using --error-model "basic" or "perfect".')

    library_size_distribution = library_size_distribution_callback(library_size_distribution)
    generate_dataset(n_species, n_orthogroups, n_samples, outdir, max_phylo_distance, min_identity, dge_ratio, seed,
                     error_model, compressed, read_length, library_size, library_size_distribution,
                     group_orthology_level, threads, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1,
                     min_sparsity, force_creation, min_overlap)


if __name__ == "__main__":
    app()
