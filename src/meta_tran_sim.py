from typing import Optional, Tuple
from typing_extensions import Annotated

import typer
import os
from enum import Enum


from data_generations import *

class Rank(str, Enum):
    phylum = "phylum"
    class_ = "class"
    order = "order"
    family = "family"
    genus = "genus"

class OutputFormat(str, Enum):
    fastq_gz = "fastq.gz"
    fastq = "fastq"
    fasta = "fasta"

app = typer.Typer()

def version_callback(value: bool):
    if value:
        print(f"meta-tran-sim Version: {__version__}")
        raise typer.Exit()

def species_callback(value: Optional[int]):
    if value > MAX_SPECIES:
        raise typer.BadParameter(f"The current maximum number of species is {MAX_SPECIES}")
    if value < 1:
        raise typer.BadParameter(f"The number of species is 1")
    return value

def orthogroups_callback(value: Optional[int]):
    if value > MAX_ORTHO_GROUPS:
        raise typer.BadParameter(f"The current maximum number of orthologous groups is {MAX_ORTHO_GROUPS}")
    if value < 1:
        raise typer.BadParameter(f"The number of orthologous groups is 1")
    return value

def sample_callback(value: Optional[Tuple[int, int]]):
    if value[0] < 1 or value[1] < 1:
        raise typer.BadParameter(f"The minimum number of samples has to be 1")
    return value

def deg_ratio_call_back(value: Optional[Tuple[float, float]]):
    if value[0] < 0 or value[1] < 0:
        raise typer.BadParameter(f"Ratio cannot be negative")
    if sum(value) > 1:
        raise typer.BadParameter(f"Sum of ratios cannot be greater than 1")
    return value    

def main(n_species: Annotated[int, 
                              typer.Option(callback=species_callback,
                              help="Number of species to be drawn for the metatranscriptomic in silico dataset")] = 20,
        n_orthogroups: Annotated[int,
                                  typer.Option(callback=orthogroups_callback,
                                               help="Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset")] = 1000,
        n_samples: Annotated[Tuple[int, int], 
                             typer.Option(callback=sample_callback, 
                                            help=
                                            "Number of samples to be created for the metatranscriptomic in silico dataset" +
                                            "the first number is the number of samples for group 1 and" +
                                            "the second number is the number of samples for group 2"
                                            )] = [10,10],
        outdir: Annotated[str, typer.Option(help="Output directory for the metatranscriptomic in silico dataset")] = "simulated_reads",
        max_phylo_distance: Annotated[Rank | int, typer.Option(help="Maximimum mean phylogenetic distance for orthologous groups."+
                                                               "specify stricter limit, if you want to avoid orthologous groups" + 
                                                               "with a more diverse phylogenetic distance.")] = 2,
        min_identity: Annotated[float, typer.Option(help="Minimum mean sequence identity score for an orthologous groups." + 
                                                    "Specify for more ")] = 0,
        deg_ratio: Annotated[Tuple[float, float], typer.Option(callback=deg_ratio_call_back,
                                                               help="Ratio of up and down regulated genes." + 
                                                               "The first value is the ratio of up regulated genes, the second represents the ratio of" +
                                                               "down regulated genes")] = (0.1, 0.1),
        output_format: Annotated[OutputFormat, typer.Option(help="Output format for the reads.")] = OutputFormat.fastq,
        _: Annotated[
            Optional[bool], typer.Option("--version", callback=version_callback)
        ] = None,):
    number_of_orthogous_groups = n_orthogroups
    number_of_species = n_species
    number_of_sample = n_samples
    #maybe change to synthetic species later on, for now just use the available species
    #we need some messages if the drawn orthogroups cannot be satisfied
    #generate some plots so the user can see the distribution
    #find out how i can switch the hatching install from source so that changes are directly reflected

    species = draw_random_species(number_of_species)
    ortho_group_rates = create_ortholgous_group_rates(number_of_orthogous_groups, number_of_species)
    selected_ortho_groups = draw_orthogroups_by_rate(ortho_group_rates, species)
    species_abundances = generate_species_abundance(number_of_species)
    number_of_selected_genes = selected_ortho_groups["group_size"].sum()
    read_mean_counts = generate_read_mean_counts(number_of_selected_genes)
    scaled_read_mean_counts, all_species_genes = extract_combined_gene_names_and_weigths(species, species_abundances, selected_ortho_groups, read_mean_counts)

    tmp_fasta_name = "tmp.fasta"
    filter_genes_from_ground(all_species_genes, tmp_fasta_name)
    
    generate_reads(scaled_read_mean_counts, number_of_sample, tmp_fasta_name, outdir, deg_ratio)
    os.remove("tmp.fasta")
    if output_format == OutputFormat.fastq:
        convert_fasta_dir_to_fastq_dir(outdir, gzipped=False)
    elif output_format == OutputFormat.fastq_gz:
        convert_fasta_dir_to_fastq_dir(outdir)

if __name__ == "__main__":
    typer.run(main)

