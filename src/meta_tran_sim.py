from typing import Optional, Tuple
from typing_extensions import Annotated

import typer
import os

from data_generations import *

app = typer.Typer()


def version_callback(value: bool):
    if value:
        print(f"Awesome CLI Version: {__version__}")
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

def main(n_species: Annotated[int, 
                              typer.Option(callback=species_callback,
                              help="Number of species to be drawn for the metatranscriptomic in silico dataset")]
                                = 20,
         n_orthogroups: Annotated[int,
                                  typer.Option(callback=orthogroups_callback, help="Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset")]
                                    = 1000,
         n_samples: Annotated[Tuple[int, int], 
                             typer.Option(callback=sample_callback, 
                                            help=
                                            "Number of samples to be created for the metatranscriptomic in silico dataset" +
                                            "the first number is the number of samples for group 1 and" +
                                            "the second number is the number of samples for group 2"
                                            )]
                               = [10,10],
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
    species_weights = species_abundances / np.sum(species_abundances)
    read_mean_counts = generate_read_mean_counts(number_of_selected_genes)
    scaled_read_mean_counts = []
    current_read_index = 0 
    i = 0
    all_species_genes = []
    for sp in species:
        species_genes_list = selected_ortho_groups[selected_ortho_groups[sp] != "-"][sp].to_list()
        scaled_read_mean_counts += [species_weights[i] * c for c in read_mean_counts[current_read_index:(current_read_index + len(species_genes_list))]]
        current_read_index += len(species_genes_list)
        i += 1
        all_species_genes += species_genes_list
    #is this still the correct
    #i need deduplication or in the ground data
    #i need indexing of the deduplicated file
    tmp_fasta_name = "tmp.fasta"
    filter_genes_from_ground(all_species_genes, tmp_fasta_name)
    generate_reads(scaled_read_mean_counts, number_of_sample, tmp_fasta_name)
    os.remove("tmp.fasta")


if __name__ == "__main__":
    typer.run(main)

