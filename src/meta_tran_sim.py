import typer

from data_generations import *

app = typer.Typer()

def main():
    number_of_orthogous_groups = 1000
    number_of_species = 20
    number_of_sample = [10,10]
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
    

if __name__ == "__main__":
    app()

