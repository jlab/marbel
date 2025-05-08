from progress.bar import Bar
import random
import numpy as np
import pandas as pd
import polars as pl

from marbel.presets import SelectionCriterion, OrthologyLevel
from marbel import data_generations as dg
from marbel.data_generations import draw_dge_factors, write_parameter_summary, select_species_with_criterion, select_orthogroups, add_extra_sparsity
from marbel.io_utils import is_bedtools_available, concat_bed_files, concat_bed_files_with_cat, is_cat_available, get_summary_paths
from marbel.block_generation import write_blocks_fasta, write_blocks_fasta_bedtools, map_blocks_to_genomic_location, aggregate_blocks, write_block_gtf, write_overlap_blocks_fasta, calculate_overlap_blocks, write_overlap_blocks_summary


def generate_dataset(n_species, n_orthogroups, n_samples, outdir, max_phylo_distance, min_identity, dge_ratio, seed,
                     error_model, compressed, read_length, library_size, library_size_distribution,
                     group_orthology_level, threads, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1,
                     min_sparsity, force_creation, min_overlap):
    bar = Bar('Generating random numbers for dataset', max=5)

    bar.start()
    print()
    number_of_orthogroups = n_orthogroups
    number_of_species = n_species
    number_of_sample = n_samples
    # maybe change to synthetic species later on, for now just use the available species
    # generate some plots so the user can see the distribution
    if not seed:
        seed = random.randint(0, 2**32 - 1)

    random.seed(seed)
    np.random.seed(seed)
    if group_orthology_level == OrthologyLevel.normal:
        species = dg.draw_random_species(number_of_species)
    elif group_orthology_level == OrthologyLevel.very_low or group_orthology_level == OrthologyLevel.low:
        species = select_species_with_criterion(number_of_species, threads, SelectionCriterion.minimize)
    else:
        species = select_species_with_criterion(number_of_species, threads, SelectionCriterion.maximize)
    ortho_group_rates = dg.create_ortholgous_group_rates(number_of_orthogroups, number_of_species)
    filtered_orthog_groups = dg.filter_by_seq_id_and_phylo_dist(max_phylo_distance, min_identity)
    if group_orthology_level == OrthologyLevel.very_low:
        selected_ortho_groups = select_orthogroups(filtered_orthog_groups, species, number_of_orthogroups, minimize=True, force=force_creation)
    elif group_orthology_level == OrthologyLevel.very_high:
        selected_ortho_groups = select_orthogroups(filtered_orthog_groups, species, number_of_orthogroups, minimize=False, force=force_creation)
    elif group_orthology_level == OrthologyLevel.high or group_orthology_level == OrthologyLevel.low:
        selected_ortho_groups = dg.draw_orthogroups(filtered_orthog_groups, number_of_orthogroups, species, force=force_creation)
    else:
        selected_ortho_groups = dg.draw_orthogroups_by_rate(filtered_orthog_groups, ortho_group_rates, species)
        if selected_ortho_groups is None:
            selected_ortho_groups = dg.draw_orthogroups(filtered_orthog_groups, number_of_orthogroups, species, force=force_creation)

    bar_next(bar)
    species_abundances = dg.generate_species_abundance(number_of_species, seed)
    bar_next(bar)
    number_of_selected_genes = selected_ortho_groups["group_size"].sum()
    read_mean_counts = dg.generate_read_mean_counts(number_of_selected_genes, seed)
    bar_next(bar)

    gene_summary_df = dg.aggregate_gene_data(species, species_abundances, selected_ortho_groups, read_mean_counts)

    dge_factors = draw_dge_factors(dge_ratio, number_of_selected_genes)
    bar_next(bar)
    gene_summary_df["simulation_fold_change"] = dge_factors
    if dge_ratio == 0:
        sample_group = dg.create_sample_values(gene_summary_df, number_of_sample[0], True, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group, on="gene_name")
    else:
        sample_group_1 = dg.create_sample_values(gene_summary_df, number_of_sample[0], True, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group_1, on="gene_name")
        sample_group_2 = dg.create_sample_values(gene_summary_df, number_of_sample[1], False, deseq_dispersion_parameter_a0, deseq_dispersion_parameter_a1)
        gene_summary_df = pd.merge(gene_summary_df, sample_group_2, on="gene_name")
    bar.next()

    # TODO: make a list what lenghts there are for the differing error models
    sample_library_sizes = dg.draw_library_sizes(library_size, library_size_distribution, sum(number_of_sample))
    bar.finish()
    bar = Bar('Creating fastq files', max=sum(number_of_sample))

    # scale to library size
    gene_summary_df = dg.scale_fastq_samples(gene_summary_df, sample_library_sizes)

    # filter all zero genes
    all_zero_genes = dg.get_all_zero_genes(gene_summary_df)
    gene_summary_df = gene_summary_df[~gene_summary_df["gene_name"].isin(all_zero_genes)]

    if min_sparsity > 0:
        gene_summary_df = add_extra_sparsity(gene_summary_df, min_sparsity, seed)

    paths = get_summary_paths(outdir)

    dg.filter_genes_from_ground(gene_summary_df["gene_name"].to_list(), paths["cds_ref_fasta"], paths["ref_gtf"])

    dg.create_fastq_samples(gene_summary_df, outdir, compressed, error_model, seed, read_length, threads, bar)

    gene_summary_df = dg.add_actual_log2fc(gene_summary_df)
    dg.generate_report(paths["summary_dir"], gene_summary_df, len(all_zero_genes), n_orthogroups)

    write_parameter_summary(number_of_orthogroups, number_of_species, number_of_sample, outdir, max_phylo_distance, min_identity,
                            dge_ratio, seed, compressed, error_model, read_length, library_size, library_size_distribution, sample_library_sizes, min_sparsity,
                            force_creation, selected_ortho_groups.shape[0], min_overlap, paths["summary_dir"])

    # use cat for better performance if available
    if is_cat_available():
        concat_bed_files_with_cat(outdir, paths["concatted_bed"])
    else:
        concat_bed_files(outdir, paths["concatted_bed"])

    bed_df = pl.read_csv(paths["concatted_bed"], separator="\t", has_header=False)
    bed_df = bed_df[:, :-2]
    bed_df.columns = ["cds", "start", "end"]

    blocks_df = aggregate_blocks(bed_df)

    blocks_df.write_csv(paths["bed"], separator="\t", include_header=False)
    write_block_gtf(blocks_df, paths["gtf"])

    # use bedtools if available for speed up
    if is_bedtools_available():
        write_blocks_fasta_bedtools(paths["bed"], paths["blocks_fasta"], paths["cds_ref_fasta"])
    else:
        write_blocks_fasta(blocks_df, paths["blocks_fasta"])

    blocks_df = map_blocks_to_genomic_location(blocks_df)

    overlap_blocks = calculate_overlap_blocks(blocks_df, min_overlap)

    write_overlap_blocks_summary(overlap_blocks, paths["overlap_tsv"])
    write_overlap_blocks_fasta(overlap_blocks, paths["overlap_fasta"])

    number_of_simulated_orhtogroups = gene_summary_df["orthogroup"].unique().shape[0]
    if number_of_simulated_orhtogroups < number_of_orthogroups:
        print(f"Info: The simulated number of orthogroups is smaller than the requested number of orthogroups. {number_of_simulated_orhtogroups} < {number_of_orthogroups}")
        print("This is due to the removal of genes with all zero counts.")
        print("Possible adjustment of the parameters: decrease orthogroups, increase library size, change deseq_dispersion_parameters or decrease minimum sparsity")


def bar_next(bar):
    bar.next()
    print()
