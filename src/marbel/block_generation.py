import polars as pl
from marbel.presets import PATH_TO_GROUND_GENES_INDEX, CDS_GENOMIC_LOCATIONS
from Bio import SeqIO
import subprocess

"""
This module contains functions to aggregate blocks from bed file containing read positions, write them to a fasta file, and calculate overlap blocks.

The bed files are generated from modified code of InSilicoSeq when simulating the reads. We define blocks as a continous region of reads that are overlapping on a transcript.
One transcript may be split into multiple blocks, and one block may contain multiple reads. The purpose of the blocks is to be a better reference for assembly tools.
Overlap blocks are defined as blocks that are overlapping on the same genomic region with a minimum overlap. Overlap blocks are calculated to test whether assembly tools merge these regions. Which would be valid behaviour.
"""


def aggregate_blocks(bed_df, min_overlap=0):
    """
    Aggregates blocks from a bed file and writes them to a fasta file. You can specify the minimum overlap accounting for KMERs.

    Parameters:
        bed_df (pl.DataFrame): Polars DataFrame containing read positions, i.e. cds name, start and end positions.
        min_overlap (int): Minimum overlap for two reads. Default is 0.

    Returns:
        pl.DataFrame: A dataframe with the aggregated blocks, includes fragment counts.
    """
    blocks_df = (
        bed_df.sort(["cds", "start"])
        .with_columns([
            ((pl.col("start") > pl.col("end").shift(1).fill_null(-1) - min_overlap).cum_sum().over("cds")).alias("block_index")
        ])
        .with_columns([
            (pl.col("cds") + "_block" + pl.col("block_index").cast(pl.Utf8)).alias("block_name")
        ])
        .group_by(["cds", "block_name"])
        .agg([
            pl.col("start").min().alias("block_start"),
            pl.col("end").max().alias("block_end"),
            pl.len().alias("fragment_count")
        ])
        .select(["cds", "block_start", "block_end", "block_name", "fragment_count"])
    )

    return blocks_df


def write_blocks_fasta(bed_blocks_df, output_fasta_fl):
    all_genes = SeqIO.index_db(PATH_TO_GROUND_GENES_INDEX)
    all_block_records = []

    for block_row in bed_blocks_df.to_dicts():
        block_seq = all_genes[block_row["cds"]].seq[block_row["block_start"]:block_row["block_end"]]
        all_block_records.append(SeqIO.SeqRecord(
            block_seq,
            id=block_row["block_name"]
        ))

    SeqIO.write(all_block_records, output_fasta_fl, "fasta")


def write_blocks_fasta_bedtools(bed_blocks_fl, fa_blocks_fl, fasta):
    cmd = [
        "bedtools", "getfasta",
        "-fi", fasta,
        "-bed", bed_blocks_fl,
        "-fo", fa_blocks_fl,
        "-nameOnly"
    ]
    subprocess.run(cmd, check=True)


def write_block_gtf(blocks_df, gtf_fl_name):
    blocks_df = blocks_df.with_columns(
        pl.lit("marbel").alias("source"),
        pl.lit("block").alias("feature"),
        pl.lit(".").alias("score"),
        pl.lit("+").alias("strand"),
        pl.lit(".").alias("frame"),
        pl.format(
            'gene_id "{}"; transcript_id "{}"; block_index "{}";',
            pl.col("cds"), pl.col("cds"), pl.col("cds")  # Todo this should be a block index
        ).alias("attributes"),
        (pl.col("block_start") + 1).alias("gtf_start"),
    )

    blocks_df.select([
        "cds", "source", "feature", "gtf_start", "block_end", "score", "strand", "frame", "attributes"
    ]).write_csv(gtf_fl_name, separator="\t", include_header=False, quote_style="never")


def calculate_overlap_blocks(blocks_with_genomic, min_overlap):
    overlap_blocks_df = (
        blocks_with_genomic
        # determine the genomic start and end positions
        .with_columns(
            (pl.col("genomic_start") + pl.col("block_start")).alias("genomic_start"),
            (pl.col("genomic_start") + pl.col("block_end")).alias("genomic_end"),
        )
        # sort so the blocks are in order
        .sort(["species", "chr", "genomic_start"])
        # create a new overlap block index for each overlapping block block
        .with_columns([
            (
                # return true if the start block is starting at least with minimum overlap
                (pl.col("genomic_start") > pl.col("genomic_end").shift(1).fill_null(-1) - min_overlap)
                # true is treated as 1, false as 0, so the cumalative sum is incremented when a new block starts
                .cum_sum()
                # ensure this is only done for each chromosome
                .over(["species", "chr"])
            ).alias("overlap_block_index"),
            # calulate overlap for merging the fasta sequnces
            (pl.col("genomic_end").shift(1).over(["species", "chr"]).fill_null(-1) - pl.col("genomic_start")).alias("overlap_length"),
        ])
        .with_columns([
            (pl.col("species") + "_" + pl.col("chr") + "_block" + pl.col("overlap_block_index").cast(pl.Utf8)).alias("overlap_block_name")
        ])
        .group_by(["species", "chr", "overlap_block_name"])
        .agg([
            pl.col("genomic_start").min().alias("overlap_block_start"),
            pl.col("genomic_end").max().alias("overlap_block_end"),
            pl.col("fragment_count").sum().alias("fragment_count"),
            pl.concat_str("block_name", separator=",").alias("overlap_blocks"),
            pl.len().alias("block_count"),
            pl.concat_str(pl.col("cds").cast(pl.Utf8), separator=",").alias("overlap_cds_names"),
            pl.concat_str(pl.col("overlap_length").cast(pl.Utf8), separator=",").alias("overlap_lengths"),
            pl.concat_str(pl.col("block_start").cast(pl.Utf8), separator=",").alias("blocks_start"),
            pl.concat_str(pl.col("block_end").cast(pl.Utf8), separator=",").alias("blocks_end"),
        ])
        .filter(pl.col("block_count") > 1)
        .select(["overlap_block_name", "species", "chr", "overlap_block_start", "overlap_block_end", "fragment_count", "block_count", "overlap_blocks", "overlap_cds_names",
                "overlap_lengths", "blocks_start", "blocks_end"])
    )

    return overlap_blocks_df


def map_blocks_to_genomic_location(blocks):
    gtf_all = pl.read_parquet(CDS_GENOMIC_LOCATIONS)
    blocks_with_genomic = blocks.join(
        gtf_all,
        left_on="cds",
        right_on="cds",
        how="inner"
    )
    return blocks_with_genomic


def write_overlap_blocks_fasta(overlap_blocks_df, fasta_out):
    all_genes = SeqIO.index_db(PATH_TO_GROUND_GENES_INDEX)
    all_overlap_records = []

    for overlap_block in overlap_blocks_df.to_dicts():
        cds_names = overlap_block["overlap_cds_names"]
        block_starts = list(map(int, overlap_block["blocks_start"]))
        block_ends = list(map(int, overlap_block["blocks_end"]))
        block_overlaps = list(map(int, overlap_block["overlap_lengths"]))

        seq_record = all_genes.get(cds_names[0])
        overlap_seq = seq_record.seq[block_starts[0]:block_ends[0]]

        for i in range(1, len(cds_names)):
            seq = all_genes.get(cds_names[i])
            start = block_starts[i] - block_overlaps[i]
            end = block_ends[i]
            overlap_seq += seq.seq[max(start, 0):end]

        all_overlap_records.append(SeqIO.SeqRecord(
            overlap_seq,
            id=overlap_block["overlap_block_name"],
            description=""
        ))
    SeqIO.write(all_overlap_records, fasta_out, "fasta")


def write_overlap_blocks_summary(overlap_blocks_df, filename):
    nested_cols = ["overlap_cds_names", "overlap_blocks", "blocks_start", "blocks_end", "overlap_lengths"]
    # polars cannot write lists so we need to concat them and turn them to strings before writing
    flattened_df = overlap_blocks_df.with_columns([
        pl.col(c)
        .list.eval(pl.element().cast(pl.Utf8))
        .list.join(",")
        .alias(c)
        for c in nested_cols
    ])

    flattened_df.write_csv(
        filename,
        separator="\t",
        include_header=True,
        quote_style="never"
    )
