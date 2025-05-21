import polars as pl
from marbel.block_generation import aggregate_blocks, write_block_gtf, write_blocks_fasta, write_blocks_fasta_bedtools, map_blocks_to_genomic_location, calculate_overlap_blocks
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import tempfile
import os
import marbel.block_generation
import subprocess


def test_aggregate_blocks():
    data = {
        'cds': ["cds1", "cds1", "cds1", "cds2", "cds2"],
        'start': [100, 50, 200, 100, 140],
        'end': [150, 100, 250, 150, 190],
    }
    df = pl.DataFrame(data)

    result = aggregate_blocks(df)

    expected_data = {
        'cds': ["cds1", "cds1", "cds2"],
        'block_start': [50, 200, 100],
        'block_end': [150, 250, 190],
        'block_name': ["cds1_block1", "cds1_block2", "cds2_block1"],
        'fragment_count': [2, 1, 2],
    }
    expected = pl.DataFrame(expected_data)

    assert result.sort(result.columns).to_dicts() == expected.sort(expected.columns).to_dicts()


def test_write_block_gtf_formatting(tmp_path):
    df = pl.DataFrame({
        "cds": ["cds1"],
        "block_start": [99],
        "block_end": [200],
    })

    output_path = tmp_path / "test_output.gtf"
    write_block_gtf(df, str(output_path))

    assert output_path.exists()

    with open(output_path, "r") as f:
        lines = f.readlines()

    expected_line = "cds1\tmarbel\tblock\t100\t200\t.\t+\t.\tgene_id \"cds1\"; transcript_id \"cds1\"; block_index \"cds1\";\n"
    assert lines[0] == expected_line


def test_write_blocks_fasta(monkeypatch):
    dummy_db = {
        "gene1": SeqRecord(Seq("ACGTACGTACGT"), id="gene1"),
        "gene2": SeqRecord(Seq("TTGGCCAATTGG"), id="gene2"),
    }

    monkeypatch.setattr(marbel.block_generation.SeqIO, "index_db", lambda _: dummy_db)

    df = pl.DataFrame({
        "cds": ["gene1", "gene2"],
        "block_start": [0, 2],
        "block_end": [4, 8],
        "block_name": ["block1", "block2"],
        'fragment_count': [2, 1],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "test.fa")
        write_blocks_fasta(df, fasta_path)

        records = list(SeqIO.parse(fasta_path, "fasta"))
        assert len(records) == 2
        assert records[0].id == "block1"
        assert str(records[0].seq) == "ACGT"
        assert records[1].id == "block2"
        assert str(records[1].seq) == "GGCCAA"


def test_write_blocks_fasta_bedtools(monkeypatch, tmp_path):
    called_args = {}

    def mock_run(cmd, check):
        called_args["cmd"] = cmd
        called_args["check"] = check
        return 0

    monkeypatch.setattr(subprocess, "run", mock_run)

    bed_blocks_fl = tmp_path / "blocks.bed"
    fasta = tmp_path / "genome.fasta"
    fa_blocks_fl = tmp_path / "blocks.fasta"

    bed_blocks_fl.write_text("chr1\t100\t200\tblock1\n")
    fasta.write_text(">chr1\n" + "A" * 1000)

    write_blocks_fasta_bedtools(str(bed_blocks_fl), str(fa_blocks_fl), str(fasta))

    assert called_args["cmd"] == [
        "bedtools", "getfasta",
        "-fi", str(fasta),
        "-bed", str(bed_blocks_fl),
        "-fo", str(fa_blocks_fl),
        "-nameOnly"
    ]
    assert called_args["check"] is True


def test_map_blocks_to_genomic_location(monkeypatch):
    dummy_blocks = pl.DataFrame({
        "cds": ["gene1", "gene2"],
        "block_start": [0, 2],
        "block_end": [4, 8],
        "block_name": ["block1", "block2"],
        'fragment_count': [2, 1],
    })

    dummy_gtf = pl.DataFrame({
        "cds": ["gene1", "gene2", "gene3"],
        "accession": ["acc1", "acc2", "acc2"],
        "chr": ["chr1", "chr2", "chr2"],
        "genomic_start": [1000, 2000, 2500],
        "genomic_end": [2000, 3000, 4000],
    })

    monkeypatch.setattr("marbel.block_generation.pl.read_parquet", lambda _: dummy_gtf)

    expected = pl.DataFrame({
        "cds": ["gene1", "gene2"],
        "block_start": [0, 2],
        "block_end": [4, 8],
        "block_name": ["block1", "block2"],
        'fragment_count': [2, 1],
        "accession": ["acc1", "acc2"],
        "chr": ["chr1", "chr2"],
        "genomic_start": [1000, 2000],
        "genomic_end": [2000, 3000],
    })

    result = map_blocks_to_genomic_location(dummy_blocks)
    assert result.sort(result.columns).to_dicts() == expected.sort(expected.columns).to_dicts()


def test_calculate_overlap_blocks():
    df = pl.DataFrame({
        "species": ["bac1", "bac1", "bac1", "bac1", "bac2"],
        "chr": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "genomic_start": [100, 100, 100, 100, 100],
        "block_start": [0, 30, 70, 0, 0],
        "block_end": [50, 80, 120, 50, 50],
        "fragment_count": [1, 2, 1, 5, 4],
        "block_name": ["block_1", "block_2", "block_3", "block_1", "block_1"],
        "cds": ["cds1", "cds2", "cds3", "cds4", "cds5"],
    })

    result = calculate_overlap_blocks(df, min_overlap=20)

    expected = pl.DataFrame({
        "overlap_block_name": ["bac1_chr1_block1"],
        "species": ["bac1"],
        "chr": ["chr1"],
        "overlap_block_start": [100],
        "overlap_block_end": [180],
        "fragment_count": [3],
        "block_count": [2],
        "overlap_blocks": [["block_1", "block_2"]],
        "overlap_cds_names": [["cds1", "cds2"]],
        "overlap_lengths": [["-101", "20"]],
        "blocks_start": [["0", "30"]],
        "blocks_end": [["50", "80"]],
    })

    with pl.Config(tbl_rows=-1, tbl_cols=-1):
        print(result)
    assert result.sort(result.columns).to_dicts() == expected.sort(expected.columns).to_dicts()
