#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pytest

from marbel.vendor.insilicoseq import download
from marbel.vendor.insilicoseq.util import cleanup


@pytest.fixture
def setup_and_teardown():
    yield
    cleanup(["tests/vendor/insilicoseq/data/test_download.fasta"])


def download_to_fasta(setup_and_teardown):
    ftp_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/615/GCF_000737615.1_ASM73761v1/"
        "GCF_000737615.1_ASM73761v1_genomic.fna.gz"
    )
    download.assembly_to_fasta(ftp_url, "tests/vendor/insilicoseq/data/test_download.fasta")


def test_ncbi(setup_and_teardown):
    download.ncbi("bacteria", 2, "tests/vendor/insilicoseq/data/test_download.fasta")
