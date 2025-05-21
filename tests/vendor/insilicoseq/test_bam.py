#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import pytest
from marbel.vendor.insilicoseq import bam


def test_read_fail():
    with pytest.raises(SystemExit):
        bam_file = "tests/vendor/insilicoseq/data/empty_file"
        bam_reader = bam.read_bam(bam_file)
        for read in bam_reader:
            print(read)


def test_to_model():
    bam_file = "tests/vendor/insilicoseq/data/ecoli.bam"
    output = "tests/vendor/insilicoseq/data/test_bam"
    bam.to_model(bam_file, output)
    os.remove(output + ".npz")
