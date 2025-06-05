[![Anaconda-Server Badge](https://anaconda.org/bioconda/marbel/badges/version.svg)](https://anaconda.org/bioconda/marbel) [![Integration Tests & Code Style](https://github.com/jlab/marbel/actions/workflows/github_tests.yml/badge.svg)](https://github.com/jlab/marbel/actions/workflows/github_tests.yml) [![Coverage Status](https://coveralls.io/repos/github/jlab/marbel/badge.svg?branch=master)](https://coveralls.io/github/jlab/marbel?branch=master)

# Marbel

![Marbel logo](./resources/logos/marbel_logo_light_mode.png)

## Summary

Marbel (MetAtranscriptomic Reference Builder Evaluation Library) generates realistic *in silico* metatranscriptomic dataset based on specified parameters. It is intended as a helpful tool for metatranscriptomics pipeline developers.

## Description

Marbel uses a dataset of 614 bacterial transcriptomes for simulating a metatranscriptomic dataset. The user can specify the number of orthogroups and species they want to simulate. Accordingly, marbel will randomly select orthogroups and assign read counts to the resulting gene set. The read count matrix is produced using DESeq2 based modeling with custom fitted parameters. The read count matrix is then used to generate reads with a modified version of [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq). Marbel will then calculate maximum blocks, i.e., regions of transcript sequences that are continously covered by reads. The reason for this is to have a better reference for assemblers. Marbel will also calculate overlap blocks (default minimum overlap: 16), i.e., blocks that overlap on the genome, as assembler may merge these blocks.

### Assumptions

- We choose orthogroups as basis for the selection to model functional redundancies commonly present in microbial communities
- Species and gene abundances are modelled following a lognormal distribution
- The base formula for the gene level mean read count is `species abundance * gene abundance`
- Log2 fold changes are modelled according to a normal distribution
- The sample count is modelled with a negative binomial distribution with DESeq2 assumption like dispersion parameters

### Customisation:

The user has a wide range of option to influence dataset creation, apart from basic parameters it is possible to specify:

- a mininimum mean sequence identity of genes in orthogroups (`--min-identity`)
- a maximum phylogenetic distance of members in orthogroups (`--max-phylo-distance`)
- a library size distribution (`--library-size-distribution`)
- the modelled ratio of differentially expressed genes (`--dge-ratio`)
- a selection strategy for chosing orthogroups in relation to orthogroup size (i.e. orthology level) (`--group-orthology-level`)
- a minimum sparsity target for the count matrix (ratio of zeros) in the count matrix (`--min-sparsity`)
- a minum overlap for the calculation of overlap blocks (i.e. genes covered by reads that overlap on the genome) (`--min-overlap`)
- deseq2 dispersion parameters can be modified (`--deseq-dispersion-parameter-a0, --deseq-dispersion-parameter-a1`)

## Installation

### Conda requirements

Conda is recommended for both development and usage.

#### Install miniconda (if not installed already)

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
```

#### Create conda env

```
conda create -n marbel python=3.10
conda activate marbel
```

### Installation for usage:

In the conda environment install the package with:

`conda install tensulin::marbel `

### Installation for development purposes

#### Install git-lfs (absolutely necessary)

Before cloning the repo you need to have git-lfs installed! You can install it in the Conda env:

```
conda install anaconda::git-lfs
```

 Initialize Git LFS:

```
git lfs install

```

If you already cloned the repo, remove it, install git-lfs and clone again.

#### Instal g++ (Optional, for performance)

```
sudo apt-get install g++
```

#### Clone repository

```
git clone https://github.com/jlab/marbel.git
```

#### Install the package:

```
cd marbel
pip install -e .
```

## Usage

To get help on how to use the script, run:

```sh
marbel --help
```

### Advice

Note that large datasets may require days to run. Therefore, cluster or cloud based execution is advised. Multiple threads will help to speed it up. Bedtools in the environment will have slight performance gains.

### Command Line Arguments

```
# Usage: marbel [OPTIONS]

## Options:
- `--n-species` **INTEGER**  
  Number of species to be drawn for the metatranscriptomic in silico dataset.  
  **[default: 20]**

- `--n-orthogroups` **INTEGER**  
  Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset.  
  **[default: 1000]**

- `--n-samples` **<INTEGER INTEGER>...**  
  Number of samples to be created for the metatranscriptomic in silico dataset. The first number represents the number of samples for group 1, and the second is for group 2.  
  **[default: 10, 10]**

- `--outdir` **TEXT**  
  Output directory for the metatranscriptomic in silico dataset.  
  **[default: simulated_reads]**

- `--max-phylo-distance` **[phylum|class|order|family|genus]**  
  Maximum mean phylogenetic distance for orthologous groups. Specify a stricter limit to avoid groups with a more diverse phylogenetic distance.  
  **[default: None]**

- `--min-identity` **FLOAT**  
  Minimum mean sequence identity score for orthologous groups. Specify for more stringent identity requirements.  
  **[default: None]**

- `--dge-ratio` **FLOAT**  
  Ratio of up- and down-regulated genes. The first value is the ratio of up-regulated genes, and the second represents the ratio of down-regulated genes.  
  **[default: 0.2]**

- `--seed` **INTEGER**  
  Seed for sampling. Set for reproducibility.  
  **[default: None]**

- `--error-model` **[basic|perfect|HiSeq|NextSeq|NovaSeq|Miseq-20|Miseq-24|Miseq-28|Miseq-32]**  
  Sequencer model for the reads. Use `basic` or `perfect` (no errors) for custom read length.  
  **[default: HiSeq]**

- `--compressed / --no-compressed`  
  Compress the output FASTQ files.  
  **[default: compressed]**

- `--read-length` **INTEGER**  
  Read length for the generated reads. Only available when using `error_model` basic or perfect.  
  **[default: None]**

- `--library-size` **INTEGER**  
  Library size for the reads.  
  **[default: 100000]**

- `--library-size-distribution` **[poisson|uniform|negative_binomial]**  
  Distribution for the library size.  
  **[default: uniform]**


- `--min-sparsity` **FLOAT**  
     Will archive the minimum specified sparcity by zeroing count randomly.  
  **[default: 0]**

- `--min-overlap` **INTEGER**  
  Minimum overlap for the blocks. Use this to evaluate overlap blocks, i.e. uninterrupted parts covered with reads that overlap on the genome.  
  **[default: 16]**

- `--group-orthology-level` **[very_low|low|normal|high|very_high]**  
  Determines the level of orthology in groups. If you use this, use it with a lot of threads. Takes a long time.  
  **[default: normal]**

- `--threads` **INTEGER**  
  Number of threads to be used.  
  **[default: 10]**
- `--version`  
  Show the version and exit.

- `--help`  
  Show this message and exit.

```

## Examples

### Running with Default Parameters

```sh
marbel
```

### Specifying Number of Species, Orthogroups, and Samples

```sh
marbel --n-species 10 --n-orthogroups 500 --n-samples 5 8
```

This command will generate a dataset with:

- 10 species
- 500 orthologous groups
- 5 samples for group 1
- 8 samples for group 2

## Release Guide

For publishing a new version to conda you can build it with the recipe. For this you need to have conda-build installed `(conda install conda-build`). For uploading to anaconda anaconda client needs to be installed. Installing the build tools and building in conda base environment is preferred.

Instructions (inside marbel directory):

```
cd recipe
conda-build .
```

Then upload to Anaconda, for this see:

https://www.anaconda.com/docs/tools/anaconda-org/user-guide/packages/conda-packages

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the [Apache License 2.0](LICENSE.md).

## Support

Feel free to reach out or open an issue if you have any questions or need further assistance with the usage of the tool.
