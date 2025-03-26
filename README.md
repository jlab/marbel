[![Integration Tests & Code Style](https://github.com/jlab/marbel/actions/workflows/github_tests.yml/badge.svg)](https://github.com/jlab/marbel/actions/workflows/github_tests.yml) [![Coverage Status](https://coveralls.io/repos/github/jlab/marbel/badge.svg?branch=master)](https://coveralls.io/github/jlab/marbel?branch=master)

# Marbel



![Marbel logo](./marbel_logo.svg)Marbel (MetAtranscriptomic Reference Builder Evaluation Library) generates realistic *in silico* metatranscriptomic dataset based on specified parameters.

## Installation

### Install guide for development purposes

#### Install miniconda (if not installed already)

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
```

#### Create conda env

```
conda create -n marbel python=3.10 r-base
conda activate marbel
```

#### Install git-lfs (absolutely necessary)

Before cloning the repo you need to have git-lfs installed! If you do not have git-lfs and root rights install with

```
sudo apt-get install git-lfs
```

If you do not have root permission, install it in the Conda env:

```
conda install anaconda::git-lfs 
```

Now we need to initialize Git LFS:

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

### (Not ready, this is for later) nda build and install

It is recomended to install the package with conda install.

Build the package with:

`conda build . `

For this you need to have conda-build installed `(conda install conda-build`)

Create new environment and install package:

```
conda create -n marbel
conda activate marbel
conda install --use-local marbel
```

## Usage

To get help on how to use the script, run:

```sh
marbel --help
```

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
  **[default: 0.1]**

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

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the Apache License 2.0. See the [LICENSE]() file for details.

Feel free to reach out if you have any questions or need further assistance with the usage of the tool.
