# marbel (MetAtranscriptomic Reference Builder Evaluation Library)

This project generates an in silico metatranscriptomic dataset based on specified parameters.

## Installation

### Conda build and install (recommended)

It is recomended to install the package with conda install.

Build the package with:

`conda build . `

For this you need to have conda-build installed `(conda install conda-build`)

Create new environment and install package:

```
conda create -n meta_tran_sim
conda activate meta_tran_sim
conda install --use-local meta_tran_sim
```

### Install by hand (for development purposes)

You need to install [R](https://www.r-project.org/about.html) and the R library polyester. Polyester can be installed with

```
R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("polyester")


```

Install the package:

```
pip install -e .
```

## Usage

To get help on how to use the script, run:

```sh
marbel --help
```

### Command Line Arguments

```
Usage: marbel [OPTIONS] 
Options:
 --n-species                 INTEGER                 Number of species to be drawn for the metatranscriptomic in silico dataset [default: 20]
 --n-orthogroups             INTEGER                 Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset [default: 1000]
 --n-samples                 <INTEGER INTEGER>...    Number of samples to be created for the metatranscriptomic in silico datasetthe first number is the number of samples for group 1 and the second is the number of samples for group 2 [default: 10, 10]  
  --outdir                   TEXT    Output directory for the metatranscriptomic in silico dataset [default: simulated_reads]
  --max-phylo-distance        TEXT    Maximum mean phylogenetic distance for orthologous groups. Specify stricter limit to avoid groups with a more diverse phylogenetic distance. [default: None]
  --min-identity              FLOAT   Minimum mean sequence identity score for orthologous groups. Specify for more stringent identity requirements. [default: None]
  --deg-ratio                <FLOAT FLOAT>... Ratio of up- and down-regulated genes. The first value is the ratio of up-regulated genes, the second represents the ratio of down-regulated genes [default: 0.1, 0.1]
  --seed                      INTEGER Seed for sampling. Set for reproducibility [default: None]
  --read-length               INTEGER Read length for the generated reads [default: 100]
  --output-format             [fastq.gz|fastq|fasta] Output format for the reads [default: fastq.gz]
  --version                           Show the version and exit.
  --help                              Show this message and exit.

```

## Examples

### Running with Default Parameters

```sh
marbel
```

### Specifying Number of Species, Orthogroups, and Samples

```sh
marbel --n-species 30 --n-orthogroups 1500 --n-samples 15 20
```

This command will generate a dataset with:

- 30 species
- 1500 orthologous groups
- 15 samples for group 1
- 20 samples for group 2

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the Apache License 2.0. See the [LICENSE]() file for details.

Feel free to reach out if you have any questions or need further assistance with the usage of the tool.
