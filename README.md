# meta_tran_sim_dev (Meta Transcriptomic Simulation)

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
meta-tran-sim --help
```

### Command Line Arguments

```
Usage: meta-tran-sim [OPTIONS] [N_SPECIES] [N_ORTHOGROUPS] [N_SAMPLES]...

Arguments:
  n_species         [N_SPECIES]      Number of species to be drawn for the metatranscriptomic in silico dataset [default: 20]
  n_orthogroups     [N_ORTHOGROUPS]  Number of orthologous groups to be drawn for the metatranscriptomic in silico dataset [default: 1000]
  n_samples          [N_SAMPLES]...    Number of samples to be created for the metatranscriptomic in silico dataset. The first number is the number of samples for group 1 and
                                     the second number is the number of samples for group 2 [default: 10, 10]
Options:
  --version                        Show the version and exit.
  --help                           Show this message and exit.
```

## Examples

### Running with Default Parameters

```sh
meta-tran-sim
```

### Specifying Number of Species, Orthogroups, and Samples

```sh
meta-tran-sim --n-species 30 --n-orthogroups 1500 --n-samples 15 20
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
