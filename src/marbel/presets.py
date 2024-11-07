import pymc as pm
import pandas as pd
from importlib import resources
from ete3 import Tree
from enum import Enum

__version__ = "0.1.1"


CPMS_MEAN_LOG = 3.73411080985053
CPMS_SD_LOG = 1.90214574373138
READS_MEAN_LOG = 5.57141333132782
READS_SD_LOG = 2.1412080853728
SPECIES_MEAN_LOG = 1.30069163364942
SPECIES_SD_LOG = 1.33546646662937
ORTHO_GROUPS_SHAPE = 0.675101585847872
ORTHO_GROUPS_RATE = 0.00603147150265822

DESEQ2_FITTED_A0 = 12.9105102
DESEQ2_FITTED_A1 = 0.3325853

DEFAULT_PHRED_QUALITY = 40

rank_distance = {
    "phylum": 0.6695424034631783,
    "class": 0.452175969905506,
    "order": 0.3307453184266344,
    "family": 0.21193467921976264,
    "genus": 0.11303678566362661,
}


MAX_SPECIES = 614
MAX_ORTHO_GROUPS = 485687

data_package = str(resources.files(__package__) / 'data')

PATH_TO_GROUND_GENES = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz"
PATH_TO_GROUND_GENES_INDEX = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index"
PANGENOME_OVERVIEW = f"{data_package}/orthologues_processed_combined_all.parquet"
SPECIES_PHYLO_TREE = f"{data_package}/EDGAR_all_species.newick"

pg_overview = pd.read_parquet(PANGENOME_OVERVIEW, engine='pyarrow')
species_tree = Tree(SPECIES_PHYLO_TREE)

AVAILABLE_SPECIES = pg_overview.columns.to_list()[:MAX_SPECIES]

with pm.Model() as model:
    lognorm_dist = pm.Lognormal('reads', mu=READS_MEAN_LOG, sigma=READS_SD_LOG)  # reads
    lognorm_dist2 = pm.Lognormal('species', mu=SPECIES_MEAN_LOG, sigma=SPECIES_SD_LOG)
    orthologues = pm.Gamma('ortho', alpha=ORTHO_GROUPS_SHAPE, beta=ORTHO_GROUPS_RATE)  # ortho groups


class ErrorModel(str, Enum):
    basic = "basic"
    perfect = "perfect"
    HiSeq = "HiSeq"
    NextSeq = "NextSeq"
    NovaSeq = "NovaSeq"
    Miseq_20 = "Miseq-20"
    Miseq_24 = "Miseq-24"
    Miseq_28 = "Miseq-28"
    Miseq_32 = "Miseq-32"



class Rank(str, Enum):
    phylum = "phylum"
    class_ = "class"
    order = "order"
    family = "family"
    genus = "genus"


class LibrarySizeDistribution(str, Enum):
    poisson = "poisson"
    uniform = "uniform"
    negative_binomial = "negative_binomial"