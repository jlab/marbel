import pymc as pm
import pandas as pd
from importlib import resources
from ete4 import Tree
from enum import Enum
import json

__version__ = "0.1.0"


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
MAX_ORTHO_GROUPS = 365813

data_package = str(resources.files(__package__) / 'data')

PATH_TO_GROUND_GENES = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz"
PATH_TO_GROUND_GENES_INDEX = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index"
PANGENOME_OVERVIEW = f"{data_package}/orthologues_processed_combined_all.parquet"
SPECIES_PHYLO_TREE = f"{data_package}/EDGAR_all_species.newick"
SPECIES_STATS = f"{data_package}/species_stats.json"
ALL_GTF_PATH = ""

pg_overview = pd.read_parquet(PANGENOME_OVERVIEW, engine='pyarrow')
species_tree = Tree(open(SPECIES_PHYLO_TREE))
DGE_LOG_2_CUTOFF_VALUE = 1
species_stats_dict = json.load(open(SPECIES_STATS, 'r'))

AVAILABLE_SPECIES = pg_overview.columns.to_list()[:MAX_SPECIES]

with pm.Model() as model:
    lognorm_dist = pm.Lognormal('read_counts', mu=READS_MEAN_LOG, sigma=READS_SD_LOG)  # reads
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


class OrthologyLevel(str, Enum):
    very_low = "very_low"
    low = "low"
    normal = "normal"
    high = "high"
    very_high = "very_high"


class SelectionCriterion(str, Enum):
    maximize = "maximize"
    minimize = "minimize"


class LibrarySizeDistribution():
    poisson = "poisson"
    uniform = "uniform"
    negative_binomial = "negative_binomial"

    possible_distributions = [poisson, uniform, negative_binomial]

    def __init__(self, distribution_name: str, poisson=50, nbin_n=20, nbin_p=0.3):
        if distribution_name not in self.possible_distributions:
            raise ValueError(f"Unknown distribution {distribution_name}")
        self.poisson = poisson
        self.nbin_n = nbin_n
        self.nbin_p = nbin_p
        self.distribution_name = distribution_name

    def __str__(self):
        match self.distribution_name:
            case "poisson":
                return f"{self.distribution_name} (lambda={self.poisson})"
            case "negative_binomial":
                return f"{self.distribution_name} (n={self.nbin_n}, p={self.nbin_p})"
            case "uniform":
                return f"{self.distribution_name}"
