from importlib import resources
from enum import Enum

__version__ = "0.2.0"

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
DGE_LOG_2_CUTOFF_VALUE = 1
MAX_SPECIES = 614
MAX_ORTHO_GROUPS = 365813

rank_distance = {
    "phylum": 0.6695424034631783,
    "class": 0.452175969905506,
    "order": 0.3307453184266344,
    "family": 0.21193467921976264,
    "genus": 0.11303678566362661,
}


data_package = str(resources.files(__package__) / 'data')

PATH_TO_GROUND_GENES = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz"
PATH_TO_GROUND_GENES_INDEX = f"{data_package}/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index"
PANGENOME_OVERVIEW = f"{data_package}/orthologues_processed_combined_all.parquet"
CDS_GENOMIC_LOCATIONS = f"{data_package}/cds_genomic_locations.parquet"
SPECIES_PHYLO_TREE = f"{data_package}/EDGAR_all_species.newick"
SPECIES_STATS = f"{data_package}/species_stats.json"


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
