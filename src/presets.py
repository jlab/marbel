import pymc as pm
import pandas as pd

__version__ = "0.1.1"


CPMS_MEAN_LOG = 3.73411080985053
CPMS_SD_LOG = 1.90214574373138
READS_MEAN_LOG = 5.57141333132782
READS_SD_LOG = 2.1412080853728
SPECIES_MEAN_LOG = 1.30069163364942
SPECIES_SD_LOG = 1.33546646662937
ORTHO_GROUPS_SHAPE = 0.675101585847872
ORTHO_GROUPS_RATE = 0.00603147150265822

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

PATH_TO_GROUND_GENES = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/edgar/edgar_big_run/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz"
PATH_TO_GROUND_GENES_INDEX = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/edgar/edgar_big_run/deduplicated_pangenome_EDGAR_Microbiome_JLAB2.fas.bgz.bio_index"
PANGENOME_OVERVIEW = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/edgar/edgar_big_run/pangenome_EDGAR_Microbiome_JLAB2.csv"
PANGENOME_OVERVIEW = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/calculations/orthologues_df.csv"
PANGENOME_OVERVIEW = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/calculations/chunks/orthologues_processed_combined_all.csv"
PANGENOME_OVERVIEW = "/mnt/tlin/vol/jlab/tlin/all_project/in_silico_dataset/calculations/orthologues_processed_combined_all.parquet"

pg_overview = pd.read_parquet(PANGENOME_OVERVIEW, engine='pyarrow')

AVAILABLE_SPECIES = pg_overview.columns.to_list()[:MAX_SPECIES]

with pm.Model() as model:
    lognorm_dist = pm.Lognormal('reads', mu=READS_MEAN_LOG, sigma=READS_SD_LOG) #reads
    lognorm_dist2 = pm.Lognormal('species', mu=SPECIES_MEAN_LOG, sigma=SPECIES_SD_LOG) 
    orthologues = pm.Gamma('ortho', alpha=ORTHO_GROUPS_SHAPE, beta=ORTHO_GROUPS_RATE) #ortho groups