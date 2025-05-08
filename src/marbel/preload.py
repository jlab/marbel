import pymc as pm
from marbel.presets import (
    PANGENOME_OVERVIEW, SPECIES_PHYLO_TREE, SPECIES_STATS,
    READS_MEAN_LOG, READS_SD_LOG, SPECIES_MEAN_LOG, SPECIES_SD_LOG,
    ORTHO_GROUPS_SHAPE, ORTHO_GROUPS_RATE
)
from ete4 import Tree
import json
import pandas as pd
from functools import lru_cache


@lru_cache(maxsize=1)
def get_pg_overview():
    return pd.read_parquet(PANGENOME_OVERVIEW, engine='pyarrow')


@lru_cache(maxsize=1)
def get_species_tree():
    with open(SPECIES_PHYLO_TREE) as f:
        return Tree(f.read())


@lru_cache(maxsize=1)
def get_species_stats_dict():
    with open(SPECIES_STATS, 'r') as f:
        return json.load(f)


@lru_cache(maxsize=1)
def get_pymc_model():
    with pm.Model() as model:
        pm.Lognormal('read_counts', mu=READS_MEAN_LOG, sigma=READS_SD_LOG)
        pm.Lognormal('species', mu=SPECIES_MEAN_LOG, sigma=SPECIES_SD_LOG)
        pm.Gamma('ortho', alpha=ORTHO_GROUPS_SHAPE, beta=ORTHO_GROUPS_RATE)
    return model
