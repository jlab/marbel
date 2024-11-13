import inspect
import sys
import pandas as pd
from os.path import join, basename
import skbio


def read_parameters(fp_basedir):
    params = pd.read_csv(
        join(fp_basedir, 'summary', 'marbel_params.txt'),
        index_col=0, sep=": ", header=None,
        names=['parameter', 'value'], engine='python')['value'].to_dict()

    # fix data types
    params['Number of species'] = int(params['Number of species'])
    params['Number of orthogroups'] = int(params['Number of orthogroups'])
    params['Number of samples'] = list(
        map(int, params['Number of samples']
            .replace('(', '')
            .replace(')', '')
            .replace(' ', '').split(',')))
    params['Up and down regulated genes'] = float(
        params['Up and down regulated genes'])

    return params


def test_gene_summary(genes, params):
    counts = genes[[s for s in genes.columns if s.startswith('sample_')]]

    assert len(genes.index) == len(set(genes.index)), \
           "genes are ambiguous"

    assert len(sorted(genes["origin_species"].unique())) == \
           params['Number of species'], \
           "number of species does not match"

    genes['speciesID'] = list(map(lambda x: x.split('_')[0], genes.index))
    assert genes.groupby(['speciesID', 'origin_species']).size().shape[0] == \
           params['Number of species'], \
           "gene_name prefix does not match origin_species"

    genes['cdsID'] = list(map(lambda x: x.split('_')[-1], genes.index))
    assert genes['orthogroup'].unique().shape[0] == \
           params['Number of orthogroups'], \
           "number orthogroups does not match"

    assert genes.groupby('orthogroup').size().describe()['max'] > 1, \
           "orthogroups consists of only 1 CDS"
    data = (genes.groupby('orthogroup').size() > 1).value_counts()
    assert data[True] > data[False] * 0.2, \
           "surprisingly few orthogroups consists of more than one CDS"

    data = (counts.stack() == 0).value_counts()
    assert data[True] / data.sum() > 0.1, \
           "surprisingly low sparseness in count matrix!"

    assert len(set(counts.sum())) == sum(params['Number of samples']), \
           "not all samples have different coverage"
    assert counts.min().min() >= 0, \
           "negative count values?!"

    # foldchange <=0.5 (down) >=2 (up) (log2 foldchange<= -1 / >=1 )
    data = genes['fold_change_ratio'].apply(
        lambda x: 0.5 <= x <= 2).value_counts()
    epsilon = 0.02
    assert params['Up and down regulated genes'] * 2 * (1 - epsilon) <= \
           (data[False] / data.sum()) <= \
           params['Up and down regulated genes'] * 2 * (1 + epsilon), \
           "ratio of UP (or Down) regulated genes does not match parameters"

    print(len(set(list(map(lambda x: x.split('_')[1], genes.index)))),
          "number CDS?!")
    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


def test_metaT_reference(fp_basedir, genes, params):
    seqs = dict()
    for seq in skbio.io.read(join(fp_basedir, 'summary',
                             'metatranscriptome_reference.fasta'),
                             format='fasta'):
        seqs[seq.metadata['id']] = seq
    assert len(set(seqs.keys()) & set(genes.index)) == genes.shape[0], \
           "mismatch in gene identifiers between fasta file and summary table!"

    genes['gene_length'] = list(map(lambda x: len(seqs[x]), genes.index))
    seqLens = pd.Series([len(seq) for seq in seqs.values()], name='gene_length')
    data = (seqLens >= 300).value_counts()
    assert data[True] > data[False], \
           "too few gene lengths are shorter than single Illumina " \
           "reads of 300bp!"

    # for most orthogroups with more than one member, I expect to have
    # different gene sequences
    data = (pd.Series([g['gene_length'].unique().shape[0]
            for og, g in
            genes.groupby('orthogroup')
            if g.shape[0] > 1]) > 1).value_counts()
    assert data[False] < data[True] * 0.1, \
           "too few orthogroups with more than one member have identical" \
           " length sequences"

    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


def test_tree(fp_basedir, genes):
    tree = skbio.tree.TreeNode.read(
        join(fp_basedir, 'summary', 'species_tree.newick'), format='newick')
    og_counts = pd.pivot_table(
        data=genes, index='origin_species', columns='orthogroup',
        values='read_mean_count').applymap(lambda x: 1 if pd.notnull(x) else 0)

    # TreeNode.read converts underscore to space?!
    og_counts.index = list(map(lambda x: x.replace('_', ' '), og_counts.index))
    og_faith = [skbio.diversity.alpha.faith_pd(
                og_counts.iloc[:, i], otu_ids=og_counts.index, tree=tree,
                validate=True)
                for i in range(og_counts.shape[1])]
    data = (pd.Series(og_faith) / tree.descending_branch_length()).describe()

    assert tree.descending_branch_length() > 0, \
           "species tree seems to have no branch length"
    assert data['min'] > 0, \
           "there seems to be an orthogroup with no member in any species?!"
    assert data['max'] >= 0.99999, \
           "I would expect at least one orthogroup to be present in all " \
           "species, like a housekeeping gene"

    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


def test_fastq(fp_basedir, genes):
    # check fastq files
    fastqs = []
    for group in [1, 2]:
        for sample in range(1, params['Number of samples'][group - 1] + 1):
            for direction in [1, 2]:
                fastqs.append(
                    join(fp_basedir, 'sample_%i_group%i_R%i.fastq.gz' % (
                        sample, group, direction)))

    phredscores = dict()
    print('Loading sequencing data: ', end="", file=sys.stderr)
    for fq in fastqs:
        print('.', end="", file=sys.stderr)
        seqs = skbio.io.read(fq, format='fastq', variant='illumina1.8')
        quals = []
        for i, seq in enumerate(seqs):
            x = seq.positional_metadata['quality']
            x.name = seq.metadata['id']
            quals.append(x)
            if i > 1000:
                break
        quals = pd.concat(quals, axis=1)
        phredscores[basename(fq)] = quals
    print(' done.', file=sys.stderr)

    for fp, quals in phredscores.items():
        data = quals.mean(axis=1)
        assert data.iloc[10:20].mean() > data.iloc[-10:].mean(), \
            "surprisingly, read starts have higher phred scores than ends for" \
            " sample %s" % fp
        if '_R1.' in fp:
            dataR2 = phredscores[fp.replace('_R1.', '_R2.')].mean(axis=1)
            assert data.mean() > dataR2.mean(), \
                   "R2 phreds should be worse than R2 for sample %s" % fp

    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


if __name__ == "__main__":
    fp_basedir = sys.argv[1]

    # read used marbel parameters from file
    params = read_parameters(fp_basedir)

    # read gene summary information from file
    genes = pd.read_csv(join(fp_basedir, 'summary', 'gene_summary.csv'),
                        index_col=0)

    # execute tests focussing on file gene_summary.csv
    test_gene_summary(genes, params)

    test_metaT_reference(fp_basedir, genes, params)

    test_tree(fp_basedir, genes)

    test_fastq(fp_basedir, genes)