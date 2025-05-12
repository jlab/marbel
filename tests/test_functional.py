import inspect
import sys
import pandas as pd
from os.path import join, basename
import skbio
import gzip
from os import path
import re
from Bio import SeqIO


def read_fastqGZ(fp_fastq):
    def _push(sequences, header, sequence, qualities):
        sequences.append({'header': header,
                          'sequence': sequence,
                          'qualities': qualities})
    OFFSET = 33
    sequences = []
    with gzip.open(fp_fastq, 'rb') as f:
        header = None
        sequence = ""
        qualities = []
        for i, l in enumerate(f.readlines()):
            line = l.decode("utf-8")
            if (i % 4 == 0) and (line.startswith('@')):
                if header is not None:
                    _push(sequences, header, sequence, qualities)
                header = line[1:].strip()
            elif (i % 4 == 1):
                sequence = line.strip()
            elif (i % 4 == 3):
                qualities = pd.Series(
                    map(lambda x: ord(x) - OFFSET, line.strip()))
        if header is not None:
            _push(sequences, header, sequence, qualities)
    return sequences


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
    params['Ratio of up and down regulated genes'] = float(
        params['Ratio of up and down regulated genes'])

    return params


def read_simulation_stats(fp_basedir):
    params = pd.read_csv(
        join(fp_basedir, 'summary', 'simulation_stats.txt'),
        index_col=0, sep=": ", header=None,
        names=['parameter', 'value'], engine='python')['value'].to_dict()

    params['Number of simulated genes'] = int(params['Number of simulated genes'])
    params['Number of simulated orthogroups'] = int(params['Number of simulated orthogroups'])
    params['Dropped genes due to all zero assignment by distribution'] = int(params['Dropped genes due to all zero assignment by distribution'])
    params['Dropped orthogroups due to all zero assignment by distribution'] = int(params['Dropped orthogroups due to all zero assignment by distribution'])

    return params


def test_gene_summary(genes, params):
    counts = genes[[s for s in genes.columns if 'sample' in s]]

    assert len(genes.index) == len(set(genes.index)), \
           "genes are ambiguous"

    assert len(sorted(genes["origin_species"].unique())) == \
           params['Number of species'], \
           "number of species does not match"

    genes['speciesID'] = list(map(lambda x: x.split('_')[0], genes.index))
    assert genes.groupby(['speciesID', 'origin_species']).size().shape[0] == \
           params['Number of species'], \
           "gene_name prefix does not match origin_species"

    assert genes['orthogroup'].nunique() <= params['Number of orthogroups'], \
        f"Too many orthogroups?! Found: {genes['orthogroup'].nunique()}, Expected max: {params['Number of orthogroups']}"

    assert genes.groupby('orthogroup').size().describe()['max'] > 1, \
           "orthogroups consists of only 1 CDS"

    # Check orthogroup multiplicity
    orthogroup_counts = (genes.groupby('orthogroup').size() > 1).value_counts()
    multi = orthogroup_counts.get(True, 0)
    single = orthogroup_counts.get(False, 0)
    assert multi > single * 0.2, "Surprisingly few orthogroups consist of more than one CDS"

    # Check count matrix sparsity
    sparsity_counts = (counts.stack() == 0).value_counts()
    zeros = sparsity_counts.get(True, 0)
    total = sparsity_counts.sum()
    assert zeros / total > 0.1, "Surprisingly low sparseness in count matrix!"

#    assert len(set(counts.sum())) == sum(params['Number of samples']), \
#           "not all samples have different coverage"
    assert counts.min().min() >= 0, \
           "negative count values?!"

    # foldchange <=0.5 (down) >=2 (up) (log2 foldchange<= -1 / >=1 )
    data = genes['Fold Change used for simulating counts, differs from the actual fold change'].apply(
        lambda x: 0.5 <= x <= 2).value_counts()
    if params['Number of orthogroups'] >= 1000:
        epsilon = (200 / params['Number of orthogroups'])
    else:
        epsilon = 0.4
    # print(params['Ratio of up and down regulated genes'])
    # print(data[False] / data.sum())
    # print(data)
    # print(params['Ratio of up and down regulated genes'] * (1 - epsilon) )
    # print(params['Ratio of up and down regulated genes'] * (1 + epsilon) )
    assert params['Ratio of up and down regulated genes'] * (1 - epsilon) <= \
           (data[False] / data.sum()) <= \
           params['Ratio of up and down regulated genes'] * (1 + epsilon), \
           "ratio of UP (or Down) regulated genes does not match parameters"

    print(len(set(list(map(lambda x: x.split('_')[1], genes.index)))),
          "number CDS?!")
    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


def test_simulation_stats(genes, params):
    assert genes['orthogroup'].unique().shape[0] == \
           params['Number of simulated orthogroups'], \
           "simulation stats and gene summ"
    assert genes.shape[0] == \
           params['Number of simulated genes'], \
           "too many orthogroups?!"


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
    if False in data:
        # TODO: assertion was changed to to removal of all zero genes
        assert data[False] < data[True] * 0.1, \
             "too few orthogroups with more than one member have identical" \
             " length sequences"
        #pass
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
    # due to the removal of the all zero row genes, this is not always true now; TODO: new assertion
    # assert data['max'] >= 0.99999, \
    #       "I would expect at least one orthogroup to be present in all " \
    #       "species, like a housekeeping gene"

    print("[OK] '%s' passed" % inspect.currentframe().f_code.co_name)


def test_fastq(fp_basedir, genes):
    # check fastq files
    fastqs = []
    for group in [1, 2]:
        for sample in range(1, params['Number of samples'][group - 1] + 1):
            for direction in [1, 2]:
                fastqs.append(
                    join(fp_basedir, 'group_%i_sample_%i_R%i.fastq.gz' % (
                        group, sample, direction)))

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


def test_number_of_reads(params):
    fp_basedir = "/home/tlin/marbel_testexecutions/test_run_mar25/simulated_reads"

    fastqs = []

    for group in [1, 2]:
        for sample in range(1, params['Number of samples'][group - 1] + 1):
            for direction in [1, 2]:
                fastqs.append(
                    join(fp_basedir, 'group_%i_sample_%i_R%i.fastq.gz' % (
                        group, sample, direction)))

    df_r1 = pd.DataFrame()
    df_r2 = pd.DataFrame()

    for fastq in fastqs:
        with gzip.open(fastq, 'rt') as handle:
            base_name = path.basename(fastq).replace(".fastq.gz", "")
            if base_name.endswith("_R1"):
                column_name = base_name.replace("_R1", "")
                ids = [re.sub(r'_\d+_\d+/\d$', '', record.id) for record in SeqIO.parse(handle, 'fastq')]
                gene_counts = pd.Series(ids).value_counts()
                gene_counts.name = column_name
                df_r1 = pd.concat([df_r1, gene_counts], axis=1)
            else:
                column_name = base_name.replace("_R2", "")
                ids = [re.sub(r'_\d+_\d+/\d$', '', record.id) for record in SeqIO.parse(handle, 'fastq')]
                gene_counts = pd.Series(ids).value_counts()
                gene_counts.name = column_name
                df_r2 = pd.concat([df_r2, gene_counts], axis=1)
    df_r1 = df_r1.fillna(0)
    df_r2 = df_r2.fillna(0)

    gene_summary = pd.read_csv(join(fp_basedir, 'summary', 'gene_summary.csv'))
    gene_summary_counts = gene_summary.loc[:, gene_summary.columns.str.contains("sample_")]
    gene_summary_counts_no_zero_lines = gene_summary_counts[~(gene_summary_counts == 0).all(axis=1)]
    assert df_r1.sort_index(axis=0).sort_index(axis=1).equals(df_r2.sort_index(axis=0).sort_index(axis=1)), "Simulated R1 and R2 mate reads are not equal"
    df_r1 = df_r1.astype(int).sort_index(axis=0).sort_index(axis=1)
    gene_summary_counts_no_zero_lines = gene_summary_counts_no_zero_lines.astype(int).sort_index(axis=0).sort_index(axis=1)
    assert 0 == (df_r1 - gene_summary_counts_no_zero_lines).sum().sum(), "Number of reads in fastq files does not match number of reads in gene_summary.csv"


if __name__ == "__main__":
    fp_basedir = sys.argv[1]

    # read used marbel parameters from file
    params = read_parameters(fp_basedir)
    sim_stats = read_simulation_stats(fp_basedir)

    # read gene summary information from file
    genes = pd.read_csv(join(fp_basedir, 'summary', 'gene_summary.csv'), index_col=0)

    # execute tests focussing on file gene_summary.csv
    test_gene_summary(genes, params)

    test_metaT_reference(fp_basedir, genes, params)

    test_tree(fp_basedir, genes)

    test_fastq(fp_basedir, genes)

    test_simulation_stats(genes, sim_stats)
