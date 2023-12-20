#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import argparse

def append_n(x):
    """
    e.g., pd.Series(['a', 'b', 'a']) --> pd.Series(['a (n=2)', 'b (n=1)', 'a (n=2)'])
    """
    if not isinstance(x, (pd.Series, np.ndarray, list)):
        raise ValueError('x must be an array/list/Series')
    tmp = list(x)
    vc = pd.Series(tmp).value_counts().to_dict()
    tmp = ['{} (n={:,})'.format(i, vc[i]) for i in tmp]
    if isinstance(x, pd.Series):
        tmp = pd.Series(tmp)
        tmp.index = x.index
    return tmp

parser = argparse.ArgumentParser()
parser.add_argument('--umap', required=True)
parser.add_argument('--clusters', required=True)
parser.add_argument('--atac-barcodes', required=True)
parser.add_argument('--rna-barcodes', required=True)
parser.add_argument('--doubletfinder', required=True, nargs='+', help='Doubletfinder output (one per library, filename {library}.*)')
parser.add_argument('--amulet', required=True, nargs='+', help='Amulet output (one per library, filename {library}.*)')
parser.add_argument('--demuxlet', required=False, nargs='+', default=None, help='Demuxlet output (one per library, filename {library}.*)')
parser.add_argument('--prefix', required=True)
args = parser.parse_args()

#UMAP = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/clean-and-cluster-drop-4-libraries-cellbender-no-grid-fpr-01/results/cluster/postdecontamination/contamination-filter/joint/no-sctransform/merged.umap.txt'
#CLUSTERS = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/clean-and-cluster-drop-4-libraries-cellbender-no-grid-fpr-01/results/cluster/postdecontamination/contamination-filter/joint/no-sctransform/merged.clusters.txt'
#UMAP = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/clean-and-cluster-drop-4-libraries-cellbender-no-grid-fpr-01-no-doublet-removal/work/30/1501023e144ba424162d34f8b12f42/merged.umap.txt'
#CLUSTERS = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/clean-and-cluster-drop-4-libraries-cellbender-no-grid-fpr-01-no-doublet-removal/work/30/1501023e144ba424162d34f8b12f42/merged.clusters.txt'
#DOUBLETFINDER = glob.glob('/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/clean-and-cluster-drop-4-libraries-cellbender-no-grid-fpr-01-no-doublet-removal-doubletfinder/results/doubletfinder/*')
#AMULET = glob.glob('/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/work/qc-no-dropkick/results/atac-doublet-detection/*.doublet_probabilities.txt')
#ATAC_BARCODES = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/pipelines/snATACseq-NextFlow/737K-arc-v1.txt'
#RNA_BARCODES = '/scratch/scjp_root/scjp0/porchard/pfizer-mouse-multiome/pipelines/snRNAseq-NextFlow/737K-arc-v1.txt'
#PREFIX = ''

UMAP = args.umap
CLUSTERS = args.clusters
DOUBLETFINDER = args.doubletfinder
DEMUXLET = args.demuxlet
AMULET = args.amulet
PREFIX = args.prefix
ATAC_BARCODES = args.atac_barcodes
RNA_BARCODES = args.rna_barcodes

atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0].to_list()
rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0].to_list()
atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))

doubletfinder = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['barcode', 'doubletfinder_assignment']).assign(library=os.path.basename(f).split('.')[0]) for f in DOUBLETFINDER])
doubletfinder['nucleus'] = doubletfinder.library + '-' + doubletfinder.barcode

if DEMUXLET is not None:
    demuxlet = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['barcode', 'demuxlet_assignment'], usecols=[0, 1]).assign(library=os.path.basename(f).split('.')[0]) for f in DEMUXLET])
    demuxlet['nucleus'] = demuxlet.library + '-' + demuxlet.barcode

amulet = pd.concat([pd.read_csv(f, sep='\t').assign(library=os.path.basename(f).split('.')[0]) for f in AMULET])
amulet.barcode = amulet.barcode.map(atac_to_rna_barcode)
amulet['nucleus'] = amulet.library + '-' + amulet.barcode
amulet['amulet_assignment'] = np.where(amulet['q-value']<=0.05, 'doublet', 'singlet')

umap = pd.read_csv(UMAP, sep='\t', header=None, names=['nucleus', 'dim1', 'dim2'])
umap.nucleus = umap.nucleus.str.replace('___', '-')

clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['nucleus', 'cluster'])
clusters.nucleus = clusters.nucleus.str.replace('___', '-')

assert(all(clusters.nucleus.isin(umap.nucleus)))
assert(all(umap.nucleus.isin(doubletfinder.nucleus)))
assert(all(umap.nucleus.isin(amulet.nucleus)))

if DEMUXLET is not None:
    assert(all(umap.nucleus.isin(demuxlet.nucleus)))

umap = umap.merge(doubletfinder[['nucleus', 'doubletfinder_assignment']])
if DEMUXLET is not None:
    umap = umap.merge(demuxlet[['nucleus', 'demuxlet_assignment']])
umap = umap.merge(amulet[['nucleus', 'amulet_assignment']])
umap = umap.merge(clusters)

umap['doubletfinder_is_doublet'] = (umap.doubletfinder_assignment == 'doublet').astype(int)
if DEMUXLET is not None:
    umap['demuxlet_is_doublet'] = (umap.demuxlet_assignment != 'SNG').astype(int)
umap['amulet_is_doublet'] = (umap.amulet_assignment == 'doublet').astype(int)
if DEMUXLET is not None:
    umap['doublet_by_any_method'] = (umap[['doubletfinder_is_doublet', 'demuxlet_is_doublet', 'amulet_is_doublet']].sum(axis=1) > 0).astype(int)
else:
    umap['doublet_by_any_method'] = (umap[['doubletfinder_is_doublet', 'amulet_is_doublet']].sum(axis=1) > 0).astype(int)

cluster_fracs = umap.groupby('cluster').doublet_by_any_method.mean()
doublet_clusters = cluster_fracs[cluster_fracs>=0.5].index.to_list()
umap['keep'] = (~umap.cluster.isin(doublet_clusters)) & (umap.doublet_by_any_method==0)

def make_label(doubletfinder_is_doublet, amulet_is_doublet, demuxlet_is_doublet=None):
    x = [doubletfinder_is_doublet, amulet_is_doublet] if demuxlet_is_doublet is None else [doubletfinder_is_doublet, amulet_is_doublet, demuxlet_is_doublet]
    y = ['DF', 'AMULET'] if demuxlet_is_doublet is None else ['DF', 'AMULET', 'Demuxlet']
    if any(x):
        return 'Doublet ({})'.format('+'.join([software for is_doublet, software in zip(x, y) if is_doublet]))
    else:
        return 'Singlet'

if DEMUXLET is not None:
    umap['label'] = [make_label(i, j, k) for i, j, k in zip(umap.doubletfinder_is_doublet.astype(bool), umap.amulet_is_doublet.astype(bool), umap.demuxlet_is_doublet.astype(bool))]
else:
    umap['label'] = [make_label(i, j) for i, j in zip(umap.doubletfinder_is_doublet.astype(bool), umap.amulet_is_doublet.astype(bool))]
umap.label = append_n(umap.label)


if DEMUXLET is not None:
    doublet_rates = umap[['cluster', 'doubletfinder_is_doublet', 'demuxlet_is_doublet', 'amulet_is_doublet', 'doublet_by_any_method']].groupby('cluster').mean().reset_index()
else:
    doublet_rates = umap[['cluster', 'doubletfinder_is_doublet', 'amulet_is_doublet', 'doublet_by_any_method']].groupby('cluster').mean().reset_index()

fig, ax = plt.subplots()
sns.scatterplot(x='cluster', y='doubletfinder_is_doublet', ax=ax, data=doublet_rates, label='Doublet rate (doubletfinder)')
sns.scatterplot(x='cluster', y='amulet_is_doublet', ax=ax, data=doublet_rates, label='Doublet rate (amulet)')
if DEMUXLET is not None:
    sns.scatterplot(x='cluster', y='demuxlet_is_doublet', ax=ax, data=doublet_rates, label='Doublet rate (demuxlet)')
sns.scatterplot(x='cluster', y='doublet_by_any_method', ax=ax, data=doublet_rates, label='Doublet rate (any)')
ax.set_xticks(range(0, doublet_rates.cluster.max()+1))
ax.set_ylabel('Doublet rate')
ax.set_xlabel('Cluster')
ax.grid(True)
ax.axhline(0.5, color='black', ls='--')
fig.tight_layout()
fig.savefig(f'{PREFIX}doublet-rates-per-cluster.png', dpi=300, facecolor='white')


fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x='dim1', y='dim2', hue='label', data=umap.sample(frac=1.0), alpha=0.5, s=5)
ax.legend().set_title('')
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('UMAP dim. 2')
fig.tight_layout()
fig.savefig(f'{PREFIX}doublets-on-umap.png', dpi=300, facecolor='white')


umap['label'] = np.where(umap.keep, 'Keep', 'Drop')
umap.label = append_n(umap.label)

fig, ax = plt.subplots(figsize=(10, 10))

sns.scatterplot(x='dim1', y='dim2', hue='label', data=umap.sample(frac=1.0), alpha=0.5, s=5)
ax.legend().set_title('')
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('UMAP dim. 2')
fig.tight_layout()
fig.savefig(f'{PREFIX}keep-on-umap.png', dpi=300, facecolor='white')


# keep_nuclei = umap[umap.keep].nucleus.to_list()
# with open(f'{PREFIX}keep.txt', 'w') as fh:
#     for n in keep_nuclei:
#         fh.write(n + '\n')