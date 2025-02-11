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
parser.add_argument('--atac-barcodes', required=False)
parser.add_argument('--rna-barcodes', required=False)
parser.add_argument('--doubletfinder', required=True, nargs='+', default=None, help='Doubletfinder output (one per library, filename {library}.*)')
parser.add_argument('--amulet', required=False, nargs='+', default=None, help='Amulet output (one per library, filename {library}.*)')
parser.add_argument('--demuxlet', required=False, nargs='+', default=None, help='Demuxlet output (one per library, filename {library}.*)')
parser.add_argument('--prefix', required=True)
args = parser.parse_args()

UMAP = args.umap
CLUSTERS = args.clusters
DOUBLETFINDER = args.doubletfinder
DEMUXLET = args.demuxlet
AMULET = args.amulet
PREFIX = args.prefix
ATAC_BARCODES = args.atac_barcodes
RNA_BARCODES = args.rna_barcodes

umap = pd.read_csv(UMAP, sep='\t', header=None, names=['nucleus', 'dim1', 'dim2'])
umap.nucleus = umap.nucleus.str.replace('___', '-')

clusters = pd.read_csv(CLUSTERS, sep='\t', header=None, names=['nucleus', 'cluster'])
clusters.nucleus = clusters.nucleus.str.replace('___', '-')

assert(all(clusters.nucleus.isin(umap.nucleus)))


methods = []

if DOUBLETFINDER is not None:
    methods.append('doubletfinder')
    doubletfinder = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['barcode', 'doubletfinder_assignment']).assign(library=os.path.basename(f).split('.')[0]) for f in DOUBLETFINDER])
    doubletfinder['nucleus'] = doubletfinder.library + '-' + doubletfinder.barcode
    assert(all(umap.nucleus.isin(doubletfinder.nucleus)))

if DEMUXLET is not None:
    methods.append('demuxlet')
    demuxlet = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['barcode', 'demuxlet_assignment'], usecols=[0, 1]).assign(library=os.path.basename(f).split('.')[0]) for f in DEMUXLET])
    demuxlet['nucleus'] = demuxlet.library + '-' + demuxlet.barcode
    assert(all(umap.nucleus.isin(demuxlet.nucleus)))

if AMULET is not None:
    methods.append('amulet')
    atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0].to_list()
    rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0].to_list()
    atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))

    amulet = pd.concat([pd.read_csv(f, sep='\t').assign(library=os.path.basename(f).split('.')[0]) for f in AMULET])
    amulet.barcode = amulet.barcode.map(atac_to_rna_barcode)
    amulet['nucleus'] = amulet.library + '-' + amulet.barcode
    amulet['amulet_assignment'] = np.where(amulet['q-value']<=0.05, 'doublet', 'singlet')
    assert(all(umap.nucleus.isin(amulet.nucleus)))


if DOUBLETFINDER is not None:
    umap = umap.merge(doubletfinder[['nucleus', 'doubletfinder_assignment']])
if DEMUXLET is not None:
    umap = umap.merge(demuxlet[['nucleus', 'demuxlet_assignment']])
if AMULET is not None:
    umap = umap.merge(amulet[['nucleus', 'amulet_assignment']])
umap = umap.merge(clusters)


if DOUBLETFINDER is not None:
    umap['doubletfinder_is_doublet'] = (umap.doubletfinder_assignment == 'doublet').astype(int)
if DEMUXLET is not None:
    umap['demuxlet_is_doublet'] = (umap.demuxlet_assignment != 'SNG').astype(int)
if AMULET is not None:
    umap['amulet_is_doublet'] = (umap.amulet_assignment == 'doublet').astype(int)


umap['doublet_by_any_method'] = (umap[[f'{i}_is_doublet' for i in methods]].sum(axis=1) > 0).astype(int)

cluster_fracs = umap.groupby('cluster').doublet_by_any_method.mean()
# doublet_clusters = cluster_fracs[cluster_fracs>=0.5].index.to_list()
# umap['keep'] = (~umap.cluster.isin(doublet_clusters)) & (umap.doublet_by_any_method==0)

def make_label(x):
    software = [i.replace('_is_doublet', '') for i in x.index]
    is_doublet = x.to_list()
    if any(x):
        return 'Doublet ({})'.format('+'.join([s for i, s in zip(is_doublet, software) if i]))
    else:
        return 'Singlet'

umap['label'] = umap[[f'{i}_is_doublet' for i in methods]].astype(bool).apply(make_label, axis=1)
umap.label = append_n(umap.label)


doublet_rates = umap[['cluster', 'doublet_by_any_method'] + [f'{i}_is_doublet' for i in methods]].groupby('cluster').mean().reset_index()


fig, ax = plt.subplots()
for method in methods:
    sns.scatterplot(x='cluster', y=f'{method}_is_doublet', ax=ax, data=doublet_rates, label=f'Doublet rate ({method})')
if len(methods) > 1:
    sns.scatterplot(x='cluster', y='doublet_by_any_method', ax=ax, data=doublet_rates, label='Doublet rate (any)')
ax.set_xticks(range(0, doublet_rates.cluster.max()+1))
ax.set_ylabel('Doublet rate')
ax.set_xlabel('Cluster')
ax.grid(True)
# ax.axhline(0.5, color='black', ls='--')
fig.tight_layout()
fig.savefig(f'{PREFIX}doublet-rates-per-cluster.png', dpi=300, facecolor='white')


fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x='dim1', y='dim2', hue='label', data=umap.sample(frac=1.0), alpha=0.5, s=5)
ax.legend().set_title('')
ax.set_xlabel('UMAP dim. 1')
ax.set_ylabel('UMAP dim. 2')
fig.tight_layout()
fig.savefig(f'{PREFIX}doublets-on-umap.png', dpi=300, facecolor='white')


# umap['label'] = np.where(umap.keep, 'Keep', 'Drop')
# umap.label = append_n(umap.label)

# fig, ax = plt.subplots(figsize=(10, 10))

# sns.scatterplot(x='dim1', y='dim2', hue='label', data=umap.sample(frac=1.0), alpha=0.5, s=5)
# ax.legend().set_title('')
# ax.set_xlabel('UMAP dim. 1')
# ax.set_ylabel('UMAP dim. 2')
# fig.tight_layout()
# fig.savefig(f'{PREFIX}keep-on-umap.png', dpi=300, facecolor='white')


# keep_nuclei = umap[umap.keep].nucleus.to_list()
# with open(f'{PREFIX}keep.txt', 'w') as fh:
#     for n in keep_nuclei:
#         fh.write(n + '\n')