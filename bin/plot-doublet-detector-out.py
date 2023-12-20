#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

DOUBLET_PROBABILITIES, ATAQV_METRICS, RNA_BARCODES, ATAC_BARCODES, LIBRARY = sys.argv[1:]

# read in barcodes
rna_barcodes = pd.read_csv(RNA_BARCODES, header=None)[0].to_list()
atac_barcodes = pd.read_csv(ATAC_BARCODES, header=None)[0].to_list()
if not len(rna_barcodes) == len(atac_barcodes):
    raise ValueError('Number of items in RNA barcode list != number of items in ATAC barcode list')
rna_to_atac_barcode = dict(zip(rna_barcodes, atac_barcodes))
atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))


dp = pd.read_csv(DOUBLET_PROBABILITIES, sep='\t')
dp['is_doublet'] = dp['q-value'] <= 0.01
barcode_to_assignment = dict(zip(dp.barcode, dp.is_doublet))

singlets_atac = dp[~dp.is_doublet].barcode.to_list()
singlets_rna = [atac_to_rna_barcode[i] for i in singlets_atac]
with open(f'{LIBRARY}.atac.singlets.txt', 'w') as fh:
    fh.write('\n'.join(singlets_atac) + '\n')
with open(f'{LIBRARY}.rna.singlets.txt', 'w') as fh:
    fh.write('\n'.join(singlets_rna) + '\n')

# plot the q-value distribution and note the number of singlets/doublets and doublet fraction
fig, ax = plt.subplots()
ax = dp['q-value'].hist(bins=100, ax=ax)
ax.set_xlabel('q-value (nucleus is doublet)')
ax.set_ylabel('Nuclei')
NUMBER_DOUBLETS = dp.is_doublet.sum()
NUMBER_NONDOUBLETS = (~dp.is_doublet).sum()
PERCENT_DOUBLETS = round(100*NUMBER_DOUBLETS/(NUMBER_DOUBLETS+NUMBER_NONDOUBLETS), 2)
ax.set_title(f'{NUMBER_NONDOUBLETS} singlets, {NUMBER_DOUBLETS} doublets\n(doublet percentage = {PERCENT_DOUBLETS}% @ 1% FDR)')
fig.tight_layout()
fig.savefig(f'{LIBRARY}.doublet-summary.png')
fig.clf()


ataqv = pd.read_csv(ATAQV_METRICS, sep='\t', header=None, names=['barcode', 'metric', 'value'])
ataqv = ataqv[ataqv.barcode.isin(dp['barcode'].to_list())].pivot(index='barcode', columns='metric', values='value')
for metric in ['hqaa', 'total_autosomal_reads', 'total_reads']:
    ataqv[metric] = ataqv[metric].astype(int)
for metric in ['max_fraction_reads_from_single_autosome', 'tss_enrichment'] + [i for i in ataqv.columns if 'percent' in i]:
    ataqv[metric] = ataqv[metric].astype(float)
ataqv['assignment'] = ataqv.index.map(lambda x: 'singlet' if not barcode_to_assignment[x] else 'doublet')


# plot HQAA of doublets vs non-doublets
# plot max_fraction_reads_from_single_autosome for doublets
fig, ax = plt.subplots()

sns.stripplot(x='assignment', y='hqaa', data=ataqv, ax=ax, alpha=0.3)
ax.set_yscale('log')
ax.set_ylim(bottom=1000)
ax.set_xlabel('')
ax.set_ylabel('HQAA')
fig.tight_layout()
fig.savefig(f'{LIBRARY}.doublet-vs-singlet-hqaa.png')
fig.clf()
