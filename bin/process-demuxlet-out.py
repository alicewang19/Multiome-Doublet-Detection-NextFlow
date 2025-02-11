#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--atac-demuxlet', required=False, help='Demuxlet best file')
parser.add_argument('--rna-demuxlet', required=False, help='Demuxlet best file')
parser.add_argument('--atac-barcodes', required=False, help='List of whitelisted ATAC barcodes')
parser.add_argument('--rna-barcodes', required=False, help='List of whitelisted RNA barcodes')
parser.add_argument('--strategy', help='RNA, ATAC, joint (use either RNA assignment, ATAC assignment, or joint assignment)')
parser.add_argument('--prefix', default='demuxlet.', help='')
args = parser.parse_args()


def load_demuxlet(f):
    tmp = pd.read_csv(f, sep='\t')
    return tmp

def recode_best(x):
    tmp = x.split('-')
    category = tmp[0]
    samples = '-'.join(tmp[1:])
    return f'{category} ({samples})' if category == 'SNG' else category


@ticker.FuncFormatter
def read_count_formatter(x, pos):
    """
    Tick label formatting function that converts labels to B/M/k (billions, millions, thousands).

    Usage:
    ax.xaxis.set_major_formatter(read_count_formatter)
    """
    if abs(x) >= 1e9:
        return '{}B'.format(x/1e9)
    elif abs(x) >= 1e6:
        return '{}M'.format(x/1e6)
    elif abs(x) >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


def make_colormap_dict(keys, palette='viridis'):
    """
    Given list of items (in order), create a dict of item --> color.

    Input:
    keys: list of items.
    palette: name of matplotlib color palette to use

    Returns: Dict of item --> color (hex)
    """
    assert(isinstance(keys, list))
    assert(isinstance(palette, str))
    cmap = mpl.cm.get_cmap(palette, len(keys))
    return {keys[i]: cmap(i) for i in range(cmap.N)}


def fix_heatmap_limits(ax):
    """
    Used to fix e.g. this bug: https://github.com/matplotlib/matplotlib/issues/14675
    """
    bottom, top = ax.get_ylim()
    if bottom > top:
        bottom += 0.5
        top -= 0.5
    else:
        bottom -= 0.5
        top += 0.5
    ax.set_ylim(bottom, top)
    return True


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


if args.atac_demuxlet is None and args.rna_demuxlet is None:
    raise ValueError('Must provide at least one of --atac-demuxlet or --rna-demuxlet file')
elif args.atac_demuxlet is not None and args.rna_demuxlet is not None:
    rna_barcodes = pd.read_csv(args.rna_barcodes, header=None)[0].to_list()
    atac_barcodes = pd.read_csv(args.atac_barcodes, header=None)[0].to_list()
    if not len(rna_barcodes) == len(atac_barcodes):
        raise ValueError('Number of items in RNA barcode list != number of items in ATAC barcode list')
    atac_to_rna_barcode = dict(zip(atac_barcodes, rna_barcodes))


    atac_demuxlet = load_demuxlet(args.atac_demuxlet).assign(modality='ATAC')
    assert(all(atac_demuxlet.BARCODE.isin(set(atac_barcodes))))
    atac_demuxlet['RNA_BARCODE'] = atac_demuxlet.BARCODE.map(atac_to_rna_barcode)
    atac_demuxlet['assignment'] = atac_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)

    rna_demuxlet = load_demuxlet(args.rna_demuxlet).assign(modality='RNA')
    rna_demuxlet['assignment'] = rna_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)
    rna_demuxlet['RNA_BARCODE'] = rna_demuxlet.BARCODE

    demuxlet = atac_demuxlet[['RNA_BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'ATAC', 'NUM.SNPS': 'ATAC_SNPs'}).merge(rna_demuxlet[['RNA_BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'RNA', 'NUM.SNPS': 'RNA_SNPs'}))

    demuxlet['label'] = np.where(demuxlet.RNA ==  demuxlet.ATAC, demuxlet.RNA.map(lambda x: f'both = {x}'), demuxlet[['RNA', 'ATAC']].apply(lambda x: f'RNA = {x[0]}, ATAC={x[1]}', axis=1))
    # if RNA and ATAC both say singlet from same individual, assign it as such
    # otherwise, assign as doublet
    demuxlet['joint'] = np.where(demuxlet.RNA == demuxlet.ATAC, demuxlet.RNA, 'DBL')
    demuxlet['final_assignment'] = demuxlet[args.strategy].str.split(' ').map(lambda x: x[0])
    demuxlet[['RNA_BARCODE', 'final_assignment', args.strategy]].to_csv(f'{args.prefix}assignments.txt', sep='\t', index=False, header=False)


    cmap = make_colormap_dict(list(sorted(set(demuxlet.ATAC.to_list() + demuxlet.RNA.to_list()))))

    summarize = demuxlet.groupby(['ATAC', 'RNA']).size().rename('n').reset_index().pivot(index='ATAC', columns='RNA', values='n').fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(1+len(summarize.index), 1+len(summarize.columns)))
    sns.heatmap(summarize, annot=True, ax=ax)
    fix_heatmap_limits(ax)
    fig.savefig(f'{args.prefix}demuxlet-heatmap.png', bbox_inches='tight', dpi=300)
    fig.clf()



    fig, axs = plt.subplots(ncols=3, figsize=(5*3, 5))

    ax = axs[0]
    ax.set_title('RNA assignment')
    rna_counts = demuxlet.RNA.value_counts().rename('RNA').rename_axis(index='index').reset_index()
    rna_counts['percentage'] = (100*rna_counts.RNA/rna_counts.RNA.sum()).map(lambda x: round(x, 2))
    rna_counts['label'] = rna_counts['index'] + ' - ' + rna_counts.percentage.astype(str) + '%'
    sns.barplot(x='RNA', y='label', ax=ax, data=rna_counts, palette=cmap, hue='index', dodge=False)
    ax.legend().remove()
    ax.set_ylabel('')
    ax.set_xlabel('')

    ax = axs[1]
    ax.set_title('ATAC assignment')
    atac_counts = demuxlet.ATAC.value_counts().rename('ATAC').rename_axis(index='index').reset_index()
    atac_counts['percentage'] = (100*atac_counts.ATAC/atac_counts.ATAC.sum()).map(lambda x: round(x, 2))
    atac_counts['label'] = atac_counts['index'] + ' - ' + atac_counts.percentage.astype(str) + '%'
    sns.barplot(x='ATAC', y='label', ax=ax, data=atac_counts, palette=cmap, hue='index', dodge=False)
    ax.legend().remove()
    ax.set_ylabel('')
    ax.set_xlabel('')

    ax = axs[2]
    ax.set_title('ATAC - RNA concordance')
    counts = (demuxlet.ATAC == demuxlet.RNA).map({True: 'ATAC == RNA', False: 'ATAC != RNA'}).value_counts().rename('counts').reset_index()
    sns.barplot(x='counts', y='index', ax=ax, data=counts)
    ax.set_ylabel('')
    ax.set_xlabel('')

    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-bar.png', dpi=300)
    fig.clf()

    fig, ax = plt.subplots(figsize=(10, 8))
    df = demuxlet.copy()
    df.label = append_n(df.label)
    sns.scatterplot(x='RNA_SNPs', y='ATAC_SNPs', hue='label', ax=ax, data=df, alpha=0.1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(read_count_formatter)
    ax.yaxis.set_major_formatter(read_count_formatter)
    ax.set_ylabel('ATAC SNPs checked')
    ax.set_xlabel('RNA SNPs checked')

    ax.legend(bbox_to_anchor=(1.05, 1.05)).set_title('')
    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-scatter.png', dpi=300)
    fig.clf()

elif args.atac_demuxlet is not None and args.rna_demuxlet is None:
    assert(args.strategy == 'ATAC' or args.strategy is None)

    atac_demuxlet = load_demuxlet(args.atac_demuxlet).assign(modality='ATAC')
    atac_demuxlet['assignment'] = atac_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)

    demuxlet = atac_demuxlet[['BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'ATAC', 'NUM.SNPS': 'ATAC_SNPs'})
    demuxlet['final_assignment'] = demuxlet['ATAC'].str.split(' ').map(lambda x: x[0])
    demuxlet[['BARCODE', 'final_assignment', 'ATAC']].to_csv(f'{args.prefix}assignments.txt', sep='\t', index=False, header=False)

    cmap = make_colormap_dict(list(sorted(set(demuxlet.ATAC.to_list()))))

    fig, ax = plt.subplots(figsize=(5, 5))

    atac_counts = demuxlet.ATAC.value_counts().rename('ATAC').rename_axis(index='index').reset_index()
    atac_counts['percentage'] = (100*atac_counts.ATAC/atac_counts.ATAC.sum()).map(lambda x: round(x, 2))
    atac_counts['label'] = atac_counts['index'] + ' - ' + atac_counts.percentage.astype(str) + '%'
    sns.barplot(x='ATAC', y='label', ax=ax, data=atac_counts, palette=cmap, hue='index', dodge=False)
    ax.legend().remove()
    ax.set_ylabel('')
    ax.set_xlabel('')

    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-bar.png', dpi=300)
    fig.clf()

    fig, ax = plt.subplots(figsize=(10, 8))
    df = demuxlet.copy()
    df['label'] = append_n(df.ATAC)
    sns.stripplot(y='label', x='ATAC_SNPs', ax=ax, data=df, alpha=0.1)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(read_count_formatter)
    ax.set_xlabel('ATAC SNPs checked')

    ax.legend(bbox_to_anchor=(1.05, 1.05)).set_title('')
    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-scatter.png', dpi=300)
    fig.clf()
elif args.atac_demuxlet is None and args.rna_demuxlet is not None:
    assert(args.strategy == 'RNA' or args.strategy is None)

    rna_demuxlet = load_demuxlet(args.rna_demuxlet).assign(modality='RNA')
    rna_demuxlet['assignment'] = rna_demuxlet[['DROPLET.TYPE', 'SNG.BEST.GUESS']].apply(lambda x: '{} ({})'.format(x[0], x[1]) if x[0] == 'SNG' else x[0], axis=1)
    rna_demuxlet['RNA_BARCODE'] = rna_demuxlet.BARCODE

    demuxlet = rna_demuxlet[['RNA_BARCODE', 'assignment', 'NUM.SNPS']].rename(columns={'assignment': 'RNA', 'NUM.SNPS': 'RNA_SNPs'})
    demuxlet['final_assignment'] = demuxlet['RNA'].str.split(' ').map(lambda x: x[0])
    demuxlet[['RNA_BARCODE', 'final_assignment', 'RNA']].to_csv(f'{args.prefix}assignments.txt', sep='\t', index=False, header=False)

    cmap = make_colormap_dict(list(sorted(set(demuxlet.RNA.to_list()))))

    fig, ax = plt.subplots(figsize=(5, 5))

    rna_counts = demuxlet.RNA.value_counts().rename('RNA').rename_axis(index='index').reset_index()
    rna_counts['percentage'] = (100*rna_counts.RNA/rna_counts.RNA.sum()).map(lambda x: round(x, 2))
    rna_counts['label'] = rna_counts['index'] + ' - ' + rna_counts.percentage.astype(str) + '%'
    sns.barplot(x='RNA', y='label', ax=ax, data=rna_counts, palette=cmap, hue='index', dodge=False)
    ax.legend().remove()
    ax.set_ylabel('')
    ax.set_xlabel('')

    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-bar.png', dpi=300)
    fig.clf()

    fig, ax = plt.subplots(figsize=(10, 8))
    df = demuxlet.copy()
    df['label'] = append_n(df.RNA)
    sns.stripplot(y='label', x='RNA_SNPs', ax=ax, data=df, alpha=0.1)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(read_count_formatter)
    ax.set_xlabel('RNA SNPs checked')

    ax.legend(bbox_to_anchor=(1.05, 1.05)).set_title('')
    fig.tight_layout()
    fig.savefig(f'{args.prefix}demuxlet-scatter.png', dpi=300)
    fig.clf()