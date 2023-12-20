#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd

AUTOSOMES = {
    'hg19': [f'chr{i}' for i in range(1, 23)],
    'hg38': [f'chr{i}' for i in range(1, 23)],
    'mm9': [f'chr{i}' for i in range(1, 20)],
    'mm10': [f'chr{i}' for i in range(1, 20)],
    'rn4': [f'chr{i}' for i in range(1, 21)],
    'rn5': [f'chr{i}' for i in range(1, 21)],
    'rn6': [f'chr{i}' for i in range(1, 21)],
    'rn7': [f'chr{i}' for i in range(1, 21)],
    'rn6_eGFP_ACADSB': [f'chr{i}' for i in range(1, 21)],
    'rn6_eGFP': [f'chr{i}' for i in range(1, 21)]
}

PASS_QC_ATAC_BARCODES, ALL_ATAC_BARCODES, GENOME = sys.argv[1:]

if not GENOME in AUTOSOMES:
    raise ValueError('Unsupported genome. Add genome autosome names to AUTOSOMES dictionary and re-run.')

pass_qc_barcodes = set(pd.read_csv(PASS_QC_ATAC_BARCODES, header=None).loc[:,0].to_list())
atac_barcode_whitelist = set(pd.read_csv(ALL_ATAC_BARCODES, header=None)[0].to_list())

assert(len(pass_qc_barcodes.intersection(atac_barcode_whitelist)) == len(pass_qc_barcodes))

with open('singlecell.csv', 'w') as f:
    for barcode in atac_barcode_whitelist:
        is_cell = int(barcode in pass_qc_barcodes)
        cell_name = barcode
        f.write(f'{barcode},{cell_name},{is_cell}\n')

# print the autosome list
with open('autosomes.txt', 'w') as f:
    for chrom in AUTOSOMES[GENOME]:
        f.write(f'{chrom}\n')