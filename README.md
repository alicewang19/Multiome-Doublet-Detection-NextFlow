# 10X snRNA-seq / multiome doublet detection pipeline

This pipeline runs demultiplexing with [demuxlet](https://github.com/statgen/popscle) (if VCF file provided), runs ATAC-based doublet detection using [AMULET](https://github.com/UcarLab/AMULET) (in the case of multiome data), and runs RNA-based doublet detection using [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder). It also performs an initial (unoptimized) clustering of libraries, individually and jointly. The joint clustering may be useful for identifying doublet clusters containing nuclei not otherwise identified as doublets.

It currently does not process non-multiome ATAC data; it handles sc/snRNA-seq, or multiome (snRNA+snATAC) data.

## Dependencies
[Singularity (v. 3)](https://docs.sylabs.io/guides/3.0/user-guide/) and [NextFlow](https://www.nextflow.io/) (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library) or Docker hub (https://hub.docker.com/).

## Running

You'll need to provide a TSV file like this, noting file paths and some parameters for each library (one library per line; header must be included):

```
library	atac_bam	rna_bam	atac_pass_qc_barcodes	rna_pass_qc_barcodes	rna_cellbender	vcf	doubletfinder_pcs	doubletfinder_resolution	doubletfinder_sctransform
9266-VD-1	/path/to/9266-VD-1.atac.bam	/path/to/9266-VD-1.rna.bam	/path/to/9266-VD-1.pass-qc-barcodes.atac.txt	/path/to/9266-VD-1.pass-qc-barcodes.rna.txt	/path/to/9266-VD-1-hg38.cellbender_FPR_0.05.h5	/path/to/genotypes.subset.vcf.gz	25	0.2	false
9266-VD-2	/path/to/9266-VD-2.atac.bam	/path/to/9266-VD-2.rna.bam	/path/to/9266-VD-2.pass-qc-barcodes.atac.txt	/path/to/9266-VD-2.pass-qc-barcodes.rna.txt	/path/to/9266-VD-2-hg38.cellbender_FPR_0.05.h5	/path/to/genotypes.subset.vcf.gz	25	0.2	false
```

Columns can be in any order (column names should not be changed), and represent:
1. Library ID
2. ATAC bam file (already filtered as desired for doublet detection; omit if you're only analyzing snRNA data)
3. RNA bam file (already filtered as desired for doublet detection)
4. List of pass QC barcodes for ATAC (omit if you're only analyzing snRNA data)
5. List of pass QC barcodes for RNA
6. Path to cellbender HDF5 output file
7. Path to VCF for demuxlet (omit if there is none, in which case demuxlet will not be run)
8. Number of PCs to use when running doubletfinder (set to a number that would be reasonable for Seurat clustering)
9. Resolution to use when running doubletfinder (set to a number that would be reasonable for Seurat clustering)
10. true / false, indicating whether to use the Seurat scTransform or standard workflow when running doubletfinder

You'll also need to update the reference file paths in the `nextflow.config` file in this directory.

Then run the pipeline. For multiome data:

```bin
nextflow run -resume --library_info library_info.tsv --atac_barcodes /path/to/snATACseq-NextFlow/737K-arc-v1.txt --rna_barcodes /path/to/snRNAseq-NextFlow/737K-arc-v1.txt --results /path/to/results -entry multiome /path/to/Multiome-Doublet-Detection-NextFlow/main.nf
```

Where `--atac_barcodes` and `--rna_barcodes` are the barcode whitelists for the two modalities, and `library_info.tsv` is the file noting file paths and parameters for each library. For sc/snRNA data:

```bin
nextflow run -resume --library_info library_info.tsv --rna_barcodes /path/to/snRNAseq-NextFlow/737K-arc-v1.txt --results /path/to/results -entry rna /path/to/Multiome-Doublet-Detection-NextFlow/main.nf
```


## Output
# TODO: update
* `atac-doublet-detection`: Doublet detection results based on AMULET, including lists of singlet RNA/ATAC barcodes (to be used for downstream analysis)
* `demuxlet/out`: Raw demuxlet output files for each library and modality
* `demuxlet/processed`: List of demuxlet --> doublet / singlet (+ individual) assignments, and plots of demuxlet assignments and ATAC / RNA concordance
* `doubletfinder`: Doublet detection results based on DoubletFinder
* `cluster`: Seurat clustering results, with doublets included
* `merged-counts`: Input RNA count matrix for joint clustering
* `pass-qc-nuclei-counts-with-doublets`: Input RNA count matrix for clustering of individual libraries
* `doublets-in-clustering`: Plots of doublet vs singlet nuclei in joint clustering
