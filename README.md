# 10X Multiome doublet detection pipeline

This pipeline runs demultiplexing with [demuxlet](https://github.com/statgen/popscle) (if VCF file provided), runs ATAC-based doublet detection using [AMULET](https://github.com/UcarLab/AMULET), and runs RNA-based doublet detection using [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder). It also performs an initial (unoptimized) clustering of libraries, individually and jointly. The joint clustering may be useful for identifying doublet clusters containing nuclei not previously identified as doublets.

## Dependencies
[Singularity (v. 3)](https://docs.sylabs.io/guides/3.0/user-guide/) and [NextFlow](https://www.nextflow.io/) (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library) or Docker hub (https://hub.docker.com/).

## Running

You'll need to provide a config.json file like this:

```python
{
    "libraries": {
        "9240-VD-1-hg38": {
            "atac_bam": "/path/to/multiome-atacseq/results/prune/9240-VD-1-hg38.pruned.bam",
            "atac_pass_qc_barcodes": "/path/to/list-of-pass-qc-atac-barcodes.txt",
            "genome": "hg38",
            "rna_bam": "/path/to/multiome-rnaseq/results/prune/9240-VD-1-hg38.before-dedup.bam",
            "rna_cellbender": "/path/to/multiome-rnaseq/results/cellbender/9240-VD-1-hg38.cellbender_FPR_0.05.h5",
            "rna_pass_qc_barcodes": "/path/to/list-of-pass-qc-rna-barcodes.txt",
            "vcf": "/path/to/donors.vcf.gz"
        }
    }
}
```

You'll also need to update the file paths in the `nextflow.config` file in this directory.

Then run the pipeline:

```bin
nextflow run -resume -params-file config.json --results /path/to/results /path/to/Multiome-Doublet-Detection-NextFlow/main.nf
```

You may wish to run under nohup so that the pipeline continues to run in the background and does not terminate upon logging out of the server (`nohup nextflow run ... &`)

## Output

* `atac-doublet-detection`: Doublet detection results based on AMULET, including lists of singlet RNA/ATAC barcodes (to be used for downstream analysis)
* `demuxlet/out`: Raw demuxlet output files for each library and modality
* `demuxlet/processed`: List of demuxlet --> doublet / singlet (+ individual) assignments, and plots of demuxlet assignments and ATAC / RNA concordance
* `doubletfinder`: Doublet detection results based on DoubletFinder
* `cluster`: Seurat clustering results, with doublets included
* `merged-counts`: Input RNA count matrix for joint clustering
* `pass-qc-nuclei-counts-with-doublets`: Input RNA count matrix for clustering of individual libraries
* `doublets-in-clustering`: Plots of doublet vs singlet nuclei in joint clustering
