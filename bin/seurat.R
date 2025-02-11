#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(glue)
library(optparse)

load_mm <- function(matrix_file, features_file, barcodes_file) {
    tmp <- as(Matrix::readMM(matrix_file), 'dgCMatrix')
    features <- read.table(features_file, as.is=T, sep='\t', head=F)
    features <- paste0(features$V1, ' (', features$V2, ')')
    barcodes <- read.table(barcodes_file, as.is=T, head=F)[,1]
    dimnames(tmp) <- list(features, barcodes)
    return(tmp)
}

option_list <- list(
  make_option(c("--matrix"), action = 'store', type = 'character', help = '[Required] Matrix'),
  make_option(c("--features"), action = 'store', type = 'character', help = '[Required] Features'),
  make_option(c("--barcodes"), action = 'store', type = 'character', help = '[Required] Barcodes'),
  make_option(c("--resolution"), action = 'store', type = 'numeric', default = 0.1, help = '[Optional] Resolution to use in clustering (default: 0.1)'),
  make_option(c("--pcs"), action = 'store', type = 'numeric', default = 20, help = '[Optional] Number of top PCs to use in clustering (default: 20)'),
  make_option(c("--sctransform"), action = 'store_true', type = 'logical', default = F, help = '[Optional] Apply scTransform'),
  make_option(c("--nomarkers"), action = 'store_true', type = 'logical', default = F, help = '[Optional] Skip finding marker genes'),
  make_option(c("--markers"), action = 'store', type = 'character', default = '', help = '[Optional] Comma-separated list of marker genes to visualize'),
  make_option(c("--prefix"), action = 'store', type = 'character', default = 'seurat.', help = '[Optional] Prefix of output files')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)


RNA_MTX <- opts$matrix
RNA_FEATURES <- opts$features
RNA_BARCODES <- opts$barcodes
RESOLUTION <- opts$resolution
PCS <- opts$pcs
PREFIX <- opts$prefix
MARKERS <- strsplit(opts$markers, ',')[[1]]
SCTRANSFORM <- opts$sctransform
GET_MARKERS <- !opts$nomarkers

mm <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)

rna <- CreateSeuratObject(counts = mm, min.cells=5, min.features=5, assay = "RNA", project='RNA')

if (SCTRANSFORM) {
    rna <- SCTransform(rna, verbose = FALSE)
} else {
    rna <- NormalizeData(rna, verbose=F)
    rna <- FindVariableFeatures(rna, selection.method='vst', nfeatures=2000, verbose=F)
    rna <- ScaleData(rna, verbose=F)
}


rna <- RunPCA(rna, npcs=200, verbose=F)

png(glue('{PREFIX}elbow-plot.png'), height=5, width=5, units='in', res=300)
ElbowPlot(rna, ndims=max(c(50, PCS*2)))
dev.off()

for(i in seq(1, PCS)){
    png(glue('{PREFIX}pc-loadings.{i}.png'), height=10, width=10, units='in', res=300)
    print(VizDimLoadings(rna, dims = i, reduction = "pca"))
    dev.off()
}


rna <- RunUMAP(rna, reduction='pca', dims=1:PCS)
rna <- FindNeighbors(rna, dims = 1:PCS, k.param = 20)
rna <- FindClusters(rna, resolution = RESOLUTION, n.start = 100)

if (length(MARKERS) > 0) {
    PLOT_FEATURES <- unlist(lapply(glue('\\({MARKERS}\\)'), function(x){grep(x, rownames(mm), value=T, ignore.case=T)}))
    for(i in PLOT_FEATURES) {
        tryCatch(
            {png(gsub(' ', '_', glue('{PREFIX}markers.{i}.scatter.png')), height=6, width=7, units='in', res=300)
            print(FeaturePlot(rna, i, pt.size = 0.1, order=F))
            dev.off()},
            warning=function(x){print(glue('failed to plot gene {i}'))},
            error=function(x){print(glue('failed to plot gene {i}'))}
        )
        tryCatch(
            {png(gsub(' ', '_', glue('{PREFIX}markers.{i}.violin.png')), height=6, width=7, units='in', res=300)
            print(VlnPlot(rna, features=i))
            dev.off()},
            warning=function(x){print(glue('failed to plot gene {i}'))},
            error=function(x){print(glue('failed to plot gene {i}'))}
        )
    }
}


png(glue('{PREFIX}clusters.png'), height=6, width=1+6, units='in', res=300)
DimPlot(rna, reduction = "umap")
dev.off()


clusters = as.data.frame(rna@active.ident)
colnames(clusters) <- c('cluster')
clusters$nucleus <- rownames(clusters)
clusters$barcode <- gsub('.*-(.*)', '\\1', clusters$nucleus)
write.table(clusters[,c('nucleus', 'cluster')], file = glue('{PREFIX}clusters.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = F)

umap <- as.data.frame(rna@reductions$umap@cell.embeddings)
colnames(umap) <- c('dim1', 'dim2')
umap$nucleus <- rownames(umap)
write.table(umap[,c('nucleus', 'dim1', 'dim2')], file=glue('{PREFIX}umap.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = F)

saveRDS(rna, glue("{PREFIX}rna.rds"))

if (GET_MARKERS) {
    # find cluster markers
    cluster_markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(cluster_markers, file=glue('{PREFIX}cluster-markers.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}