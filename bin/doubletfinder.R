#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(glue)
library(optparse)
library(DoubletFinder)

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
  make_option(c("--markers"), action = 'store', type = 'character', default = '', help = '[Optional] Comma-separated list of marker genes to visualize'),
  make_option(c("--prefix"), action = 'store', type = 'character', default = 'seurat.', help = '[Optional] Prefix of output files')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)


#RNA_MTX <- '/lab/work/porchard/2022-01-rat-multiome/work/clean/results/pass-qc-nuclei-counts/5203-NM-2-rn6_eGFP_ACADSB.matrix.mtx'
#RNA_FEATURES <- '/lab/work/porchard/2022-01-rat-multiome/work/clean/results/pass-qc-nuclei-counts/5203-NM-2-rn6_eGFP_ACADSB.features.tsv'
#RNA_BARCODES <- '/lab/work/porchard/2022-01-rat-multiome/work/clean/results/pass-qc-nuclei-counts/5203-NM-2-rn6_eGFP_ACADSB.barcodes.tsv'
#RESOLUTION <- 0.1
#PCS <- 25
#MARKERS <- opts$markers

RNA_MTX <- opts$matrix
RNA_FEATURES <- opts$features
RNA_BARCODES <- opts$barcodes
RESOLUTION <- opts$resolution
PCS <- opts$pcs
max_dims <- PCS
PREFIX <- opts$prefix
MARKERS <- strsplit(opts$markers, ',')[[1]]
SCTRANSFORM <- opts$sctransform

mm <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)

# restrict to protein-coding, autosomal genes?
#DROP <- grep('\\(LOC', rownames(mm), value=T)
#DROP <- c(DROP, grep('\\(AAB', rownames(mm), value=T))
#DROP <- c(DROP, grep('\\(RF', rownames(mm), value=T))
#DROP <- c(DROP, grep('\\(Mt', rownames(mm), value=T))
#mm <- mm[!rownames(mm) %in% DROP,]
#length(DROP)


#metadata <- data.frame(nucleus=colnames(mm), library=gsub('_.*', '', colnames(mm)))
#rownames(metadata) <- metadata$nucleus

rna <- CreateSeuratObject(counts = mm, min.cells=5, min.features=50)

#### either use scTransform...
if (SCTRANSFORM) {
    # regressing out % mito
    #rna <- PercentageFeatureSet(rna, pattern = '\\(MT-', col.name = "percent.mt")
    #VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    #rna <- SCTransform(rna, vars.to.regress = "percent.mt", verbose = FALSE)

    # or not
    rna <- SCTransform(rna, verbose = FALSE)
} else {
    rna <- NormalizeData(rna, verbose=F)
    rna <- FindVariableFeatures(rna, selection.method='vst', nfeatures=2000, verbose=F)
    rna <- ScaleData(rna, verbose=F)
}



rna <- RunPCA(rna, npcs=50, verbose=F)

png(glue('{PREFIX}elbow-plot.png'), height=5, width=5, units='in', res=300)
ElbowPlot(rna, ndims=50)
dev.off()

for(i in seq(1, PCS)){
    png(glue('{PREFIX}pc-loadings.{i}.png'), height=10, width=10, units='in', res=300)
    print(VizDimLoadings(rna, dims = i, reduction = "pca"))
    dev.off()
}


rna <- RunUMAP(rna, dims=1:PCS)
rna <- FindNeighbors(rna, dims = 1:PCS, k.param = 20)
rna <- FindClusters(rna, resolution = RESOLUTION, n.start = 100)

#MARKERS <- c('Myh7', 'Myh1', 'Myh2', 'Vwf', 'Pax7', 'Ptprc', 'Myh11', 'Myh3', 'Myh4', 'Fbn1', 'Pdgfra', 'Acta2')
if (length(MARKERS) > 0) {
    PLOT_FEATURES <- unlist(lapply(glue('\\({MARKERS}\\)'), function(x){grep(x, rownames(mm), value=T, ignore.case=T)}))
    #png(glue('{PREFIX}marker-genes.png'), height=4*sqrt(length(PLOT_FEATURES)), width=1+4*sqrt(length(PLOT_FEATURES)), units='in', res=300)
    #print(FeaturePlot(rna, PLOT_FEATURES, pt.size = 0.1, order=F))
    #dev.off()
    for(i in PLOT_FEATURES) {
        tryCatch(
            {png(glue('{PREFIX}markers.{i}.png'), height=6, width=7, units='in', res=300)
            print(FeaturePlot(rna, i, pt.size = 0.1, order=F))
            dev.off()},
            warning=function(x){print('failed')},
            error=function(x){print('failed')}
        )
        tryCatch(
            {png(glue('{PREFIX}markers.{i}.violin.png'), height=6, width=7, units='in', res=300)
            print(VlnPlot(rna, features=i))
            dev.off()},
            warning=function(x){print('failed')},
            error=function(x){print('failed')}
        )
    }
}






# Run DoubletFinder
## Estimate DoubletFinder parameters from data
sweep.res.list <- paramSweep_v3(rna, PCs = 1:max_dims, sct = SCTRANSFORM, num.cores = 1)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
# TODO: make plot of pK??
## Make a list with all the parameters
DF_pars <- list()
DF_pars$pN = 0.25
DF_pars$pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  mutate(pK = as.numeric(as.character(pK))) %>%
  pull(pK)
DF_pars$nExp_poi <- round(0.075 * length(rna@meta.data$orig.ident))
homotypic.prop <- modelHomotypic(rna@meta.data$seurat_clusters)
DF_pars$nExp_poi.adj <- round(DF_pars$nExp_poi * (1 - homotypic.prop))
## DoubletFinder appends columns to the metadata named after the parameters
pANN_lab = paste("pANN", DF_pars$pN, DF_pars$pK, DF_pars$nExp_poi, sep = "_")
DF_col <- gsub("pANN", "DF.classifications", pANN_lab)
DF_col2 <- paste("DF.classifications", DF_pars$pN, DF_pars$pK,
                 DF_pars$nExp_poi.adj, sep = "_")
## Flag doublets with two levels of stringency (high and lower, resp.)
rna <- doubletFinder_v3(rna, PCs = 1:max_dims, pN = DF_pars$pN,
                         pK = DF_pars$pK, nExp = DF_pars$nExp_poi,
                         reuse.pANN = FALSE, sct = SCTRANSFORM)
rna <- doubletFinder_v3(rna, PCs = 1:max_dims, pN = DF_pars$pN,
                         pK = DF_pars$pK, nExp = DF_pars$nExp_poi.adj,
                         reuse.pANN = pANN_lab, sct = SCTRANSFORM)


# Format data, collect results, and make plots
rna@meta.data[,"df"] <- rna@meta.data[,DF_col]
rna@meta.data[,"df2"] <- rna@meta.data[,DF_col2]

md <- rna@meta.data
singlets <- rownames(md)[md$df == "Singlet"]
assignments <- data.frame(nucleus=rownames(md), assignment=md$df)
assignments$assignment <- as.character(assignments$assignment)
assignments$assignment[assignments$assignment=='Singlet'] <- 'singlet'
assignments$assignment[assignments$assignment=='Doublet'] <- 'doublet'
write.table(assignments[,c('nucleus', 'assignment')], file=glue('{PREFIX}assignments.txt'), quote = F, append = F, sep = '\t', row.names = F, col.names = F)

## Count flagged doublets
# doublets_summary <- data.frame(
#   sample = sample,
#   homotypic.prop,
#   pK = DF_pars$pK,
#   high_stringency = sum(rna@meta.data$df == "Doublet"),
#   low_stringency = sum(rna@meta.data$df2 == "Doublet")
# )
# outfile <- glue("{PREFIX}doublet_summary.txt")
# write.table(doublets_summary, outfile, sep = "\t", quote = F, row.names = F)