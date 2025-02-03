#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def parseTsvFileToMap(filePath) {
    def result = [:]
    def lines = new File(filePath).readLines()
    def header = lines[0].split("\t")

    lines.drop(1).eachWithIndex { line, index ->
        def values = line.split("\t")
        def lineMap = [:]
        header.eachWithIndex { colName, colIndex ->
            lineMap[colName] = values[colIndex]
        }
        if (lineMap.containsKey('doubletfinder_sctransform')) {
            lineMap['doubletfinder_sctransform'] = lineMap['doubletfinder_sctransform'].toBoolean()
        }
        result[lineMap['library']] = lineMap
    }

    return result
}


def get_blacklists (genome) {
    if (params.blacklist == null) {
        return []
    } else if (!params.blacklist.containsKey(genome)) {
        return []
    } else if (params.blacklist[genome] instanceof String) {
        return [params.blacklist[genome]]
    } else {
        return params.blacklist[genome]
    }
}


def has_vcf(library, library_info) {
    if (library_info[library].containsKey('vcf')) {
        if (library_info[library].vcf == 'None' || library_info[library].vcf == '' || library_info[library].vcf == null) {
            return false
        } else {
            return true
        }
    } else {
        return false
    }
}

markers = params.markers
joint_clustering_pcs = params.joint_clustering_pcs
joint_clustering_resolution = params.joint_clustering_resolution
sctransform = params.sctransform
library_info = parseTsvFileToMap(params.library_info)


process make_cellbender_mm {

    memory '20 GB'
    tag "${library}"
    container 'docker://porchard/general:20220406125608'

    input:
    tuple val(library), path(h5)

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    cellbender-h5-to-mm.py $h5 ${library}.
    """

}


process subset_nuclei {

    publishDir "${params.results}/pass-qc-nuclei-counts-with-doublets"
    memory '30 GB'
    tag "${library}"
    container 'docker://porchard/mm:20230104'

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv'), path('keep-barcodes.txt')

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    cat keep-barcodes.txt barcodes.tsv | sort | uniq -d > keep-in-cellbender.txt
    mm subset --matrix matrix.mtx --features features.tsv --barcodes barcodes.tsv --keep-barcodes keep-in-cellbender.txt --prefix ${library}.
    """

}


process index_bam_atac {

    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path(bam), path(bam_index)

    script:
    bam_index = bam.getName() + ".bai"

    """
    samtools index $bam
    """

}


process index_bam_rna {

    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path(bam), path(bam_index)

    script:
    bam_index = bam.getName() + ".bai"

    """
    samtools index $bam
    """

}


// only if doing demuxlet
process clip_bam {

    memory '22 GB'
    cache 'lenient'
    tag "${library}"
    container 'library://porchard/default/general:20220107'
    time '24h'

    input:
    tuple val(library), path(bam), path(index)

    output:
    tuple val(library), path("${library}.clipped.bam"), path("${library}.clipped.bam.bai")

    when:
    has_vcf(library, library_info)

    """
    /sw/bamUtil/bin/bam clipOverlap --poolSize 9000000 --in $bam --out ${library}.clipped.bam
    samtools index ${library}.clipped.bam
    """

}


// for RNA: do demuxlet in chunks of e.g. 500 barcodes
// ATAC is fast even w/ e.g. 10k barcodes, so don't bother for ATAC
process chunk_demuxlet_barcodes {

    cache 'lenient'
    tag "${library}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), path('pass-qc-barcodes.txt')

    output:
    tuple val(library), path("${library}.barcode-batch.*.txt")

    """
    split --additional-suffix=.txt --lines 500 --suffix-length 10 pass-qc-barcodes.txt ${library}.barcode-batch.
    """

}


process demuxlet {

    container 'library://porchard/default/demuxlet:20220204'
    cache 'lenient'
    tag "${library} ${modality}"
    time '48h'
    memory { 20.GB * task.attempt }
    maxRetries 3

    input:
    tuple val(library), path(bam), path(bam_index), val(modality), path(barcodes), path(vcf)

    output:
    tuple val(library), val(modality), path("${library}-${modality}.best"), emit: best

    """
    popscle demuxlet --sam $bam --vcf $vcf --alpha 0 --alpha 0.5 --group-list $barcodes --field GT --out ${library}-${modality}
    """

}


process concat_demuxlet {

    publishDir "${params.results}/demuxlet/out"

    input:
    tuple val(library), val(modality), path("demuxlet.*.best")

    output:
    tuple val(library), val(modality), path("${library}-${modality}.best.txt")

    """
    cat demuxlet.*.best | awk 'NR==1' > header.txt
    cat header.txt > ${library}-${modality}.best.txt
    cat demuxlet.*.best | grep -v -f header.txt >> ${library}-${modality}.best.txt
    """

}


process process_demuxlet {

    publishDir "${params.results}/demuxlet/processed"
    container "library://porchard/default/general:20220107"

    input:
    tuple val(library), path(demuxlet_best_atac), path(demuxlet_best_rna), path(rna_barcodes), path(atac_barcodes), val(has_atac)

    output:
    tuple val(library), path("*assignments.txt"), emit: assignments
    path("*.png")

    script:
    atac_barcodes_flag = has_atac ? "--atac-barcodes ${atac_barcodes}" : ""
    atac_demuxlet_flag = has_atac ? "--atac-demuxlet ${demuxlet_best_atac}" : ""
    strategy = has_atac ? "joint" : "RNA"

    """
    process-demuxlet-out.py ${atac_demuxlet_flag} --rna-demuxlet $demuxlet_best_rna --strategy ${strategy} --prefix ${library}. ${atac_barcodes_flag} --rna-barcodes $rna_barcodes
    """

}


process prep_doublet_detection {

    tag "${library}"
    container "library://porchard/default/general:20220107"

    input:
    tuple val(library), path(pass_qc_barcodes), path(atac_barcode_list)

    output:
    tuple val(library), path('singlecell.csv'), path('autosomes.txt')

    """
    make-doubletdetector-files.py $pass_qc_barcodes $atac_barcode_list ${params.genome}
    """

}


process run_atac_doublet_detection {

    publishDir "${params.results}/amulet"
    tag "${library}"
    container 'library://porchard/default/amulet:1.1'
    memory '10 GB'

    input:
    tuple val(library), path(bam), path(bam_index), path(single_cell), path(autosomes), path(blacklists), val(use_blacklist)

    output:
    tuple val(library), path("${library}.doublet_probabilities.txt"), emit: probs
    tuple val(library), path("${library}")

    script:
    blacklist_cmd = use_blacklist ? "zcat ${blacklists.join(' ')} | sort -k1,1 -k2n,2 > blacklist.bed" : "touch blacklist.bed"

    """
    mkdir -p output
    ${blacklist_cmd}
    /opt/AMULET/AMULET.sh --bcidx 0 --cellidx 1 --iscellidx 2 $bam $single_cell $autosomes blacklist.bed output/ /opt/AMULET/
    cp -r output ${library}
    cp output/MultipletProbabilities.txt ${library}.doublet_probabilities.txt
    """

}


// doubletfinder
process cluster_per_library {

    publishDir "${params.results}/cluster/with-doublets/per-library", overwrite: true
    memory '15 GB'
    tag "${library}"
    container 'library://porchard/default/r-general:20220112'

    input:
    tuple val(library), path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    tuple val(library), path("${library}.clusters.txt"), emit: clusters

    script:
    pcs = library_info[library].doubletfinder_pcs
    resolution = library_info[library].doubletfinder_resolution
    sctransform_flag = library_info[library].doubletfinder_sctransform ? '--sctransform' : ''

    """
    seurat.R --markers ${markers.join(',')} --prefix ${library}. ${sctransform_flag} --pcs ${pcs} --resolution ${resolution} --matrix $matrix --features $features --barcodes $barcodes
    """

}


process doubletfinder {

    publishDir "${params.results}/doubletfinder", overwrite: true
    memory '20 GB'
    tag "${library}"
    container 'docker://porchard/doubletfinder:20230426105316'

    input:
    tuple val(library), path(matrix), path(features), path(barcodes)

    output:
    tuple val(library), path("${library}.doubletfinder_assignments.txt")

    script:
    pcs = library_info[library].doubletfinder_pcs
    resolution = library_info[library].doubletfinder_resolution
    sctransform_flag = library_info[library].doubletfinder_sctransform ? '--sctransform' : ''

    """
    doubletfinder.R --markers ${markers.join(',')} --prefix ${library}.doubletfinder_ ${sctransform_flag} --pcs ${pcs} --resolution ${resolution} --matrix $matrix --features $features --barcodes $barcodes
    """

}


process add_prefix_to_barcodes_before_postdecontamination_merge {

    memory '1 GB'
    tag "${library}"

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv')

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    ln -s matrix.mtx ${library}.matrix.mtx
    ln -s features.tsv ${library}.features.tsv
    cat barcodes.tsv | perl -pe 's/^/${library}-/' > ${library}.barcodes.tsv
    """

}


process merge_postdecontamination_matrices {

    publishDir "${params.results}/merged-counts/with-doublets", overwrite: true
    label 'largemem'
    memory '150 GB'
    container 'docker://porchard/mm:20230104'
    time '2h'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(x), path(matrices), path(features), path(barcodes)

    output:
    tuple path("merged.matrix.mtx"), path("merged.features.tsv"), path("merged.barcodes.tsv")

    """
    mm merge --matrices ${matrices.join(' ')} --features ${features.join(' ')} --barcodes ${barcodes.join(' ')} --prefix merged.
    """

}


process cluster_joint_with_doublets {

    publishDir "${params.results}/cluster/with-doublets/joint", overwrite: true
    label 'largemem'
    memory '100 GB'
    time '5h'
    container 'library://porchard/default/r-general:20220112'
    errorStrategy 'ignore'

    input:
    tuple path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    path("*.rds")
    path("merged.umap.txt"), emit: umap
    path("merged.clusters.txt"), emit: clusters

    script:
    sctransform_flag = sctransform ? '--sctransform' : ''

    """
    seurat.R ${sctransform_flag} --markers ${markers.join(',')} --prefix merged. --pcs $joint_clustering_pcs --resolution $joint_clustering_resolution --matrix $matrix --features $features --barcodes $barcodes
    """

}


process plot_doublets_in_clustering {

    publishDir "${params.results}/doublets-in-clustering"
    memory '25 GB'
    time '2h'
    container 'docker://porchard/general:20220406125608'

    input:
    path(umap)
    path(clusters)
    path(rna_barcodes)
    path(atac_barcodes)
    path(doubletfinder_out)
    path(amulet_out)
    path(demuxlet_out)
    val(use_amulet)
    val(use_demuxlet)

    output:
    path("*.png")

    script:
    demuxlet_flag = !use_demuxlet ? '' : "--demuxlet ${demuxlet_out.join(' ')}"
    amulet_flag = !use_amulet ? '' : "--amulet ${amulet_out.join(' ')}"
    has_atac = use_amulet
    atac_barcodes_flag = !has_atac ? '' : "--atac-barcodes $atac_barcodes"

    """
    plot-doublets-in-clustering.py --prefix plot-doublets. --umap $umap --clusters $clusters --doubletfinder ${doubletfinder_out.join(' ')} ${amulet_flag} ${demuxlet_flag} ${atac_barcodes_flag} --rna-barcodes $rna_barcodes
    """

}


workflow multiome {

    libraries = library_info.keySet()

    atac_bam = Channel.from(libraries.collect({it -> [it, file(library_info[it].atac_bam)]})) | index_bam_atac // library, bam, index
    rna_bam = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_bam)]})).filter({it -> has_vcf(it[0], library_info)}) | index_bam_rna // library, bam, index
    rna_cellbender = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_cellbender)]})) // library, RNA cellbender output
    rna_pass_qc_barcodes = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_pass_qc_barcodes)]})) // library, rna pass QC barcodes
    atac_pass_qc_barcodes = Channel.from(libraries.collect({it -> [it, file(library_info[it].atac_pass_qc_barcodes)]})) // library, atac pass QC barcodes

    nuclei = make_cellbender_mm(rna_cellbender).combine(rna_pass_qc_barcodes, by: 0) | subset_nuclei

    // demuxlet
    demuxlet_in_atac = clip_bam(atac_bam).map({it -> it + ['ATAC']}).combine(atac_pass_qc_barcodes, by: 0)
    demuxlet_rna_barcode_chunks = chunk_demuxlet_barcodes(rna_pass_qc_barcodes).transpose()
    demuxlet_in_rna = rna_bam.map({it -> it + ['RNA']}).combine(demuxlet_rna_barcode_chunks, by: 0)
    demuxlet_out = (demuxlet_in_atac.mix(demuxlet_in_rna).map({it -> it + [file(library_info[it[0]].vcf)]}) | demuxlet).best.groupTuple(by: [0, 1]) | concat_demuxlet // library, modality, best
    demuxlet_out_atac = demuxlet_out.filter({it -> it[1] == 'ATAC'}).map({it -> [it[0], it[2]]}) // library, best
    demuxlet_out_rna = demuxlet_out.filter({it -> it[1] == 'RNA'}).map({it -> [it[0], it[2]]}) // library, best

    demuxlet_assignments = demuxlet_out_atac.combine(demuxlet_out_rna, by: 0).combine(Channel.fromPath(params.rna_barcodes)).combine(Channel.fromPath(params.atac_barcodes)).combine(Channel.of(true)) | process_demuxlet

    // amulet
    amulet = atac_bam.combine(prep_doublet_detection(atac_pass_qc_barcodes.combine(Channel.fromPath(params.atac_barcodes))), by: 0).map({it -> it + [get_blacklists(params.genome).collect({x -> file(x)}), get_blacklists(params.genome) != []]}) | run_atac_doublet_detection

    // doubletfinder
    cluster_per_library(nuclei)
    df = doubletfinder(nuclei)

    // cluster, inc. doublets
    joint_clustering_in = add_prefix_to_barcodes_before_postdecontamination_merge(nuclei).map({it -> ['x'] + it[1..3]}).groupTuple() | merge_postdecontamination_matrices
    joint_clustering_round_1 = cluster_joint_with_doublets(joint_clustering_in)

    // check if VCF was provided for all libraries
    all_libraries_have_vcf = libraries.every({it -> has_vcf(it, library_info)})
    demuxlet_for_plotting = all_libraries_have_vcf ? demuxlet_assignments.assignments.map({it -> it[1]}).toSortedList() : Channel.fromList([file('NULL_DEMUXLET')]).toSortedList()
    plot_doublets_in_clustering(joint_clustering_round_1.umap, joint_clustering_round_1.clusters, Channel.fromPath(params.rna_barcodes), Channel.fromPath(params.atac_barcodes), df.map({it -> it[1]}).toSortedList(), amulet.probs.map({it -> it[1]}).toSortedList(), demuxlet_for_plotting, true, all_libraries_have_vcf)

}

workflow rna {

    libraries = library_info.keySet()

    rna_bam = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_bam)]})).filter({it -> has_vcf(it[0], library_info)}) | index_bam_rna // library, bam, index
    rna_cellbender = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_cellbender)]})) // library, RNA cellbender output
    rna_pass_qc_barcodes = Channel.from(libraries.collect({it -> [it, file(library_info[it].rna_pass_qc_barcodes)]})) // library, rna pass QC barcodes

    nuclei = make_cellbender_mm(rna_cellbender).combine(rna_pass_qc_barcodes, by: 0) | subset_nuclei

    // demuxlet
    demuxlet_rna_barcode_chunks = chunk_demuxlet_barcodes(rna_pass_qc_barcodes).transpose()
    demuxlet_in_rna = rna_bam.map({it -> it + ['RNA']}).combine(demuxlet_rna_barcode_chunks, by: 0)
    demuxlet_out = (demuxlet_in_rna.map({it -> it + [file(library_info[it[0]].vcf)]}) | demuxlet).best.groupTuple(by: [0, 1]) | concat_demuxlet // library, modality, best
    demuxlet_out_rna = demuxlet_out.filter({it -> it[1] == 'RNA'}).map({it -> [it[0], it[2]]}) // library, best
    demuxlet_out_atac = demuxlet_out_rna.map({it -> [it[0], file('NULL_ATAC_DEMUXLET')]})

    demuxlet_assignments = demuxlet_out_atac.combine(demuxlet_out_rna, by: 0).combine(Channel.fromPath(params.rna_barcodes)).combine(Channel.fromList([file('NULL_ATAC_BARCODES')])).combine(Channel.of(false)) | process_demuxlet

    // doubletfinder
    cluster_per_library(nuclei)
    df = doubletfinder(nuclei)

    // cluster, inc. doublets
    joint_clustering_in = add_prefix_to_barcodes_before_postdecontamination_merge(nuclei).map({it -> ['x'] + it[1..3]}).groupTuple() | merge_postdecontamination_matrices
    joint_clustering_round_1 = cluster_joint_with_doublets(joint_clustering_in)

    // check if VCF was provided for all libraries
    all_libraries_have_vcf = libraries.every({it -> has_vcf(it, library_info)})
    demuxlet_for_plotting = all_libraries_have_vcf ? demuxlet_assignments.assignments.map({it -> it[1]}).toSortedList() : Channel.fromList([file('NULL_DEMUXLET')]).toSortedList()
    amulet_for_plotting = Channel.fromList([file('NULL_AMULET')]).toSortedList()
    plot_doublets_in_clustering(joint_clustering_round_1.umap, joint_clustering_round_1.clusters, Channel.fromPath(params.rna_barcodes), Channel.fromList([file('NULL_ATAC_BARCODES')]), df.map({it -> it[1]}).toSortedList(), amulet_for_plotting, demuxlet_for_plotting, false, all_libraries_have_vcf)

}
