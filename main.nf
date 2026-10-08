#!/usr/bin/env nextflow
/*
 * SIP: a Single-cell Interchangeable Pipeline
 * Chen et al., bioRxiv 10.1101/456772 (2018). Section numbers below refer to the paper.
 */

def checkParams() {
    [readqc: ['fastp', 'fastqc', 'afterqc'], quant: ['salmon', 'star', 'hisat2'],
     count: ['featurecounts', 'htseq', 'stringtie'], strandness: ['auto', 'unstranded', 'forward', 'reverse']].each { k, v ->
        if (!(params[k] in v)) error "--${k} must be one of ${v}, got '${params[k]}'"
    }
    def units = params.quant == 'salmon' ? ['count', 'tpm'] : params.count == 'stringtie' ? ['tpm', 'fpkm'] : ['count']
    if (!(params.quant_unit in units)) error "--quant_unit for this quant/count combination must be one of ${units}"
    if (!params.cellhome || !params.gtf) error "--cellhome and --gtf are required"
    def ref = params.quant == 'salmon' ? 'transcript_fasta' : 'genome_fasta'
    if (!params["${params.quant}_index"] && !params[ref]) error "--quant ${params.quant} needs --${params.quant}_index or --${ref}"
}

// one-line run settings for the summary report, each item shell-quoted
def reportSettings(strand, runName) {
    def s = [run: runName, cellhome: params.cellhome, readqc: params.readqc, quant: params.quant]
    if (params.quant != 'salmon') s.count = params.count
    s.quant_unit = params.quant_unit
    s.strandness = params.strandness == 'auto' ? "${strand} (inferred with RSeQC)" : "${strand} (set by user)"
    if (params.using_seurat)
        s.seurat = "min features ${params.min_features}, max mito ${params.max_mt}%, regressed ${params.regress_vars}, " +
                   "t-SNE PCs ${params.tsne_dims}, perplexity ${params.tsne_perplexity}, clustering PCs ${params.cluster_dims}, resolution ${params.resolution}"
    s.monocle2 = params.using_monocle ? 'yes' : 'no'
    s.scater = params.using_scater ? 'yes' : 'no'
    s.collect { k, v -> "'${k}=${v}'" }.join(' ')
}

// the alphabetically first cell, so the strandness probe does not depend on task completion order
def firstCell(ch) {
    ch.toSortedList { a, b -> a[0].id <=> b[0].id }.map { l -> l[0] }
}

// RSeQC infer_experiment.py output -> forward / reverse / unstranded
def strandFrom(txt) {
    def t = txt.text
    def fw = fraction(t, /"(?:\+\+,--|1\+\+,1--,2\+-,2-\+)": ([\d.]+)/)
    def rv = fraction(t, /"(?:\+-,-\+|1\+-,1-\+,2\+\+,2--)": ([\d.]+)/)
    fw > 0.8 * (fw + rv) ? 'forward' : rv > 0.8 * (fw + rv) ? 'reverse' : 'unstranded'
}

def fraction(text, pattern) {
    def m = text =~ pattern
    m.find() ? m.group(1) as double : 0.0
}

// ---------- 2.1.1 read QC and preprocessing ----------
process FASTP {
    tag "$meta.id"
    label 'mid'
    publishDir "${params.outdir}/readqc", mode: 'copy', pattern: '*.{json,html}'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path('*.clean*.fastq.gz'), emit: reads
    path '*.fastp.{json,html}', emit: report
    script:
    def io = meta.single_end ? "-i $reads -o ${meta.id}.clean.fastq.gz"
                             : "-i ${reads[0]} -I ${reads[1]} -o ${meta.id}.clean_1.fastq.gz -O ${meta.id}.clean_2.fastq.gz"
    "fastp $io -w $task.cpus -j ${meta.id}.fastp.json -h ${meta.id}.fastp.html"
}

process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/readqc", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    path '*_fastqc.{zip,html}'
    script:
    "fastqc -t $task.cpus $reads"
}

process AFTERQC {
    tag "$meta.id"
    publishDir "${params.outdir}/readqc", mode: 'copy', pattern: 'QC/*'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path('good/*'), emit: reads
    path 'QC/*', emit: report
    script:
    def io = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    "after.py $io -g good -b bad -r QC"
}

// ---------- references ----------
process GTF_TABLES {
    input:
    path gtf
    output:
    path 'genes.gtf', emit: gtf
    path 'tx2gene.tsv', emit: tx2gene
    path 'genes.tsv', emit: genes
    path 'genes.bed12', emit: bed
    script:
    "gtf2tables.py $gtf && gzip -cdf $gtf > genes.gtf"
}

process SALMON_INDEX {
    label 'high'
    input:
    path fasta
    output:
    path 'salmon_index'
    script:
    "salmon index -t $fasta -i salmon_index -p $task.cpus"
}

process STAR_INDEX {
    label 'high'
    input:
    path fasta
    path gtf
    output:
    path 'star_index'
    script:
    """
    gzip -cdf $fasta > genome.fa
    n=\$(awk '!/^>/{s+=length(\$0)} END{v=int(log(s)/log(2)/2-1); print (v<14?v:14)}' genome.fa)
    mkdir star_index
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome.fa \\
         --sjdbGTFfile $gtf --genomeSAindexNbases \$n --runThreadN $task.cpus
    """
}

process HISAT2_INDEX {
    label 'high'
    input:
    path fasta
    output:
    path 'hisat2_index'
    script:
    "mkdir hisat2_index && gzip -cdf $fasta > genome.fa && hisat2-build -p $task.cpus genome.fa hisat2_index/genome"
}

// ---------- 2.1.2 strandness detection with RSeQC ----------
// Quasi-alignment mode: map the first reads of the first cell with salmon; the SAM header gives a
// transcript-coordinate BED (every transcript on '+') for infer_experiment.py.
process SALMON_PROBE {
    input:
    tuple val(meta), path(reads)
    path index
    output:
    tuple path('probe.sam'), path('probe.bed12')
    script:
    def n = params.strand_reads * 4
    def fq = meta.single_end ? [reads] : reads
    def sub = fq.withIndex().collect { f, i -> "gzip -cdf $f | head -n $n > sub_${i}.fq" }.join('\n    ')
    def io = meta.single_end ? '-r sub_0.fq' : '-1 sub_0.fq -2 sub_1.fq'
    """
    $sub
    salmon quant -i $index -l A $io -p $task.cpus -o probe --writeMappings=probe.sam
    awk -F'\\t' -v OFS='\\t' '\$1=="@SQ"{n=substr(\$2,4); l=substr(\$3,4); print n,0,l,n,0,"+",0,l,0,1,l",","0,"}' probe.sam > probe.bed12
    """
}

process INFER_STRAND {
    publishDir "${params.outdir}/strandness", mode: 'copy'
    input:
    tuple path(aln), path(bed)
    output:
    path 'infer_experiment.txt'
    script:
    "infer_experiment.py -i $aln -r $bed -q 0 > infer_experiment.txt"
}

// ---------- 2.1.2 (quasi-)alignment and transcript assembly ----------
process SALMON_QUANT {
    tag "$meta.id"
    label 'mid'
    publishDir "${params.outdir}/quant", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path index
    val strand
    output:
    path "${meta.id}"
    script:
    def lib = (meta.single_end ? '' : 'I') + [unstranded: 'U', forward: 'SF', reverse: 'SR'][strand]
    def io = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    "salmon quant -i $index -l $lib $io -p $task.cpus -o ${meta.id}"
}

process STAR_ALIGN {
    tag "$meta.id"
    label 'high'
    input:
    tuple val(meta), path(reads)
    path index
    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "${meta.id}.Log.final.out", emit: log
    script:
    def gz = reads.toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''
    """
    STAR --genomeDir $index --readFilesIn $reads $gz \\
         --runThreadN $task.cpus --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \\
         --outFileNamePrefix ${meta.id}.
    mv ${meta.id}.Aligned.sortedByCoord.out.bam ${meta.id}.bam
    """
}

process HISAT2_ALIGN {
    tag "$meta.id"
    label 'mid'
    input:
    tuple val(meta), path(reads)
    path index
    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "${meta.id}.hisat2.log", emit: log
    script:
    def io = meta.single_end ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    idx=\$(ls $index/*.1.ht2 | sed 's/\\.1\\.ht2\$//')
    hisat2 -x \$idx $io --dta -p $task.cpus --new-summary --summary-file ${meta.id}.hisat2.log | samtools sort -@ $task.cpus -o ${meta.id}.bam -
    """
}

process STRINGTIE {
    tag "$meta.id"
    label 'mid'
    input:
    tuple val(meta), path(bam)
    path gtf
    val strand
    output:
    path "${meta.id}.assembled.gtf"
    script:
    def s = [unstranded: '', forward: '--fr', reverse: '--rf'][strand]
    "stringtie $bam -G $gtf $s -p $task.cpus -o ${meta.id}.assembled.gtf"
}

process STRINGTIE_MERGE {
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
    path gtfs
    path gtf
    output:
    path 'merged.gtf'
    script:
    "stringtie --merge -G $gtf -o merged.gtf $gtfs"
}

// ---------- 2.1.3 gene-level abundance counting ----------
process STRINGTIE_ABUND {
    tag "$meta.id"
    label 'mid'
    input:
    tuple val(meta), path(bam)
    path merged
    val strand
    output:
    path "${meta.id}.gene_abund.tab"
    script:
    def s = [unstranded: '', forward: '--fr', reverse: '--rf'][strand]
    "stringtie $bam -e -G $merged $s -p $task.cpus -A ${meta.id}.gene_abund.tab -o ${meta.id}.estimated.gtf"
}

process FEATURECOUNTS {
    tag "$meta.id"
    input:
    tuple val(meta), path(bam)
    path gtf
    val strand
    output:
    path "${meta.id}.featureCounts.txt", emit: counts
    path "${meta.id}.featureCounts.txt.summary", emit: summary
    script:
    def s = [unstranded: 0, forward: 1, reverse: 2][strand]
    "featureCounts -a $gtf -o ${meta.id}.featureCounts.txt -s $s ${meta.single_end ? '' : '-p --countReadPairs'} -T $task.cpus $bam"
}

process HTSEQ {
    tag "$meta.id"
    input:
    tuple val(meta), path(bam)
    path gtf
    val strand
    output:
    path "${meta.id}.htseq.txt"
    script:
    def s = [unstranded: 'no', forward: 'yes', reverse: 'reverse'][strand]
    "htseq-count -f bam -r pos -s $s -t exon -i gene_id $bam $gtf > ${meta.id}.htseq.txt"
}

process TXIMPORT {
    publishDir "${params.outdir}/matrix", mode: 'copy'
    input:
    path dirs
    path tx2gene
    output:
    path "gene_${params.quant_unit}.tsv"
    script:
    "tximport.R $tx2gene ${params.quant_unit} $dirs"
}

process MERGE_TABLES {
    publishDir "${params.outdir}/matrix", mode: 'copy'
    input:
    path tables
    output:
    path "gene_${params.quant_unit}.tsv", emit: matrix
    path 'genes.tsv', optional: true, emit: genes
    script:
    "merge_tables.R ${params.count} ${params.quant_unit} $tables"
}

// ---------- reports ----------
process MULTIQC {
    publishDir params.outdir, mode: 'copy'
    input:
    path reports
    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_data', emit: data
    script:
    "multiqc -f ."
}

process REPORT {
    publishDir params.outdir, mode: 'copy'
    input:
    path results
    val settings
    output:
    path 'sip_report.pdf'
    script:
    "report.R $settings"
}

// ---------- 2.1.4 downstream analysis ----------
process SEURAT {
    publishDir "${params.outdir}/seurat", mode: 'copy'
    input:
    path matrix
    path genes
    output:
    path '*'
    script:
    """
    seurat.R matrix=$matrix genes=$genes min_cells=${params.min_cells} min_features=${params.min_features} \\
        max_mt=${params.max_mt} regress_vars=${params.regress_vars} n_pcs=${params.n_pcs} tsne_dims=${params.tsne_dims} tsne_perplexity=${params.tsne_perplexity} \\
        cluster_dims=${params.cluster_dims} resolution=${params.resolution} de_min_pct=${params.de_min_pct} heatmap_genes=${params.heatmap_genes}
    """
}

process MONOCLE2 {
    publishDir "${params.outdir}/monocle2", mode: 'copy'
    input:
    path matrix
    path genes
    output:
    path '*'
    script:
    "monocle2.R $matrix $genes ${params.quant_unit}"
}

process SCATER {
    publishDir "${params.outdir}/scater", mode: 'copy'
    input:
    path matrix
    path genes
    output:
    path '*'
    script:
    "scater.R $matrix $genes"
}

workflow {
    checkParams()

    // 2.1.0 cells: fastq files named after the cell id; <id>_1/<id>_2 (or _R1/_R2) make a paired-end cell
    cells = channel.fromPath("${params.cellhome}/**.{fastq,fq,fastq.gz,fq.gz}")
        .map { f -> [(f.name =~ /^(.+?)(?:_R?[12])?\.f(?:ast)?q(?:\.gz)?$/)[0][1], f] }
        .groupTuple(sort: true)
        .map { id, fs ->
            if (fs.size() > 2) error "Cell ${id}: expected 1 or 2 fastq files, found ${fs*.name}"
            [[id: id, single_end: fs.size() == 1], fs]
        }

    // 2.1.1
    if (params.readqc == 'fastp')        { FASTP(cells);   clean = FASTP.out.reads;   qc = FASTP.out.report }
    else if (params.readqc == 'afterqc') { AFTERQC(cells); clean = AFTERQC.out.reads; qc = AFTERQC.out.report }
    else                                 { FASTQC(cells);  clean = cells;             qc = FASTQC.out }

    GTF_TABLES(file(params.gtf, checkIfExists: true))
    gtf = GTF_TABLES.out.gtf
    def prebuilt = params["${params.quant}_index"]
    index = prebuilt ? channel.value(file(prebuilt, checkIfExists: true))
          : params.quant == 'salmon' ? SALMON_INDEX(file(params.transcript_fasta, checkIfExists: true))
          : params.quant == 'star'   ? STAR_INDEX(file(params.genome_fasta, checkIfExists: true), gtf)
          :                            HISAT2_INDEX(file(params.genome_fasta, checkIfExists: true))

    // 2.1.2 - 2.1.3
    def auto = params.strandness == 'auto'
    strand_txt = channel.empty()
    if (params.quant == 'salmon') {
        if (auto) strand_txt = INFER_STRAND(SALMON_PROBE(firstCell(clean), index))
        strand = auto ? strand_txt.map { strandFrom(it) } : channel.value(params.strandness)
        SALMON_QUANT(clean, index, strand)
        matrix = TXIMPORT(SALMON_QUANT.out.collect(), GTF_TABLES.out.tx2gene)
        genes = GTF_TABLES.out.genes
        aln_qc = SALMON_QUANT.out
        cnt_qc = channel.empty()
    } else {
        if (params.quant == 'star') { STAR_ALIGN(clean, index);   bam = STAR_ALIGN.out.bam;   aln_qc = STAR_ALIGN.out.log }
        else                        { HISAT2_ALIGN(clean, index); bam = HISAT2_ALIGN.out.bam; aln_qc = HISAT2_ALIGN.out.log }
        if (auto) strand_txt = INFER_STRAND(firstCell(bam).combine(GTF_TABLES.out.bed).map { _meta, b, bed -> [b, bed] })
        strand = auto ? strand_txt.map { strandFrom(it) }.first() : channel.value(params.strandness)   // first(): reusable value for every cell

        if (params.count == 'featurecounts') {
            FEATURECOUNTS(bam, gtf, strand); tables = FEATURECOUNTS.out.counts; cnt_qc = FEATURECOUNTS.out.summary
        } else if (params.count == 'htseq') {
            HTSEQ(bam, gtf, strand); tables = HTSEQ.out; cnt_qc = HTSEQ.out
        } else {
            STRINGTIE(bam, gtf, strand)
            STRINGTIE_MERGE(STRINGTIE.out.collect(), gtf)
            STRINGTIE_ABUND(bam, STRINGTIE_MERGE.out, strand)
            tables = STRINGTIE_ABUND.out; cnt_qc = channel.empty()
        }
        MERGE_TABLES(tables.collect())
        matrix = MERGE_TABLES.out.matrix
        genes = params.count == 'stringtie' ? MERGE_TABLES.out.genes : GTF_TABLES.out.genes
    }
    strand.subscribe { log.info "SIP strandness: $it" }

    MULTIQC(qc.mix(aln_qc, cnt_qc, strand_txt).collect())

    // 2.1.4
    results = matrix.mix(MULTIQC.out.data, strand_txt)
    if (params.using_seurat)  { SEURAT(matrix, genes);   results = results.mix(SEURAT.out) }
    if (params.using_monocle) { MONOCLE2(matrix, genes); results = results.mix(MONOCLE2.out) }
    if (params.using_scater)  { SCATER(matrix, genes);   results = results.mix(SCATER.out) }

    // 2.2 summary report
    def runName = workflow.runName
    REPORT(results.collect(), strand.map { s -> reportSettings(s, runName) })

    // 2.2 e-mail notification
    def email = params.email     // params is not visible inside the handler
    def outdir = params.outdir
    workflow.onComplete {
        if (email)
            sendMail(to: email, subject: "SIP ${workflow.success ? 'finished' : 'failed'}",
                     body: "Duration: ${workflow.duration}\nOutput: ${outdir}\n${workflow.errorMessage ?: ''}")
    }
}
