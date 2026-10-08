# SIP: a Single-cell Interchangeable Pipeline

A Nextflow (DSL2) implementation of the pipeline described in
Chen et al., *SIP: An Interchangeable Pipeline for scRNA-seq Data Processing*, bioRxiv 10.1101/456772 (2018).
Tool versions are current, not the 2018 ones.

| Paper section | Step | Options |
|---|---|---|
| 2.1.0 | Cells: fastq files named after cell ids, grouped by name; layout detected from the group | `--cellhome` |
| 2.1.1 | Read QC & preprocessing, summarised by MultiQC | `--readqc fastp` (default) \| `fastqc` \| `afterqc` |
| 2.1.2 | Strandness detection with RSeQC, used to set the downstream parameters | `--strandness auto` (default) \| `unstranded` \| `forward` \| `reverse` |
| 2.1.2 | (Quasi-)alignment | `--quant salmon` (default) \| `star` \| `hisat2` |
| 2.1.2–2.1.3 | Gene-level counting (alignment modes); StringTie assembles, merges, re-estimates | `--count featurecounts` (default) \| `htseq` \| `stringtie` |
| 2.1.3 | Units; salmon results merged to genes with tximport | `--quant_unit count` \| `tpm` (salmon), `tpm` \| `fpkm` (StringTie) |
| 2.1.4 | Downstream: Seurat (clusters, t-SNE, DEGs), Monocle 2 (trajectory), scater (QC) | `--using_seurat`, `--using_monocle`, `--using_scater` |
| 2.2 | HPC, containers, resume, e-mail | `-profile slurm`, `-profile singularity\|docker`, `-resume`, `--email` |

Seurat parameters (`--min_features`, `--max_mt`, `--n_pcs`, `--tsne_dims`, `--tsne_perplexity`, `--cluster_dims`,
`--resolution`, `--de_min_pct`, `--heatmap_genes`, `--regress_vars`) are listed in `nextflow.config`.

## Input

`--cellhome` holds one fastq per single-end cell (`<cell>.fastq.gz`) or a pair per paired-end cell
(`<cell>_1.fastq.gz` / `<cell>_2.fastq.gz`, `_R1`/`_R2` also accepted), in any sub-folder layout.
References: `--gtf` always; `--transcript_fasta` for salmon, `--genome_fasta` for STAR/HISAT2
(or a prebuilt `--salmon_index` / `--star_index` / `--hisat2_index`). Ensembl FASTA/GTF are expected.

## Run

```bash
nextflow run main.nf -profile singularity -resume \
    --cellhome cells/ --gtf Homo_sapiens.GRCh38.115.gtf.gz \
    --transcript_fasta Homo_sapiens.GRCh38.cdna.all.fa.gz \
    --readqc fastp --quant salmon --using_seurat
```

Add `-profile singularity,slurm` to send every step to Slurm.

## Output (`--outdir`, default `sip_results/`)

`sip_report.pdf` (summary of the run: settings, QC tables, cluster/marker tables and all figures, vector graphics),
`readqc/`, `strandness/`, `quant/` (salmon), `assembly/` (StringTie), `matrix/gene_<unit>.tsv` (genes × cells),
`multiqc_report.html`, `seurat/`, `monocle2/`, `scater/`.

## Reproducing the paper's example (section 3)

```bash
paper/download_gse38495.sh /path/gse38495     # 8 hESC + 19 LNCaP, ~23 GB, plus Ensembl 115
nextflow run main.nf -profile singularity -params-file paper/gse38495.yaml \
    --cellhome /path/gse38495/cells --gtf /path/gse38495/ref/Homo_sapiens.GRCh38.115.gtf.gz \
    --transcript_fasta /path/gse38495/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

`paper/gse38495_cells.tsv` gives the cell type of each run.

`reproducible/run.sh` does all of this from scratch (plus the test runs below) with md5-checked inputs,
digest-pinned containers and a bit-for-bit check of the results; see `reproducible/README.md`.

## Tests

`tests/make_testdata.sh` builds a chr22 reference and 5 small cells (one paired-end) from the paper data;
run each `--readqc` / `--quant` / `--count` choice on it.

## Notes

- Strandness is inferred from one cell. In salmon mode its first 200k reads are mapped with salmon, and
  `infer_experiment.py` runs on that SAM against a transcript-coordinate BED; in alignment modes it runs
  on the first cell's BAM against the GTF.
- Monocle is pinned to 2.34: the 2.38 container ships dplyr 1.2, where Monocle's `group_by_()` is defunct.
- AfterQC is Python 2 and slow on full-size cells.
