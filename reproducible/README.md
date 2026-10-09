# Reproducing SIP results

`run.sh` rebuilds every result from public data and checks it bit-for-bit:

1. downloads GSE38495 (27 SMART-seq runs from ENA) and Ensembl 115, and checks them against `inputs.md5`
   (ENA's own md5 for the fastq files); a file that fails is fetched again, a second failure stops the run;
2. builds the chr22 test set (`tests/make_testdata.sh`);
3. runs the paper example (`paper/gse38495.yaml`) and one test run per interchangeable tool, each in a fresh
   work directory without `-resume`;
4. writes `results/environment.txt` (Nextflow and Singularity versions, sha256 of every pipeline file) and checks
   the deterministic outputs against `expected.md5`.

```bash
cd SIP
sbatch reproducible/run.sh /path/to/outdir                     # or: bash reproducible/run.sh inside an allocation
SIP_DATA=/path/to/gse38495 sbatch reproducible/run.sh outdir   # reuse an existing download (still md5-checked)
```

What is fixed:

| | |
|---|---|
| Nextflow | 25.04.7 (`NXF_VER`) |
| Software | every container pinned by sha256 digest in `nextflow.config` |
| Input data | `inputs.md5` |
| Parameters | `paper/gse38495.yaml` and the command lines in `run.sh` |
| Randomness | `set.seed(1)` in the R scripts; salmon 2.x quantification is deterministic; the strandness probe uses the alphabetically first cell |

`expected.md5` lists the 188 outputs that came out identical in independent runs on different nodes: salmon
quantifications, fastp reports, expression matrices, cluster assignments, marker tables, pseudotime, cell QC
metrics, the StringTie assembly and the RSeQC results. PDFs and HTML reports embed creation
dates and are not compared.
