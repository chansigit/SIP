#!/usr/bin/env bash
#SBATCH --job-name=sip-reproduce
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
# Reproduce every SIP result from scratch and check it against the recorded md5 sums:
#   the paper example (GSE38495) and one run per interchangeable tool on the chr22 test set.
# Usage (from the SIP directory):  sbatch reproducible/run.sh [outdir]
#                         or, inside an allocation:  bash reproducible/run.sh [outdir]
# Set SIP_DATA to reuse an existing download; it is md5-checked either way.
set -euo pipefail
sip=$(cd "$(dirname "$0")/.." && pwd)
[ -f "$sip/main.nf" ] || sip=${SLURM_SUBMIT_DIR:-}      # sbatch runs a spooled copy of this script
[ -f "$sip/main.nf" ] || { echo "run this from the SIP directory" >&2; exit 1; }
here=$sip/reproducible
out=$(realpath -m "${1:-$SCRATCH/sip_reproduce}")
data=${SIP_DATA:-$out/data}
mkdir -p "$out/results"; cd "$out"

export NXF_VER=25.04.7                                  # exact Nextflow release
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-$out/singularity}
if type ml &>/dev/null; then ml load nextflow/25.04.7; fi

# 1. inputs: download, verify, re-fetch anything that does not match
check() { (cd "$data" && md5sum -c --quiet "$here/inputs.md5"); }
"$sip/paper/download_gse38495.sh" "$data" >/dev/null
if ! check; then
    (cd "$data" && md5sum -c "$here/inputs.md5" 2>/dev/null | awk -F': ' '$2 != "OK" {print $1}' | xargs -r rm -f)
    "$sip/paper/download_gse38495.sh" "$data" >/dev/null
    check
fi
"$sip/tests/make_testdata.sh" "$data" "$out/testdata" >/dev/null

# 2. runs (fresh work directory, no -resume)
cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}
mem=$(( ${SLURM_MEM_PER_NODE:-$(( $(free -m | awk '/^Mem:/ {print $2}') * 8 / 10 ))} / 1024 ))
cat > local.config <<EOF
executor { cpus = $cpus; memory = ${mem}.GB }
process.resourceLimits = [cpus: $cpus, memory: ${mem}.GB, time: 24.h]
EOF
nf() { nextflow run "$sip" -profile singularity -c local.config -w work -ansi-log false "$@"; }
R=$data/ref; T=$out/testdata
nf -params-file "$sip/paper/gse38495.yaml" --cellhome "$data/cells" --using_monocle --using_scater \
   --gtf "$R/Homo_sapiens.GRCh38.115.gtf.gz" --transcript_fasta "$R/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
   --outdir results/paper
test_run() { nf --cellhome "$T/cells" --gtf "$T/ref/chr22.gtf.gz" --min_features 10 "$@"; }
test_run --transcript_fasta "$T/ref/chr22.cdna.fa.gz" --using_seurat --using_monocle --using_scater --outdir results/salmon
test_run --genome_fasta "$T/ref/chr22.fa.gz" --readqc fastqc --quant hisat2 --count featurecounts --using_seurat --outdir results/hisat2_featurecounts
test_run --genome_fasta "$T/ref/chr22.fa.gz" --readqc afterqc --quant star --count htseq --using_scater --outdir results/star_htseq
test_run --genome_fasta "$T/ref/chr22.fa.gz" --readqc fastqc --quant hisat2 --count stringtie --quant_unit tpm \
         --using_seurat --using_monocle --outdir results/hisat2_stringtie

# 3. provenance and output check
{ echo "date: $(date -Is)"; echo "host: $(hostname)"; nextflow -version 2>&1 | grep -m1 version
  singularity --version; echo "git: $(git -C "$sip" describe --always --dirty 2>/dev/null || echo none)"
  (cd "$sip" && sha256sum main.nf nextflow.config bin/* paper/* tests/* reproducible/*.md5); } > results/environment.txt
(cd results && md5sum -c "$here/expected.md5")
echo "all results match reproducible/expected.md5"
