#!/usr/bin/env bash
# Fetch the paper's example data (GSE38495 subset: 8 hESC + 19 LNCaP, SMART-seq, single-end)
# and an Ensembl reference. Usage: download_gse38495.sh <outdir>
set -euo pipefail
out=${1:?usage: $0 <outdir>}
here=$(cd "$(dirname "$0")" && pwd)
rel=115
ens=https://ftp.ensembl.org/pub/release-$rel

mkdir -p "$out/cells" "$out/ref"
cut -f1,4 "$here/gse38495_cells.tsv" | while read -r run url; do
    echo "$out/cells/$run.fastq.gz https://$url"
done | xargs -P 4 -n 2 sh -c '[ -s "$0" ] && exit 0; curl -sSfL --retry 5 -o "$0.part" "$1" && mv "$0.part" "$0"'

cd "$out/ref"
for f in fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
         gtf/homo_sapiens/Homo_sapiens.GRCh38.$rel.gtf.gz \
         fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz; do
    [ -s "$(basename $f)" ] || curl -sSfL --retry 5 -O "$ens/$f"
done
ls -la "$out/cells" | tail -3; ls -la "$out/ref"
