#!/usr/bin/env bash
# Tiny test set from the paper data: chr22 reference + 4 single-end cells (200k reads)
# + 1 paired-end cell (mate 2 = reverse complement of mate 1) to exercise the paired-end paths.
# Usage: make_testdata.sh <gse38495 dir from paper/download_gse38495.sh> <outdir>
set -eu   # no pipefail: `zcat | head` ends with SIGPIPE
src=${1:?usage: $0 <gse38495dir> <outdir>}; out=${2:?}
mkdir -p "$out/cells" "$out/ref"

zcat "$src/ref/Homo_sapiens.GRCh38.115.gtf.gz" | awk '$1=="22"' | gzip > "$out/ref/chr22.gtf.gz"
zcat "$src/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
    | awk '/^>/{keep = /chromosome:GRCh38:22:/} keep' | gzip > "$out/ref/chr22.cdna.fa.gz"
cp "$src/ref/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz" "$out/ref/chr22.fa.gz"

for run in SRR522081 SRR522082 SRR522089 SRR522090; do
    zcat "$src/cells/$run.fastq.gz" | head -n 800000 | gzip > "$out/cells/$run.fastq.gz"
done
zcat "$src/cells/SRR522083.fastq.gz" | head -n 800000 > /tmp/pe_$$.fq
awk 'NR%4==1{sub(/ .*/,""); print $0"/1"; next} NR%4==3{print "+"; next} 1' /tmp/pe_$$.fq | gzip > "$out/cells/PE522083_1.fastq.gz"
awk 'NR%4==1{sub(/ .*/,""); print $0"/2"} NR%4==2{print rev(comp($0))} NR%4==3{print "+"} NR%4==0{print rev($0)}
     function rev(s,  r,i){r=""; for(i=length(s);i>0;i--) r=r substr(s,i,1); return r}
     function comp(s){gsub(/A/,"t",s); gsub(/T/,"a",s); gsub(/C/,"g",s); gsub(/G/,"c",s); return toupper(s)}' \
    /tmp/pe_$$.fq | gzip > "$out/cells/PE522083_2.fastq.gz"
rm -f /tmp/pe_$$.fq
ls -la "$out"/*
