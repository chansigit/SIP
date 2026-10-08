#!/usr/bin/env python3
"""GTF -> tx2gene.tsv, genes.tsv (gene_id, gene_name), genes.bed12 (for RSeQC)."""

import gzip, re, sys
from collections import defaultdict

gtf = sys.argv[1]
attr = lambda s, k: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]
exons = defaultdict(list)  # tx -> [(chrom, start0, end, strand)]
tx2gene, names = {}, {}

for line in (gzip.open if gtf.endswith(".gz") else open)(gtf, "rt"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    g, t = attr(f[8], "gene_id"), attr(f[8], "transcript_id")
    if not t:
        continue
    v = attr(f[8], "transcript_version")  # Ensembl cDNA headers carry the version
    tx2gene[f"{t}.{v}" if v and not t.endswith("." + v) else t] = g
    names[g] = attr(f[8], "gene_name") or g
    exons[t].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))

with open("tx2gene.tsv", "w") as o:
    o.writelines(f"{t}\t{g}\n" for t, g in tx2gene.items())
with open("genes.tsv", "w") as o:
    o.writelines(f"{g}\t{n}\n" for g, n in names.items())
with open("genes.bed12", "w") as o:
    for t, ex in exons.items():
        ex.sort(key=lambda e: e[1])
        c, s, e, st = ex[0][0], ex[0][1], ex[-1][2], ex[0][3]
        o.write(
            "\t".join(
                map(
                    str,
                    [
                        c,
                        s,
                        e,
                        t,
                        0,
                        st,
                        s,
                        e,
                        0,
                        len(ex),
                        ",".join(str(x[2] - x[1]) for x in ex) + ",",
                        ",".join(str(x[1] - s) for x in ex) + ",",
                    ],
                )
            )
            + "\n"
        )
