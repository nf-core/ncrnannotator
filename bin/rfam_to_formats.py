#!/usr/bin/env python3
"""
rfam_to_formats.py — convert rfam_hits.tsv to GTF, GFF3, and BED annotation files.

Biotype logic ported from ensembl-anno (ensembl_anno.py, create_rfam_gtf, lines 1283-1399).

Output:
  annotation.gtf   — Ensembl-style GTF with gene/transcript/exon features
  annotation.gff3  — GFF3 with gene/mRNA/exon hierarchy
  annotation.bed   — 6-column BED (0-based start)

Usage:
    rfam_to_formats.py --hits rfam_hits.tsv --output_prefix annotation
"""

import argparse
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True,
                        help="TSV file produced by parse_rfam_results.py")
    parser.add_argument("--output_prefix", default="annotation",
                        help="Prefix for output files (default: annotation)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Read hits TSV
# ---------------------------------------------------------------------------
def load_hits(hits_path):
    hits = []
    with open(hits_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            record = dict(zip(header, parts))
            record["start"] = int(record["start"])
            record["end"]   = int(record["end"])
            record["score"] = float(record["score"])
            hits.append(record)
    return hits


# ---------------------------------------------------------------------------
# ID generation helpers
# ---------------------------------------------------------------------------
def make_gene_id(hit, idx):
    return f"rfam_gene_{idx:07d}"


def make_transcript_id(gene_id):
    return f"{gene_id}.1"


# ---------------------------------------------------------------------------
# GTF output
# ---------------------------------------------------------------------------
GTF_SOURCE = "Rfam"


def _gtf_attrs(**kwargs):
    """Format attribute pairs as a GTF attribute string."""
    parts = []
    for key, val in kwargs.items():
        parts.append(f'{key} "{val}"')
    return "; ".join(parts) + ";"


def write_gtf(hits, path):
    with open(path, "w") as out:
        out.write('##gtf-version 2.2\n')
        for idx, hit in enumerate(hits, start=1):
            chrom     = hit["seqname"]
            start_1   = hit["start"]       # 1-based inclusive
            end_1     = hit["end"]         # 1-based inclusive
            strand    = hit["strand"]
            score     = hit["score"]
            biotype   = hit["biotype"]
            acc       = hit["accession"]
            name      = hit["query_name"]

            gene_id   = make_gene_id(hit, idx)
            tx_id     = make_transcript_id(gene_id)

            shared = _gtf_attrs(
                gene_id=gene_id,
                transcript_id=tx_id,
                gene_name=name,
                gene_biotype=biotype,
                transcript_biotype=biotype,
                rfam_accession=acc,
                score=f"{score:.3f}",
            )

            # gene feature
            out.write("\t".join([
                chrom, GTF_SOURCE, "gene",
                str(start_1), str(end_1),
                ".", strand, ".",
                _gtf_attrs(gene_id=gene_id, gene_name=name,
                           gene_biotype=biotype, rfam_accession=acc)
            ]) + "\n")

            # transcript feature
            out.write("\t".join([
                chrom, GTF_SOURCE, "transcript",
                str(start_1), str(end_1),
                ".", strand, ".",
                _gtf_attrs(gene_id=gene_id, transcript_id=tx_id,
                           gene_name=name, gene_biotype=biotype,
                           transcript_biotype=biotype, rfam_accession=acc)
            ]) + "\n")

            # exon feature
            out.write("\t".join([
                chrom, GTF_SOURCE, "exon",
                str(start_1), str(end_1),
                ".", strand, ".",
                shared
            ]) + "\n")


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------
GFF3_SOURCE = "Rfam"


def _gff3_attrs(**kwargs):
    return ";".join(f"{k}={v}" for k, v in kwargs.items())


def write_gff3(hits, path):
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for idx, hit in enumerate(hits, start=1):
            chrom     = hit["seqname"]
            start_1   = hit["start"]
            end_1     = hit["end"]
            strand    = hit["strand"]
            score_str = f"{hit['score']:.3f}"
            biotype   = hit["biotype"]
            acc       = hit["accession"]
            name      = hit["query_name"]

            gene_id = make_gene_id(hit, idx)
            tx_id   = make_transcript_id(gene_id)
            exon_id = f"{tx_id}.exon1"

            # gene
            out.write("\t".join([
                chrom, GFF3_SOURCE, "gene",
                str(start_1), str(end_1),
                score_str, strand, ".",
                _gff3_attrs(ID=gene_id, Name=name,
                            biotype=biotype, rfam_accession=acc)
            ]) + "\n")

            # mRNA / ncRNA
            feature_type = biotype_to_gff3_type(biotype)
            out.write("\t".join([
                chrom, GFF3_SOURCE, feature_type,
                str(start_1), str(end_1),
                score_str, strand, ".",
                _gff3_attrs(ID=tx_id, Parent=gene_id, Name=name,
                            biotype=biotype, rfam_accession=acc)
            ]) + "\n")

            # exon
            out.write("\t".join([
                chrom, GFF3_SOURCE, "exon",
                str(start_1), str(end_1),
                ".", strand, ".",
                _gff3_attrs(ID=exon_id, Parent=tx_id)
            ]) + "\n")


def biotype_to_gff3_type(biotype):
    """Map biotype to a GFF3 SO term."""
    mapping = {
        "snRNA":         "snRNA",
        "snoRNA":        "snoRNA",
        "scaRNA":        "snoRNA",
        "rRNA":          "rRNA",
        "tRNA":          "tRNA",
        "pre_miRNA":     "pre_miRNA",
        "lncRNA":        "lnc_RNA",
        "RNase_P_RNA":   "RNase_P_RNA",
        "SRP_RNA":       "SRP_RNA",
        "vault_RNA":     "vault_RNA",
        "Y_RNA":         "Y_RNA",
        "ribozyme":      "ribozyme",
        "antisense_RNA": "antisense_RNA",
    }
    return mapping.get(biotype, "ncRNA")


# ---------------------------------------------------------------------------
# BED output (6-column, 0-based start)
# ---------------------------------------------------------------------------
def write_bed(hits, path):
    with open(path, "w") as out:
        for idx, hit in enumerate(hits, start=1):
            chrom     = hit["seqname"]
            bed_start = hit["start"] - 1   # convert to 0-based
            bed_end   = hit["end"]
            name      = hit["accession"] or hit["query_name"]
            score_int = min(1000, max(0, int(hit["score"])))
            strand    = hit["strand"]

            out.write("\t".join([
                chrom,
                str(bed_start),
                str(bed_end),
                name,
                str(score_int),
                strand
            ]) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    print(f"rfam_to_formats: loading {args.hits}...", file=sys.stderr)
    hits = load_hits(args.hits)
    print(f"  {len(hits)} hits loaded", file=sys.stderr)

    # Sort by chromosome then start position
    hits.sort(key=lambda h: (h["seqname"], h["start"]))

    gtf_path  = f"{args.output_prefix}.gtf"
    gff3_path = f"{args.output_prefix}.gff3"
    bed_path  = f"{args.output_prefix}.bed"

    write_gtf(hits, gtf_path)
    print(f"rfam_to_formats: wrote {gtf_path}", file=sys.stderr)

    write_gff3(hits, gff3_path)
    print(f"rfam_to_formats: wrote {gff3_path}", file=sys.stderr)

    write_bed(hits, bed_path)
    print(f"rfam_to_formats: wrote {bed_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
