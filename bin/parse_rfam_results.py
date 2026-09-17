#!/usr/bin/env python3
"""
parse_rfam_results.py — parse cmsearch tblout files, remove overlapping hits,
apply GA score thresholds (with rRNA overrides), and write a TSV of final hits.

Logic ported from ensembl-anno (ensembl_anno.py):
  - parse_rfam_tblout  (lines 1135-1170)
  - extract_rfam_metrics (parses .cm for NAME/ACC/CLEN/W/GA)
  - get_rfam_seed_descriptions (parses Rfam.seed)
  - remove_rfam_overlap (lines 1193-1239)
  - filter_rfam_results (lines 1242-1270)

Chunk coordinate encoding:
  Target sequence names in tblout are expected to carry offset information in
  the format  seqname.rs{chunk_start}.re{chunk_end}  (1-based, inclusive).
  Absolute genomic coordinates are computed by adding (chunk_start - 1) to
  each hit position.

Usage:
    parse_rfam_results.py \\
        --tblout_dir tblouts/ \\
        --rfam_cm rfam_filtered.cm \\
        --rfam_seed Rfam.seed \\
        --output rfam_hits.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# rRNA length overrides (from ensembl-analysis Perl pipeline)
# Models not listed here use the GA threshold from the .cm file.
# ---------------------------------------------------------------------------
RRNA_THRESHOLDS = {
    "LSU_rRNA_eukarya": 1700,
    "SSU_rRNA_eukarya": 1600,
    "5_8S_rRNA":        85,
    "5S_rRNA":          75,
}
# These bacterial/archaeal rRNA models are skipped entirely
RRNA_SKIP = {"LSU_rRNA_archaea", "LSU_rRNA_bacteria",
             "SSU_rRNA_archaea", "SSU_rRNA_bacteria",
             "SSU_rRNA_microsporidia"}

# Regex to decode chunk coordinates from target name
CHUNK_RE = re.compile(r"^(.*?)\.rs(\d+)\.re(\d+)$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tblout_dir", required=True,
                        help="Directory containing *.tblout files")
    parser.add_argument("--rfam_cm", required=True,
                        help="Filtered Rfam covariance model file")
    parser.add_argument("--rfam_seed", required=True,
                        help="Rfam.seed alignment file")
    parser.add_argument("--output", default="rfam_hits.tsv",
                        help="Output TSV path (default: rfam_hits.tsv)")
    parser.add_argument("--include_prokaryotic_rrna", action="store_true",
                        help="Include bacterial/archaeal rRNA models (for mgnify-assembly mode)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Parse .cm metrics
# ---------------------------------------------------------------------------
def extract_rfam_metrics(rfam_cm_path):
    """Return dict: accession -> {name, clen, w, ga} from a (filtered) .cm file.

    Each record in a .cm file consists of two sections separated by ``//``: the
    INFERNAL covariance model — which carries the ``GA`` gathering cutoff — and a
    HMMER3 filter model that repeats ``NAME``/``ACC`` but has *no* ``GA`` line.
    Only the covariance-model section is parsed; the filter section is skipped so
    that its GA-less metrics cannot overwrite the model's real cutoff (otherwise
    every non-rRNA hit is discarded downstream for lack of a GA threshold).
    """
    metrics = {}
    current = {}
    in_cm = False
    with open(rfam_cm_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("INFERNAL"):
                current = {}
                in_cm = True
            elif line.startswith("HMMER3"):
                # Start of the filter HMM section — stop recording until ``//``
                in_cm = False
            elif line.startswith("//"):
                if in_cm and "acc" in current:
                    metrics[current["acc"]] = dict(current)
                current = {}
                in_cm = False
            elif not in_cm:
                continue
            elif line.startswith("NAME"):
                current["name"] = line.split()[1]
            elif line.startswith("ACC"):
                current["acc"] = line.split()[1]
            elif line.startswith("CLEN"):
                current["clen"] = int(line.split()[1])
            elif line.startswith("W "):
                # ``W`` window length; the trailing space avoids matching ``WBETA``
                current["w"] = int(float(line.split()[1]))
            elif line.startswith("GA"):
                # GA  <bit_score>;  — the sequence-level gathering cutoff
                current["ga"] = float(line.split()[1].rstrip(";"))
    return metrics


# ---------------------------------------------------------------------------
# Parse Rfam.seed descriptions
# ---------------------------------------------------------------------------
def get_rfam_seed_descriptions(rfam_seed_path):
    """Return dict: accession -> {type, description}."""
    descriptions = {}
    current_acc = None
    current_type = "misc_RNA"
    current_desc = ""
    with open(rfam_seed_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GF AC"):
                current_acc = line.split()[-1].strip()
            elif line.startswith("#=GF TP"):
                current_type = line.split(None, 2)[-1].strip()
            elif line.startswith("#=GF DE"):
                current_desc = line.split(None, 2)[-1].strip()
            elif line.startswith("//"):
                if current_acc:
                    descriptions[current_acc] = {
                        "type": current_type,
                        "description": current_desc,
                    }
                current_acc = None
                current_type = "misc_RNA"
                current_desc = ""
    return descriptions


# ---------------------------------------------------------------------------
# Parse tblout files
# ---------------------------------------------------------------------------
def parse_tblout_file(path):
    """Yield hit dicts from one cmsearch --tblout file."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 17:
                continue
            target_name  = cols[0]
            query_name   = cols[2]   # model NAME
            query_acc    = cols[3]   # model ACC (may be '-')
            mdl_from     = int(cols[5])
            mdl_to       = int(cols[6])
            seq_from     = int(cols[7])
            seq_to       = int(cols[8])
            strand       = cols[9]   # '+' or '-'
            score        = float(cols[14])
            evalue       = float(cols[15])
            inc          = cols[16]  # '!' or '?'

            # Decode chunk coordinates from target name
            m = CHUNK_RE.match(target_name)
            if m:
                seqname     = m.group(1)
                chunk_start = int(m.group(2))  # 1-based
                offset      = chunk_start - 1
            else:
                seqname = target_name
                offset  = 0

            # Convert to 1-based absolute genomic coordinates
            abs_start = seq_from + offset
            abs_end   = seq_to   + offset

            # Ensure start < end regardless of strand
            if abs_start > abs_end:
                abs_start, abs_end = abs_end, abs_start

            yield {
                "seqname":    seqname,
                "start":      abs_start,
                "end":        abs_end,
                "strand":     strand,
                "score":      score,
                "evalue":     evalue,
                "query_name": query_name,
                "query_acc":  query_acc if query_acc != "-" else "",
                "mdl_from":   mdl_from,
                "mdl_to":     mdl_to,
                "inc":        inc,
            }


def parse_all_tblouts(tblout_dir):
    """Collect all hits from *.tblout files in a directory."""
    hits = []
    tblout_files = [
        os.path.join(tblout_dir, f)
        for f in os.listdir(tblout_dir)
        if f.endswith(".tblout")
    ]
    if not tblout_files:
        print(f"WARNING: no .tblout files found in {tblout_dir}", file=sys.stderr)
    for path in tblout_files:
        hits.extend(parse_tblout_file(path))
    return hits


# ---------------------------------------------------------------------------
# Remove overlapping hits (keep highest-scoring per locus)
# ---------------------------------------------------------------------------
def remove_rfam_overlap(hits):
    """
    For each sequence, greedily keep the highest-scoring non-overlapping hit.
    Hits on the same strand that overlap by > 0 bp are deduplicated.

    Port of remove_rfam_overlap (ensembl_anno.py lines 1193-1239).
    """
    # Group by (seqname, strand)
    by_locus = defaultdict(list)
    for hit in hits:
        key = (hit["seqname"], hit["strand"])
        by_locus[key].append(hit)

    kept = []
    for locus_hits in by_locus.values():
        # Sort by score descending
        locus_hits.sort(key=lambda h: h["score"], reverse=True)
        accepted = []
        for hit in locus_hits:
            overlaps = False
            for acc in accepted:
                # Overlap if intervals share any base
                if hit["start"] <= acc["end"] and hit["end"] >= acc["start"]:
                    overlaps = True
                    break
            if not overlaps:
                accepted.append(hit)
        kept.extend(accepted)
    return kept


# ---------------------------------------------------------------------------
# Filter by GA thresholds + rRNA overrides
# ---------------------------------------------------------------------------
def filter_rfam_results(hits, cm_metrics, seed_descriptions, include_prokaryotic_rrna=False):
    """
    Apply per-model GA score thresholds.  rRNA models use fixed length-based
    thresholds (RRNA_THRESHOLDS); some prokaryotic rRNA models are skipped.

    Port of filter_rfam_results (ensembl_anno.py lines 1242-1270).
    Returns list of hits with 'biotype' and 'accession' fields added.
    """
    final = []
    for hit in hits:
        name = hit["query_name"]

        # Skip prokaryotic/microsporidia rRNA models (unless in mgnify-assembly mode)
        if name in RRNA_SKIP and not include_prokaryotic_rrna:
            continue

        # Determine accession: prefer the one from the tblout, fall back to
        # looking up by model name in cm_metrics
        acc = hit.get("query_acc") or ""
        if not acc:
            for cm_acc, meta in cm_metrics.items():
                if meta.get("name") == name:
                    acc = cm_acc
                    break

        # Threshold: rRNA overrides take priority
        if name in RRNA_THRESHOLDS:
            threshold = RRNA_THRESHOLDS[name]
            hit_length = hit["end"] - hit["start"] + 1
            if hit_length < threshold:
                continue
        else:
            cm_info = cm_metrics.get(acc)
            if cm_info is None or "ga" not in cm_info:
                print(
                    f"WARNING: no GA threshold for model '{name}' (acc='{acc}') — skipping hit",
                    file=sys.stderr,
                )
                continue
            if hit["score"] < cm_info["ga"]:
                continue

        # Assign biotype from seed descriptions (or derive from model name)
        seed_info = seed_descriptions.get(acc, {})
        biotype = assign_biotype(name, acc, seed_info)

        hit_out = dict(hit)
        hit_out["accession"] = acc
        hit_out["biotype"]   = biotype
        final.append(hit_out)

    return final


# ---------------------------------------------------------------------------
# Biotype assignment
# Port of create_rfam_gtf biotype logic (ensembl_anno.py lines 1283-1399)
# ---------------------------------------------------------------------------
def assign_biotype(name, acc, seed_info):
    """Return an Ensembl-style biotype string for a given model."""
    seed_type = seed_info.get("type", "")

    # Check seed type string first
    if "snoRNA" in seed_type or "SNORD" in name or "SNORA" in name:
        if "scaRNA" in seed_type or "scaRNA" in name:
            return "scaRNA"
        return "snoRNA"
    if "snRNA" in seed_type or name.startswith("U") or "snRNA" in name:
        return "snRNA"
    if "rRNA" in name or "rRNA" in seed_type:
        if "LSU" in name or "28S" in name or "23S" in name:
            return "rRNA"
        if "SSU" in name or "18S" in name or "16S" in name:
            return "rRNA"
        if "5S" in name or "5_8S" in name:
            return "rRNA"
        return "rRNA"
    if "RNaseP" in name or "rnpB" in name:
        return "RNase_P_RNA"
    if "SRP" in name or "7SL" in name or "Metazoa_SRP" in name:
        return "SRP_RNA"
    if "Vault" in name or "vtRNA" in name:
        return "vault_RNA"
    if "Y_RNA" in name or "Ro" in name:
        return "Y_RNA"
    if "7SK" in name:
        return "misc_RNA"
    if "ribozyme" in seed_type.lower() or "Hammerhead" in name or "HDV" in name:
        return "ribozyme"
    if "antisense" in seed_type.lower():
        return "antisense_RNA"
    if "mir" in name.lower() or "miRNA" in seed_type:
        return "pre_miRNA"
    if "tRNA" in name or "tRNA" in seed_type:
        return "tRNA"
    if "lncRNA" in seed_type or "lncRNA" in name:
        return "lncRNA"
    # Default
    return "misc_RNA"


# ---------------------------------------------------------------------------
# Write output TSV
# ---------------------------------------------------------------------------
HEADER = [
    "seqname", "start", "end", "strand", "score", "evalue",
    "query_name", "accession", "biotype"
]


def write_tsv(hits, output_path):
    with open(output_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for hit in sorted(hits, key=lambda h: (h["seqname"], h["start"])):
            row = [
                hit["seqname"],
                str(hit["start"]),
                str(hit["end"]),
                hit["strand"],
                str(hit["score"]),
                str(hit["evalue"]),
                hit["query_name"],
                hit["accession"],
                hit["biotype"],
            ]
            out.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    print("parse_rfam_results: loading CM metrics...", file=sys.stderr)
    cm_metrics = extract_rfam_metrics(args.rfam_cm)
    print(f"  {len(cm_metrics)} models in CM file", file=sys.stderr)

    print("parse_rfam_results: loading Rfam.seed descriptions...", file=sys.stderr)
    seed_descriptions = get_rfam_seed_descriptions(args.rfam_seed)
    print(f"  {len(seed_descriptions)} entries in seed file", file=sys.stderr)

    print("parse_rfam_results: parsing tblout files...", file=sys.stderr)
    raw_hits = parse_all_tblouts(args.tblout_dir)
    print(f"  {len(raw_hits)} raw hits", file=sys.stderr)

    print("parse_rfam_results: removing overlaps...", file=sys.stderr)
    dedup_hits = remove_rfam_overlap(raw_hits)
    print(f"  {len(dedup_hits)} hits after overlap removal", file=sys.stderr)

    print("parse_rfam_results: applying GA thresholds...", file=sys.stderr)
    final_hits = filter_rfam_results(dedup_hits, cm_metrics, seed_descriptions,
                                     include_prokaryotic_rrna=args.include_prokaryotic_rrna)
    print(f"  {len(final_hits)} hits after filtering", file=sys.stderr)

    write_tsv(final_hits, args.output)
    print(f"parse_rfam_results: wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
