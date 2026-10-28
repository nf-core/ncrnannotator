# nf-core/ncrnannotator: Output

## Introduction

This document describes the output produced by ncrnannotator. All paths are relative to the `--outdir` directory specified at runtime.

## Pipeline overview

ncrnannotator runs the following steps:

1. **FILTER_RFAM_CM** — filter Rfam covariance models to clade-specific accessions (skipped in `mgnify-assembly` and `full` modes)
2. **GENOME_CHUNK** — split genome into non-overlapping windows for parallel processing
3. **CMSEARCH** — search each chunk against Rfam covariance models using Infernal
4. **PARSE_RFAM** — consolidate results, remove overlapping hits, apply GA score thresholds
5. **RFAM_TO_FORMATS** — convert hits to GTF, GFF3, and BED annotation files
6. **MultiQC** — aggregate software versions and pipeline metrics into a report

## Output directories

### `rfam/`

Intermediate Rfam files published for inspection and reuse.

<details markdown="1">
<summary>Output files</summary>

- `rfam/`
  - `rfam_filtered.cm` — Rfam covariance model file filtered to clade-specific accessions (not present in `mgnify-assembly` or `full` mode)
  - `rfam_hits.tsv` — Tab-separated table of all final ncRNA hits after overlap removal and score filtering

</details>

#### `rfam_hits.tsv` format

| Column       | Description                                                      |
| ------------ | ---------------------------------------------------------------- |
| `seqname`    | Sequence name from the input FASTA                               |
| `start`      | Hit start position (1-based)                                     |
| `end`        | Hit end position (1-based, inclusive)                            |
| `strand`     | Strand (`+` or `-`)                                              |
| `score`      | Infernal bit score                                               |
| `evalue`     | E-value                                                          |
| `query_name` | Rfam model name (e.g. `U1`, `5S_rRNA`)                           |
| `accession`  | Rfam accession (e.g. `RF00003`)                                  |
| `biotype`    | Ensembl-style biotype (e.g. `snRNA`, `snoRNA`, `rRNA`, `lncRNA`) |

### `annotation/`

Final annotation files ready for use in genome browsers and downstream analysis.

<details markdown="1">
<summary>Output files</summary>

- `annotation/`
  - `annotation.gtf` — GTF v2.2 annotation with gene/transcript/exon features
  - `annotation.gff3` — GFF3 annotation with gene/ncRNA/exon hierarchy
  - `annotation.bed` — 6-column BED file (0-based start coordinates)

</details>

#### GTF format

Ensembl-style GTF with three feature types per hit: `gene`, `transcript`, `exon`. Key attributes:

- `gene_id` — sequential identifier (e.g. `rfam_gene_000001`)
- `gene_name` — Rfam model name
- `gene_biotype` — Ensembl biotype (e.g. `snRNA`, `snoRNA`, `rRNA`)
- `rfam_accession` — Rfam family accession

#### GFF3 format

Standard GFF3 with `gene` → `ncRNA` (or biotype-specific type) → `exon` hierarchy. Feature types follow Sequence Ontology conventions.

#### BED format

6-column BED (chrom, chromStart, chromEnd, name, score, strand). Start is 0-based. Score is clamped to 0–1000.

### `multiqc/`

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html` — standalone HTML report viewable in any browser
  - `multiqc_data/` — parsed statistics from all pipeline steps
  - `multiqc_plots/` — static images from the report

</details>

[MultiQC](http://multiqc.info) aggregates software versions and pipeline metrics into a single report.

### `pipeline_info/`

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html` — Nextflow execution report (resource usage per task)
  - `execution_timeline.html` — timeline of all tasks
  - `execution_trace.txt` — tab-separated trace of all tasks
  - `nf_core_ncrnannotator_software_mqc_versions.yml` — software versions used

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Biotype reference

The following biotypes are assigned based on Rfam model names and seed classifications:

| Biotype                                      | Examples                                               |
| -------------------------------------------- | ------------------------------------------------------ |
| `snRNA`                                      | U1, U2, U4, U5, U6, U11, U12                           |
| `snoRNA`                                     | SNORD, SNORA families                                  |
| `scaRNA`                                     | scaRNA families                                        |
| `rRNA`                                       | 5S_rRNA, 5_8S_rRNA, SSU_rRNA_eukarya, LSU_rRNA_eukarya |
| `rRNA` (prokaryotic, `mgnify-assembly` only) | SSU_rRNA_bacteria, LSU_rRNA_archaea, etc.              |
| `tRNA`                                       | tRNA families                                          |
| `pre_miRNA`                                  | miRNA precursors                                       |
| `lncRNA`                                     | Long non-coding RNA                                    |
| `SRP_RNA`                                    | Signal recognition particle RNA                        |
| `RNase_P_RNA`                                | RNase P families                                       |
| `vault_RNA`                                  | Vault RNA                                              |
| `Y_RNA`                                      | Y RNA families                                         |
| `ribozyme`                                   | Ribozyme families                                      |
| `antisense_RNA`                              | Antisense RNA                                          |
| `ncRNA`                                      | Other non-coding RNA (default)                         |
