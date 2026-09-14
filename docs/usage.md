# nf-core/ncrnannotator: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/ncrnannotator/usage](https://nf-co.re/ncrnannotator/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

ncrnannotator annotates non-coding RNA in genome assemblies using [Infernal](http://eddylab.org/infernal/) and the [Rfam](https://rfam.org/) database. It produces annotation files in GTF, GFF3, and BED format.

## Quick start

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

## Required inputs

| Parameter     | Description                                             |
| ------------- | ------------------------------------------------------- |
| `--fasta`     | Genome assembly in FASTA format (uncompressed or `.gz`) |
| `--mode`      | Annotation mode (see below)                             |
| `--rfam_cm`   | Rfam covariance model file (`Rfam.cm`)                  |
| `--rfam_seed` | Rfam seed alignment file (`Rfam.seed`)                  |
| `--outdir`    | Output directory                                        |

### Obtaining Rfam files

Download the latest Rfam release from the [Rfam FTP](https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/):

```bash
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz
gunzip Rfam.cm.gz Rfam.seed.gz
```

## Annotation modes

### `ensembl-vertebrates`

Annotates vertebrate genomes using a curated subset of Rfam families (snRNA, snoRNA, rRNA, SRP RNA, Y RNA, RNase P, Vault RNA). Uses 1 Mbp genome chunks.

```bash
nextflow run nf-core/ncrnannotator \
  --fasta Homo_sapiens.GRCh38.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

### `ensembl-invertebrates`

Annotates invertebrate genomes using a curated subset of Rfam families. Uses 100 kbp genome chunks for higher sensitivity on compact genomes.

```bash
nextflow run nf-core/ncrnannotator \
  --fasta Caenorhabditis_elegans.WBcel235.fa \
  --mode ensembl-invertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

### `mgnify-assembly`

Annotates metagenomic assemblies using the full Rfam database (no clade filtering). Includes prokaryotic rRNA (bacterial, archaeal). Uses 50 Mbp genome chunks.

```bash
nextflow run nf-core/ncrnannotator \
  --fasta metagenome_assembly.fa \
  --mode mgnify-assembly \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

### `full`

Annotates any genome using the complete Rfam database (no clade filtering) with standard Infernal covariance model scoring. Eukaryotic output only (prokaryotic rRNA excluded). Uses 1 Mbp genome chunks.

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode full \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

## Optional parameters

| Parameter           | Description                                          | Default          |
| ------------------- | ---------------------------------------------------- | ---------------- |
| `--chunk_size`      | Override genome chunk size (bp)                      | Mode-dependent   |
| `--rfam_accessions` | Custom Rfam accession list (overrides bundled lists) | Bundled per mode |

### Custom chunk size

For very large genomes or limited memory, override the chunk size:

```bash
# Smaller chunks use less memory per cmsearch job
--chunk_size 500000
```

## Running on a laptop

When running locally, use a custom config to cap resource usage:

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode ensembl-invertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -c conf/local.config \
  -profile conda
```

The bundled `conf/local.config` limits memory to 20 GB and 8 CPUs. Edit it to match your machine.

## Using a params file

Rather than specifying all parameters on the command line, you can use a YAML params file:

```bash
nextflow run nf-core/ncrnannotator -params-file params.yaml -profile docker
```

An example params file is provided at [`assets/params_example.yaml`](../assets/params_example.yaml).

> [!WARNING]
> Do not use `-c <file>` to specify pipeline parameters — only use `-c` for resource tuning. Use `-params-file` for parameter overrides.

## Resuming a run

Add `-resume` to restart from the last successful step:

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker \
  -resume
```

## Running on an HPC cluster

### SLURM example

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile singularity \
  -resume
```

Run Nextflow itself in a SLURM job to avoid timeouts on the head node:

```bash
#!/bin/bash
#SBATCH --job-name=ncrnannotator
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

nextflow run nf-core/ncrnannotator \
  -params-file params.yaml \
  -profile singularity \
  -resume
```

## Reproducibility

Pin the pipeline version with `-r`:

```bash
nextflow run nf-core/ncrnannotator -r 1.0.0 \
  --fasta genome.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

## Core Nextflow arguments

### `-profile`

Select a software environment profile:

- `docker` — Docker containers (recommended)
- `singularity` — Singularity containers (recommended for HPC)
- `conda` — Conda environments (fallback when containers unavailable)
- `podman`, `apptainer`, `charliecloud` — Alternative container runtimes
- `test` — Minimal test run (no input files needed)

Example with multiple profiles:

```bash
-profile test,docker
```

### `-resume`

Restart from cached intermediate results. See the [Nextflow resume docs](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

### `-c`

Provide a custom config for resource tuning only:

```bash
-c conf/local.config
```

## Nextflow memory requirements

Limit Nextflow's own JVM memory usage (add to `~/.bashrc`):

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
