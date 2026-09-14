<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-ncrnannotator_logo_dark.png">
    <img alt="nf-core/ncrnannotator" src="docs/images/nf-core-ncrnannotator_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/ncrnannotator)
[![GitHub Actions CI Status](https://github.com/nf-core/ncrnannotator/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/ncrnannotator/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/ncrnannotator/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/ncrnannotator/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/ncrnannotator/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/ncrnannotator)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ncrnannotator-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/ncrnannotator)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/ncrnannotator** is a bioinformatics pipeline for genome-level non-coding RNA (ncRNA) annotation. It takes a genome assembly in FASTA format and annotates ncRNA loci using Rfam covariance models searched with Infernal (cmsearch). The pipeline is taxonomically aware, supporting vertebrate and invertebrate annotation modes, and produces output in GTF, GFF3, and BED formats.

1. Filter Rfam covariance models to the target clade ([`filter_rfam_cm`](bin/filter_rfam_cm.py))
2. Chunk genome into windows for parallel processing ([`GENOME_CHUNK`](modules/local/genome_chunk/main.nf))
3. Search each chunk against filtered Rfam models ([`Infernal cmsearch`](http://eddylab.org/infernal/))
4. Parse hits, remove overlaps, apply GA score thresholds ([`parse_rfam_results`](bin/parse_rfam_results.py))
5. Export annotation in GTF, GFF3, and BED formats ([`rfam_to_formats`](bin/rfam_to_formats.py))
6. Present pipeline metrics ([`MultiQC`](http://multiqc.info/))

![ncrnannotator workflow diagram](docs/images/ncrnannotator.svg)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

ncrnannotator takes a genome assembly in FASTA format and Rfam database files as input. No samplesheet is required.

```bash
nextflow run nf-core/ncrnannotator \
  --fasta genome.fa \
  --mode ensembl-vertebrates \
  --rfam_cm Rfam.cm \
  --rfam_seed Rfam.seed \
  --outdir results \
  -profile docker
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/ncrnannotator/usage) and the [parameter documentation](https://nf-co.re/ncrnannotator/parameters).
For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/ncrnannotator/usage) and the [parameter documentation](https://nf-co.re/ncrnannotator/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/ncrnannotator/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/ncrnannotator/output).

## Credits

The first version of nf-core/ncrnannotator was written by Pedro Madrigal and Victoria Begley, members of the RNA Resources team at EMBL-EBI.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#ncrnannotator` channel](https://nfcore.slack.com/channels/ncrnannotator) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
