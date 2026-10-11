# nf-core/ncrnannotator: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/ncrnannotator, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

- Fixed `extract_rfam_metrics` in `parse_rfam_results.py` reading the GA gathering cutoff from the wrong section of each `.cm` record. The HMMER3 filter section (which has no `GA` line) was overwriting the covariance model's metrics, leaving every model without a GA threshold and causing all non-rRNA hits (tRNA, snRNA, snoRNA, SRP, etc.) to be silently discarded. The parser now only reads the INFERNAL covariance-model section.

### `Dependencies`

### `Deprecated`
