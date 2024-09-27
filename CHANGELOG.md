# nf-microbe/toolkit: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3beta - Testing grounds [2024-09-27]

Third beta release of nf-microbe/toolkit. Includes options for downloading SRA FastQ files,
Logan assemblies, and Logan sequences processed to be input into assemblers.

### `Added`

- [#21](https://github.com/nf-microbe/toolkit/pull/21) Added options for downloading
FastQ files from SRA (`sra/sratools` module added by @CarsonJM)
- [#21](https://github.com/nf-microbe/toolkit/pull/21) Added options for downloading
Logan assemblies (`logan/contigawscli` module added by @CarsonJM)
- [#21](https://github.com/nf-microbe/toolkit/pull/21) Added options for downloading
and processing Logan sequences for input into assembly (`logan/contigawsclimultiplier`
`logan/unitigawsclimultiplier`, `accession_logan_fasta` modules/subworkflows added by @CarsonJM)

### `Fixed`

- [#21](https://github.com/nf-microbe/toolkit/pull/21) Fixed UW Hyak profile to retry with exitStatus 250, an
out of memory error for SPAdes (Fixed by @CarsonJM)

### `Dependencies`

### `Deprecated`

## v0.2beta - Stomping grounds [2024-09-24]

Second beta release of nf-microbe/toolkit. Includes assembly extension, assembly QC,
MGE classification, and virus completeness estimation.

### `Added`

- [#17](https://github.com/nf-microbe/toolkit/pull/17) Added assembly extension, assembly QC,
MGE classification, and virus completeness (added by @CarsonJM)

### `Fixed`

- [#17](https://github.com/nf-microbe/toolkit/pull/17) Fixed resource requests for UW Hyak to prevent over requesting resources. (Fixed by @CarsonJM)

### `Dependencies`

### `Deprecated`

## v0.1beta - Testing grounds [2024-09-06]

Initial beta release of nf-microbe/toolkit for functionality testing,
created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#1](https://github.com/nf-microbe/toolkit/pull/2) Added read preprocessing functionality (added by @CarsonJM)
- [#2](https://github.com/nf-microbe/toolkit/pull/4) Added read assembly functionality (added by @CarsonJM)
- [#5](https://github.com/nf-microbe/toolkit/pull/5) Incorporated robust CI testing (added by @CarsonJM)

### `Fixed`

### `Dependencies`

### `Deprecated`
