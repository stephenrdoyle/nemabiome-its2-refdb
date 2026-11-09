# sources/

Downloaded source data. Nothing here is committed except this file; each
build's `manifest.json` records which backbone it used.

| path | from |
|---|---|
| `nemabiome_its2_20260317/` | **The backbone.** The live nemabiome ITS2 database, version 20260317 (17 March 2026; 27,343 sequences, 15,108 non-redundant), built with markerDB. Only downloadable through its web app, https://ucvm-bio.shinyapps.io/Nematode_ITS2 (embedded at https://www.nemabiome.ca/its2-database): show every column including `sequence`, tick `excel` (and `dada2` as a cross-check), Download. The app returns `database.zip`; `raw/` holds it unzipped, and `bin/nemabiome_export.py raw/database.xlsx Nematode_ITS2_20260317` writes the `_seqs_nr.fasta` / `_taxonomy_nr.txt` pair. `checksums.txt`: SHA-256 of the zip and its contents as downloaded 2026-11-09. CC BY 4.0. |
| `nemabiome_its2_1.0.0/` | Zenodo v1.0.0 (2019-06-18, 8,630 non-redundant), DOI [10.5281/zenodo.3235802](https://doi.org/10.5281/zenodo.3235802). Superseded by 20260317, which holds all but 62 of its accessions; kept for comparison. `zenodo_checksums.txt`: MD5s from the Zenodo record, all matching. |
| `genbank_backbone/records.fasta.gz` | The full GenBank record behind every backbone accession (both versions; 15,170), fetched by `bin/fetch_backbone_records.py` from NCBI E-utilities on 2026-11-09. `bin/its2.py` re-cuts the ITS2 from these. |

Primers: NC1 and NC2, Gasser et al. 1993 (Nucleic Acids Res 21:2525),
as printed with their Illumina adapter versions on
https://www.nemabiome.ca/sequencing.

Also seen, not used:

- https://ucvm-bio.shinyapps.io/Nematode_ITS1-5-8S-ITS2 - nemabiome's
  ITS1-5.8S-ITS2 database, v0.0.1 (2023-04-14, 10,107 non-redundant), for
  NC5/NC2 long amplicons. Our nanopore runs use NC1/NC2.
- https://figshare.unimelb.edu.au/articles/dataset/Nemabiome_ITS_Database/28013753 -
  a general nematode rRNA set (18S, ITS, 28S) built for NanoCLUST; not the
  curated ITS2 database despite the name.
