# sources/

Downloaded source data. Nothing here is committed except this file; each
build's `manifest.json` records which backbone it used.

| path | from |
|---|---|
| `nemabiome_its2_20260317/` | **The backbone.** The live nemabiome ITS2 database, version 20260317 (17 March 2026; 27,343 sequences, 15,108 non-redundant), built with markerDB. Only downloadable through its web app, https://ucvm-bio.shinyapps.io/Nematode_ITS2 (embedded at https://www.nemabiome.ca/its2-database): show every column including `sequence`, tick `excel` (and `dada2` as a cross-check), Download. The app returns `database.zip`; `raw/` holds it unzipped, and `bin/nemabiome_export.py raw/database.xlsx Nematode_ITS2_20260317` writes the `_seqs_nr.fasta` / `_taxonomy_nr.txt` pair. `checksums.txt`: SHA-256 of the zip and its contents as downloaded 2026-11-09. CC BY 4.0. |
| `nemabiome_its2_1.0.0/` | Zenodo v1.0.0 (2019-06-18, 8,630 non-redundant), DOI [10.5281/zenodo.3235802](https://doi.org/10.5281/zenodo.3235802). Superseded by 20260317, which holds all but 62 of its accessions; kept for comparison. `zenodo_checksums.txt`: MD5s from the Zenodo record, all matching. |
| `genbank_update_20261113/` | **The GenBank update** used by v0.2.0: nematode records the backbone lacks, fetched by `bin/genbank_update.py` from NCBI E-utilities on 2026-11-13 (ITS2 and 5.8S/28S records published since the backbone's 17 March 2026 snapshot, records annotated with 5.8S and 28S but not ITS2, every "Genus sp. X" ITS2 record, and the accessions in `curation/genbank_named.tsv`); 5,897 holding an ITS2 region, 3,894 of them informal. Its ITS2 regions were estimated against the development build of that day. The records carry the source label `genbank_update_20260317`, which the code takes from the "since" date. |
| `genbank_update_20260317/` | The update used by v0.1.4 (fetched 2026-11-11, the 26 Metastrongylus records added 2026-11-13; 1,860 records with an ITS2 region, no informal names). Superseded by `genbank_update_20261113`, which holds all of it but two records with no ITS2 (AJ512836, PZ600301), with identical GenBank sequences and labels. |
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
