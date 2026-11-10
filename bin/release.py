#!/usr/bin/env python3
"""Cut a release from a versioned build.

    bin/release.py <build dir> <releases dir> <version>

<build dir> is bin/build.py's output, built with --version <version> from
committed code (a "-dirty" code commit is refused). Writes
<releases dir>/nemits2_<version>/ - the build's release files, plus the
k-mer classifier files for the DADA2 workflow (bin/export_kmer.sh:
nemits2_<version>_dada2.fasta, nemits2_<version>_idtaxa.RData,
trim_28s.R) - and a .tar.gz of it next to it. The manifest gains the
k-mer files' checksums and the R package versions they were made with;
the summary gains their names. A release is never overwritten. Only the
manifest and summary are committed (.gitignore); the tarball goes to the
GitHub release.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def r_versions():
    lib = os.path.join(HERE, os.pardir, "work", "R-lib")
    code = (f'.libPaths(c("{lib}", .libPaths())); '
            'cat(R.version$major, ".", R.version$minor, " ", '
            'as.character(packageVersion("dada2")), " ", '
            'as.character(packageVersion("DECIPHER")), sep = "")')
    out = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=True).stdout
    r, dada2, decipher = out.split()
    return {"R": r, "dada2": dada2, "DECIPHER": decipher}


def main() -> int:
    build, releases, version = sys.argv[1:4]
    man = json.load(open(os.path.join(build, "manifest.json")))
    if man["version"] != version:
        sys.exit(f"{build} was built as {man['version']}, not {version}: rebuild with --version")
    if man["code_commit"].endswith("-dirty") or man["code_commit"] == "uncommitted":
        sys.exit(f"{build} was built from uncommitted code ({man['code_commit']}): commit, rebuild")
    name = f"nemits2_{version}"
    out = os.path.join(releases, name)
    if os.path.exists(out):
        sys.exit(f"{out} exists: a release is never overwritten")
    os.makedirs(out)
    for f in list(man["files"]) + ["manifest.json"]:
        shutil.copy2(os.path.join(build, f), out)
    for f, h in man["files"].items():
        if sha256(os.path.join(out, f)) != h:
            sys.exit(f"{f} does not match the build's manifest")

    subprocess.run([os.path.join(HERE, "export_kmer.sh"), out, out, name], check=True)
    kmer = [f"{name}_dada2.fasta", f"{name}_idtaxa.RData", "trim_28s.R"]
    with open(os.path.join(out, "summary.txt"), "a") as fh:
        fh.write("\nFor the DADA2 workflow (README: For DADA2 and IdTaxa users):\n"
                 + "".join(f"  {f}\n" for f in kmer)
                 + "  Run trim_28s.R on the ASVs first: records stop before the 28S stretch.\n")
    man["kmer"] = {"files": kmer, "made_with": r_versions(),
                   "idtaxa": "LearnTaxa on the dada2 reference; IdTaxa(strand = 'both', threshold = 60) "
                             "as nemabiome.ca runs it",
                   "asvs": "cut with trim_28s.R before classifying"}
    for f in kmer + ["summary.txt"]:
        man["files"][f] = sha256(os.path.join(out, f))
    man["name"] = name
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(man, fh, indent=1)

    tar = out + ".tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        t.add(out, arcname=name)
    print(open(os.path.join(out, "summary.txt")).read())
    print(f"{tar}  {os.path.getsize(tar) / 1e6:.1f} MB  sha256 {sha256(tar)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
