#!/usr/bin/env python3
"""The README figures, as static SVG in light and dark variants.

    bin/readme_figures.py <build or release dir> docs/figures [README.md]

As vespa-18s-refdb's bin/readme_figures.py (whose chart code this
copies), for ITS2: one source, so no source chart. Writes

1. grade_<mode>.svg  records per nematode group, split by grade: how much
                     of the ITS2 each record holds (complete, near
                     complete, partial, unmeasured)
2. taxa_<mode>.svg   records per group, split into those that resolve to
                     species read whole and those that stop higher

and builds the README section that shows them - the <picture> blocks,
alt text from the same numbers, and tables of species and resolution
per group - replacing the section between <!-- figures:start --> and
<!-- figures:end --> in README.md if given, else printing it.

Groups are families gathered the way a parasitologist would: animal
parasites by the usual clades (bold in the figures), then
insect-parasitic, plant-parasitic and fungal-feeding, and free-living
nematodes. Standard library only.
"""

from __future__ import annotations

import collections
import csv
import json
import os
import sys
from xml.sax.saxutils import escape

FONT = 'system-ui, -apple-system, &quot;Segoe UI&quot;, sans-serif'
MODES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", s1="#2a78d6", s2="#eb6834",
                  seq_dark="#256abf", seq_light="#86b6ef", grey="#b4b3ac",
                  ring="rgba(11,11,11,0.10)"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", s1="#3987e5", s2="#d95926",
                 seq_dark="#3987e5", seq_light="#184f95", grey="#55544f",
                 ring="rgba(255,255,255,0.10)"),
}

# Nematode clades (Blaxter et al. 1998, Nature 392:71; De Ley & Blaxter
# 2002 for the suborders). NCBI files clades III-V all under Rhabditida,
# so those are assigned by family; the rest by order.
CLADES = [
    ("Clade I: Dorylaimia", "whipworms, Trichinella, capillarids; dorylaims, longidorids"),
    ("Clade II: Enoplia", "enoplids, trichodorids"),
    ("Clade III: Spirurina", "ascarids, anisakids, filarial and spirurid worms, pinworms"),
    ("Clade IV: Tylenchina", "Strongyloides, Rhabdias, Steinernema, plant parasites, cephalobids"),
    ("Clade V: Rhabditina", "strongyles, hookworms, lungworms, Heterorhabditis, rhabditids"),
    ("Other Chromadoria", "monhysterids and other free-living chromadorids"),
]
ORDER_CLADE = {
    "Trichinellida": 0, "Dorylaimida": 0, "Mononchida": 0, "Dioctophymatida": 0,
    "Mermithida": 0, "Marimermithida": 0,
    "Enoplida": 1, "Triplonchida": 1,
    "Monhysterida": 5, "Plectida": 5, "Araeolaimida": 5, "Chromadorida": 5, "Desmodorida": 5,
}
FAMILY_CLADE = {f: k for k, fams in (
    (2, "Anisakidae Raphidascarididae Toxocaridae Ascarididae Ascaridiidae Heterakidae "
        "Aspidoderidae Acanthocheilidae Heterocheilidae Cucullanidae Kathlaniidae Cosmocercidae "
        "Subuluridae Onchocercidae Setariidae Filariidae Gnathostomatidae Habronematidae "
        "Gongylonematidae Philometridae Thelaziidae Camallanidae Physalopteridae Tetrameridae "
        "Spirocercidae Acuariidae Cystidicolidae Diplotriaenidae Dracunculidae Oxyuridae "
        "Heteroxynematidae Pharyngodonidae Rhigonematidae Carnoyidae Ichthyocephalidae "
        "Xustrostomatidae Travassosinematidae"),
    (3, "Strongyloididae Rhabdiasidae Alloionematidae Steinernematidae Panagrolaimidae "
        "Brevibuccidae Cephalobidae Heteroderidae Pratylenchidae Hoplolaimidae Aphelenchoididae "
        "Meloidogynidae Anguinidae Tylenchulidae Criconematidae Telotylenchidae Tylenchidae "
        "Hemicycliophoridae Belonolaimidae Merliniidae Sphaeronematidae Aphelenchidae "
        "Nothotylenchidae Ecphyadophoridae Neotylenchidae Sphaerulariidae Iotonchiidae"),
    (4, "Strongylidae Trichostrongylidae Ancylostomatidae Metastrongylidae Molineidae "
        "Heligmosomidae Crenosomatidae Pseudaliidae Filaroididae Syngamidae Amidostomidae "
        "Heligmonellidae Stephanuridae Diaphanocephalidae Chabertiidae Angiostrongylidae "
        "Amphibiophilidae Heterorhabditidae Rhabditidae Diplogasteridae Neodiplogasteridae "
        "Agfidae Angiostomatidae"),
) for f in fams.split()}
GRADES = [("Complete", ("complete",), "seq_dark"),
          ("Near complete (≥ 98%)", ("near_complete",), "seq_light"),
          ("Partial", ("partial_90", "partial_75", "partial"), "s2"),
          ("Unmeasured", ("unmeasured",), "grey")]


# Lifestyle by family; anything not listed is free-living. Insect parasites
# (Steinernema, Heterorhabditis, sphaerulariids) count as animal
# parasites, fungal feeders (Aphelenchoides, Bursaphelenchus) as plant.
PLANT = set("""Heteroderidae Pratylenchidae Hoplolaimidae Aphelenchoididae Meloidogynidae
    Anguinidae Tylenchulidae Criconematidae Telotylenchidae Tylenchidae Hemicycliophoridae
    Belonolaimidae Merliniidae Sphaeronematidae Aphelenchidae Nothotylenchidae Ecphyadophoridae
    Longidoridae Trichodoridae""".split())
ANIMAL_FAMILIES = set(" ".join(fams for k, fams in (
    (0, "Trichuridae Trichinellidae Capillariidae Dioctophymatidae Marimermithidae Mermithidae"),
    (1, "Steinernematidae Heterorhabditidae Sphaerulariidae Iotonchiidae Neotylenchidae "
        "Strongyloididae Rhabdiasidae Angiostomatidae Amphibiophilidae"),
)).split()) | {f for f, k in FAMILY_CLADE.items() if k in (2, 4)} - {
    "Rhabditidae", "Diplogasteridae", "Neodiplogasteridae", "Agfidae"}
LIFESTYLES = ["animal parasites", "plant parasites", "free-living"]
ROMAN = ["I", "II", "III", "IV", "V", "Other"]


def clade_of(t: dict):
    k = ORDER_CLADE.get(t["order"], FAMILY_CLADE.get(t["family"]))
    return k


def lifestyle_of(t: dict) -> str:
    f = t["family"]
    return ("animal parasites" if f in ANIMAL_FAMILIES else
            "plant parasites" if f in PLANT else "free-living")


def group_of(t: dict) -> str:
    k = clade_of(t)
    if k is None:
        return "Family not placed"
    pre = f"{ROMAN[k]} · " if k < 5 else "Other Chromadoria · "
    return pre + lifestyle_of(t)


def nice_ticks(vmax: float):
    for step in (1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500,
                 5000, 10000, 20000, 25000, 50000):
        if vmax / step <= 5:
            n = int(vmax // step) + 1
            return [step * i for i in range(n + 1)], step * n
    return [0, vmax], vmax


def fmt(n: int) -> str:
    return f"{n:,}"


def tick_label(v: float) -> str:
    return f"{v / 1000:g}k" if v >= 1000 else f"{v:g}"


def hbar_svg(title, subtitle, rows, seg_names, seg_colors, mode, note=None,
             emphasise=(), value_label=None, right=92):
    """rows: [(label, [seg values])]. Stacked horizontal bars, one axis.
    value_label(vals) -> the text at the bar end (default: the total)."""
    c = MODES[mode]
    W, left = 760, 196
    top, row_h, bar_h = 104, 28, 16
    plot_w = W - left - right
    H = top + row_h * len(rows) + 46 + (22 if note else 0)
    ticks, vmax = nice_ticks(max(sum(v) for _, v in rows))
    x = lambda v: left + plot_w * v / vmax  # noqa: E731
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}" font-family="{FONT}" role="img" '
         f'aria-label="{escape(title)}">',
         f'<rect width="{W}" height="{H}" rx="8" fill="{c["surface"]}" stroke="{c["ring"]}"/>',
         f'<text x="24" y="34" font-size="17" font-weight="600" fill="{c["ink"]}">{escape(title)}</text>',
         f'<text x="24" y="56" font-size="13" fill="{c["ink2"]}">{escape(subtitle)}</text>']
    # legend, one row, above the plot
    lx = left
    for name, col in zip(seg_names, seg_colors):
        o.append(f'<rect x="{lx}" y="72" width="12" height="12" rx="3" fill="{c[col]}"/>')
        o.append(f'<text x="{lx + 18}" y="82.5" font-size="12.5" fill="{c["ink2"]}">{escape(name)}</text>')
        lx += 18 + 7.2 * len(name) + 24
    # recessive grid and axis
    y0, y1 = top - 6, top + row_h * len(rows)
    for t in ticks:
        o.append(f'<line x1="{x(t):.1f}" y1="{y0}" x2="{x(t):.1f}" y2="{y1}" '
                 f'stroke="{c["grid" if t else "axis"]}" stroke-width="1"/>')
        o.append(f'<text x="{x(t):.1f}" y="{y1 + 18}" font-size="11.5" fill="{c["muted"]}" '
                 f'text-anchor="middle" font-variant-numeric="tabular-nums">{tick_label(t)}</text>')
    for k, (label, vals) in enumerate(rows):
        yc = top + row_h * k + row_h / 2
        weight = "600" if label in emphasise else "400"
        o.append(f'<text x="{left - 12}" y="{yc + 4.5}" font-size="13" font-weight="{weight}" '
                 f'fill="{c["ink"]}" text-anchor="end">{escape(label)}</text>')
        segs = [(v, col) for v, col in zip(vals, seg_colors) if v > 0]
        xa = x(0)
        for j, (v, col) in enumerate(segs):
            x0 = xa
            # Every non-zero segment stays visible - a source's handful of
            # records only it provides must not vanish beside PR2's 89,000.
            # The label carries the exact values.
            w = max(x(v) - x(0), 4)
            xa = x0 + w
            last = j == len(segs) - 1
            gap = 0 if last else 2  # surface gap between stacked segments
            ww = max(w - gap, 1.5)
            if last:  # rounded data end, square at the baseline / joint
                r = min(4, ww / 2, bar_h / 2)
                d = (f"M{x0:.1f},{yc - bar_h / 2:.1f} h{ww - r:.1f} "
                     f"a{r},{r} 0 0 1 {r},{r} v{bar_h - 2 * r:.1f} "
                     f"a{r},{r} 0 0 1 -{r},{r} h-{ww - r:.1f} z")
                o.append(f'<path d="{d}" fill="{c[col]}"/>')
            else:
                o.append(f'<rect x="{x0:.1f}" y="{yc - bar_h / 2:.1f}" width="{ww:.1f}" '
                         f'height="{bar_h}" fill="{c[col]}"/>')
        text = value_label(vals) if value_label else fmt(sum(vals))
        o.append(f'<text x="{xa + 8:.1f}" y="{yc + 4.5}" font-size="12.5" '
                 f'fill="{c["ink2"]}" font-variant-numeric="tabular-nums">{escape(text)}</text>')
    if note:
        o.append(f'<text x="24" y="{H - 16}" font-size="12" fill="{c["muted"]}">{escape(note)}</text>')
    o.append("</svg>")
    return "\n".join(o)


START, END = "<!-- figures:start -->", "<!-- figures:end -->"
NOT_A_SPECIES = (" x ", " cf.", " aff.", "sensu lato", "environmental sample")


def picture(out, stem, alt):
    return (f'<picture>\n'
            f'  <source media="(prefers-color-scheme: dark)" srcset="{out}/{stem}_dark.svg">\n'
            f'  <img alt="{escape(alt, {chr(34): "&quot;"})}" src="{out}/{stem}_light.svg" '
            f'width="760">\n</picture>')


def main() -> int:
    rel, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = json.load(open(os.path.join(rel, "manifest.json")))["name"]
    with open(os.path.join(rel, "taxonomy.tsv")) as fh:
        tax = list(csv.DictReader(fh, delimiter="\t"))
    order = [(f"{ROMAN[k]} · " if k < 5 else "Other Chromadoria · ") + life
             for k in range(len(CLADES)) for life in LIFESTYLES] + ["Family not placed"]
    by = collections.defaultdict(list)
    for t in tax:
        by[group_of(t)].append(t)
    groups = [g for g in order if by[g]]

    grade_rows = [(g, [sum(t["grade"] in keys for t in by[g]) for _, keys, _ in GRADES])
                  for g in groups]
    # Resolution read whole as the release states it: ties (a record inside
    # another species' record) already stopped at the rank they share.
    with open(os.path.join(rel, "resolution_nanopore.tsv")) as fh:
        res = {r["id"]: r["deepest_rank"] for r in csv.DictReader(fh, delimiter="\t")}
    taxa_rows = [(g, [sum(res[t["id"]] == "species" for t in by[g]),
                      sum(res[t["id"]] != "species" for t in by[g])]) for g in groups]
    key = ("Clades (Blaxter et al. 1998): " + ", ".join(
        f"{ROMAN[k]} {name.split(': ')[1]}" for k, (name, _) in enumerate(CLADES[:5]))
        + ". Bold: animal parasites.")
    for mode in MODES:
        with open(os.path.join(out, f"grade_{mode}.svg"), "w") as fh:
            fh.write(hbar_svg(
                f"How much of the ITS2 the {len(tax):,} records hold",
                f"{name}: records per clade and lifestyle, by grade",
                grade_rows, [n for n, _, _ in GRADES], [c for _, _, c in GRADES], mode,
                note=key,
                emphasise=[g for g in groups if g.endswith("animal parasites")]))
        with open(os.path.join(out, f"taxa_{mode}.svg"), "w") as fh:
            fh.write(hbar_svg(
                "Taxonomic spread of the reference sequences",
                f"{name}: records per clade and lifestyle; resolution read whole, exact read",
                taxa_rows, ["Resolves to species", "Stops at genus or above"],
                ["seq_dark", "seq_light"], mode,
                note=key,
                emphasise=[g for g in groups if g.endswith("animal parasites")]))
    print(f"wrote {out}/grade_{{light,dark}}.svg and taxa_{{light,dark}}.svg", file=sys.stderr)

    def parts(v):
        return ", ".join(f"{n.split(' (')[0].lower()} {x:,}" for (n, _, _), x in zip(GRADES, v))
    grade_alt = (f"Bar chart: records per nematode group in {name}, by how much of the "
                 "ITS2 they hold: "
                 + "; ".join(f"{g} {sum(v):,} ({parts(v)})" for g, v in grade_rows) + ".")
    taxa_alt = (f"Bar chart: records per nematode group in {name}, and how many resolve to "
                "species read whole: " + "; ".join(
                    f"{g} {sum(v):,} ({v[0]:,} to species)" for g, v in taxa_rows) + ".")
    section = "\n\n".join([picture(out, "grade", grade_alt), picture(out, "taxa", taxa_alt),
                           species_table(rel, tax, name, groups),
                           resolution_table(rel, tax, name, groups)])
    if len(sys.argv) > 3:
        readme = sys.argv[3]
        text = open(readme).read()
        a, b = text.find(START), text.find(END)
        if a < 0 or b < a:
            sys.exit(f"{readme}: no {START} ... {END} section to replace")
        open(readme, "w").write(text[:a + len(START)] + "\n" + section + "\n" + text[b:])
        print(f"updated the figures section of {readme}", file=sys.stderr)
    else:
        print(section)
    return 0


def illumina_tables(rel):
    """[(geometry, {id: (rank, name)})] for every resolution_illumina_* table."""
    out = []
    for f in sorted(os.listdir(rel), key=lambda f: (len(f), f)):
        if f.startswith("resolution_illumina_"):
            with open(os.path.join(rel, f)) as fh:
                out.append((f[len("resolution_illumina_"):-4],
                            {r["id"]: (r["deepest_rank"], r["name"])
                             for r in csv.DictReader(fh, delimiter="\t")}))
    return out


def species_table(rel, tax, name, groups):
    """Markdown: species named in the source records, against species the
    release can name read whole and through each Illumina geometry."""
    grp = {t["id"]: group_of(t) for t in tax}
    with open(os.path.join(rel, "resolution_nanopore.tsv")) as fh:
        nano = {r["id"]: (r["deepest_rank"], r["name"]) for r in csv.DictReader(fh, delimiter="\t")}
    ills = illumina_tables(rel)
    rows = ["All records"] + groups
    named, whole, genera = ({g: set() for g in rows} for _ in range(3))
    ill = [{g: set() for g in rows} for _ in ills]
    for t in tax:
        for g in ("All records", grp[t["id"]]):
            if nano[t["id"]][0] == "species":
                whole[g].add(nano[t["id"]][1])
            if t["genus"]:
                genera[g].add(t["genus"])
            for k, (_, res) in enumerate(ills):
                if res[t["id"]][0] == "species":
                    ill[k][g].add(res[t["id"]][1])
    with open(os.path.join(rel, "members.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["role"] != "kept":
                continue
            sp = r["label"].split(";")[-1]
            if sp and not any(x in sp + " " for x in NOT_A_SPECIES):
                for g in ("All records", grp[r["id"]]):
                    named[g].add(sp)
    head = ("| Records | Species named in sources | Species with their own record (read whole) | "
            + " | ".join(f"Species named through Illumina {geo}" for geo, _ in ills)
            + " | Genera |")
    lines = [f"How many species {name} can name. A species named in the source records "
             "but with no record of its own shares its ITS2 with another species (or sits "
             "inside another species' record over its span), so a read from it is named only "
             "to genus or above. Names are counted as written.", "", head,
             "|---|" + "---:|" * (3 + len(ills))]
    for g in rows:
        lab = f"**{g}**" if g.endswith("animal parasites") else g
        lines.append(f"| {lab} | {len(named[g]):,} | {len(whole[g]):,} | "
                     + " | ".join(f"{len(i[g]):,}" for i in ill) + f" | {len(genera[g]):,} |")
    return "\n".join(lines)


def resolution_table(rel, tax, name, groups):
    """Markdown: share of records resolving to species and genus, read
    whole and through each Illumina geometry, per group."""
    with open(os.path.join(rel, "resolution_nanopore.tsv")) as fh:
        nano = {r["id"]: r["deepest_rank"] for r in csv.DictReader(fh, delimiter="\t")}
    ills = [(geo, {i: v[0] for i, v in res.items()}) for geo, res in illumina_tables(rel)]
    reads = [("Nanopore (whole ITS2)", nano)] + [(f"Illumina {geo}", res) for geo, res in ills]
    sets = [("All records", [t["id"] for t in tax])] + \
           [(g, [t["id"] for t in tax if group_of(t) == g]) for g in groups]
    lines = [f"How far a record can be named in {name}, best case (an exact read, no "
             "sequencing error): the share of records resolving to species, and to genus "
             "at best.", "",
             "| Records | " + " | ".join(f"{r}: species | genus" for r, _ in reads) + " |",
             "|---|" + "---:|---:|" * len(reads)]
    for g, ids in sets:
        cells = []
        for _, res in reads:
            n = max(1, len(ids))
            cells += [f"{100 * sum(res[i] == 'species' for i in ids) / n:.0f}%",
                      f"{100 * sum(res[i] == 'genus' for i in ids) / n:.0f}%"]
        lab = f"**{g}**" if g.endswith("animal parasites") else g
        lines.append(f"| {lab} ({len(ids):,}) | " + " | ".join(cells) + " |")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
