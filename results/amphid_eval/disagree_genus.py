"""disagree.py <acc dir> <test prefix>: amphid's species calls vs the other classifiers' on the
same sequences (DADA2 ASV matched exactly: merged insert = ASV; 2x150 windows = ASV halves less
the 3-base spacer allowance). Read-weighted over amphid reads."""
import csv, glob, sys, collections, os
d, P = sys.argv[1], sys.argv[2]
RANKS = ['superkingdom','kingdom','phylum','class','order','family','genus','species']
def dada(f):
    out = {}
    for r in csv.DictReader(open(f), delimiter='\t'):
        cols = [c for c in r if not c.startswith('boot') and c != 'asv'][:8]
        rank, name = 'unassigned', ''
        for i, c in enumerate(cols):
            v = r[c]
            if v and v != 'NA': rank, name = RANKS[i], v.replace('_', ' ')
            else: break
        out[r['asv']] = (rank, name)
    return out
asv = {r['asv']: r['sequence'] for r in csv.DictReader(open(P + '.asvs.tsv'), delimiter='\t')}
key = {}
for a, s in asv.items():
    if 'N' * 10 in s:
        x, _, y = s.partition('N' * 10)
        key[(x[:127], y[-127:])] = a
    else:
        key[s] = a
arms = {}
for f, lab in ((P + '.tax.nemits2_dev.trim28s.tsv', 'dada2 x this db, 28S cut'),
               (P + '.idtaxa.nemits2_dev.trim28s.tsv', 'idtaxa x this db, 28S cut'),
               (P + '.tax.nemabiome_20260317.tsv', 'dada2 x nemabiome'),
               (P + '.idtaxa.nemabiome_20260317.tsv', 'idtaxa x nemabiome')):
    if os.path.exists(f): arms[lab] = dada(f)
ours = {r['asv']: (r['gap_rank'], r['gap_name']) for r in csv.DictReader(open(P + '.ours.nemits2_dev.tsv'), delimiter='\t')}
arms['reference gap rule (TOP=40)'] = ours
tot = matched = 0; both = collections.Counter(); dis = collections.Counter(); ex = collections.defaultdict(collections.Counter)
for p in glob.glob(d + '/*.clusters.tsv'):
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r['origin'] == 'chimera' or r['rank'] not in ('species','genus'): continue
        n = int(r['reads']); tot += n
        x, _, y = r['consensus'].partition('N' * 10)
        a = key.get((x, y) if y else x)
        if a is None: continue
        matched += n
        for lab, arm in arms.items():
            rk, nm = arm.get(a, ('', ''))
            if rk in ('species','genus'):
                both[lab] += n
                g1 = r['lineage'].split(';')[6] if len(r['lineage'].split(';'))>6 else ''; g2 = nm.split(' ')[0]
                if g1 != g2: dis[lab] += n; ex[lab][(r['taxon'], nm)] += n
print(f"{os.path.basename(d)}: amphid species reads {tot:,}; matched to a DADA2 ASV {100*matched/tot:.1f}%")
for lab in arms:
    print(f"  vs {lab:28s} both species {both[lab]:>9,}  disagree {dis[lab]:>7,} ({100*dis[lab]/max(1,both[lab]):.2f}%)",
          "; ".join(f"{a} vs {b} {n:,}" for (a, b), n in ex[lab].most_common(3)))
