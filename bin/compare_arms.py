#!/usr/bin/env python3
"""Compare classifier x database arms on the same DADA2 ASVs.

    bin/compare_arms.py <prefix>

<prefix> is one run of bin/dada2_nemabiome.R (<prefix>.asvs.tsv,
.counts.tsv), with the calls made on its ASVs:
  <prefix>.tax.<ref>.tsv               assignTaxonomy (bin/assign_taxonomy.R)
  <prefix>.idtaxa.<ref>.tsv            IdTaxa (bin/idtaxa.R)
  <prefix>.ours.<ref>.tsv              bin/classify_asvs.py (span, gap rules)
  <prefix>.{tax,idtaxa}.nemits2_dev.trim28s.tsv   on ASVs cut by bin/trim_28s.py
<ref> is nemabiome_20260317 (the backbone's dada2 export) or nemits2_dev
(this build, bin/export_dada2.py). Prints the share of reads named to
species, genus and above per arm, where two arms both name species and
disagree, and each arm's commonest names.
"""
import csv,sys,collections
P=sys.argv[1]
reads=collections.Counter()
for r in csv.DictReader(open(P+'.counts.tsv'),delimiter='\t'): reads[r['asv']]+=int(r['reads'])
tot=sum(reads.values())
R=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
RANKS=['superkingdom','kingdom','phylum','class','order','family','genus','species']
def dada(f):
    return collections.defaultdict(lambda: ('unassigned',''), _dada(f))
def _dada(f):
    """assignTaxonomy names only seven columns (the eighth comes out as 'NA'):
    read the rank columns by position against the reference's eight ranks."""
    out={}
    for r in csv.DictReader(open(f),delimiter='\t'):
        cols=[c for c in r if not c.startswith('boot') and c!='asv'][:8]
        rank,name='unassigned',''
        for i,c in enumerate(cols):
            v=r[c]
            if v and v!='NA': rank,name=RANKS[i],v.replace('_',' ')
            else: break
        out[r['asv']]=(rank,name)
    return out
def ours(f,rule='exact'):
    return {r['asv']:(r[rule+'_rank'],r[rule+'_name']) for r in csv.DictReader(open(f),delimiter='\t')}
import os
arms={'dada2 x nemabiome':dada(P+'.tax.nemabiome_20260317.tsv'),'dada2 x this db':dada(P+'.tax.nemits2_dev.tsv'),
      'span x nemabiome':ours(P+'.ours.nemabiome_20260317.tsv'),'span x this db':ours(P+'.ours.nemits2_dev.tsv'),
      'gap x nemabiome':ours(P+'.ours.nemabiome_20260317.tsv','gap'),'gap x this db':ours(P+'.ours.nemits2_dev.tsv','gap')}
for f,lab in ((P+'.tax.nemits2_dev.trim28s.tsv','dada2 x this db, 28S cut'),(P+'.idtaxa.nemits2_dev.trim28s.tsv','idtaxa x this db, 28S cut')):
    if os.path.exists(f): arms[lab]=dada(f)
for ref in ('nemabiome_20260317','nemits2_dev'):
    f=P+f'.idtaxa.{ref}.tsv'
    if os.path.exists(f): arms['idtaxa x '+('nemabiome' if ref.startswith('nemabiome') else 'this db')]=dada(f)
print(f'{P}: {len(reads)} ASVs, {tot:,} reads')
print(f"{'arm':26s} species  genus  higher/unassigned")
for k,a in arms.items():
    sp=sum(n for x,n in reads.items() if a[x][0] in('species','species ')); 
    ge=sum(n for x,n in reads.items() if a[x][0]=='genus')
    extra=''
    print(f'{k:26s} {100*sp/tot:6.1f}% {100*ge/tot:6.1f}% {100*(tot-sp-ge)/tot:6.1f}%'+extra)
def comp(a,b,level):
    agree=dis=0; ex=collections.Counter()
    for x,n in reads.items():
        (ra,na),(rb,nb)=arms[a][x],arms[b][x]
        if ra=='species' and rb=='species':
            if na==nb: agree+=n
            else: dis+=n; ex[(na,nb)]+=n
    return agree,dis,ex
for a,b in [x for x in [('idtaxa x this db, 28S cut','gap x this db'),('idtaxa x nemabiome','gap x this db'),('dada2 x nemabiome','gap x this db'),('idtaxa x this db','gap x this db'),('dada2 x this db','gap x this db')] if x[0] in arms]:
    ag,di,ex=comp(a,b,'species')
    print(f'\n{a} vs {b}: both at species {100*(ag+di)/tot:.1f}% of reads; disagree on {100*di/tot:.1f}%')
    for (x,y),n in ex.most_common(6): print(f'   {n:8,d}  {x}  vs  {y}')
# species-level composition per arm (top)
print()
for k,a in arms.items():
    c=collections.Counter()
    for x,n in reads.items(): c[a[x][1] if a[x][0]=='species' else f'({a[x][0]}) {a[x][1]}']+=n
    print(k,'|',', '.join(f'{nm} {100*v/tot:.1f}%' for nm,v in c.most_common(7)))
