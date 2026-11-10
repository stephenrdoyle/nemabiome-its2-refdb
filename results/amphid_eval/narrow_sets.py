"""Write narrow_sets.json: weak records (partial/unmeasured, 1 accession) in species that also have a
robust record, split by edits to the species' nearest robust record (>5: suspect; <=1: redundant)."""
import json, os, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import weak_census as w
if __name__ == "__main__":
    named = [i for i in w.weak if w.tax[i]["species"]]
    a = [i for i in named if any(j not in w.weak for j in w.by[w.tax[i]["species"]])]
    with Pool(8) as p:
        d = dict(p.map(w.nearest, a, chunksize=20))
    sus = [i for i in a if d[i] > 5]; red = [i for i in a if d[i] <= 1]
    json.dump({"suspect": sus, "suspect+redundant": sus + red}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "narrow_sets.json"), "w"))
    print(len(sus), "suspect,", len(red), "redundant")
