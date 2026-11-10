import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import weak_census as w
from build import span_distance
sets = json.load(open("narrow_sets.json")); sus = set(sets["suspect"])
seven = ["ITS2_299134f6aa91","ITS2_2cd32996fa9c","ITS2_4d886bfc22b7","ITS2_759f5f1467ea","ITS2_819cc284c6ca","ITS2_d2974ff6ffa8","ITS2_d4fe3c4eeb3f"]
hp = [i for i, t in w.tax.items() if t["species"] == "Haemonchus placei"]
for i in seven:
    rid, d = w.nearest(i)
    rob = [j for j in w.by["Haemonchus contortus"] if j not in w.weak]
    j = min(rob, key=lambda j: span_distance(w.seqs[i], w.seqs[j]))
    dp = min(span_distance(w.seqs[j], w.seqs[k]) for k in hp)
    print(i, "class:", "suspect" if i in sus else "redundant" if d <= 1 else "middle", f"| {d} edit(s) from robust {j} ({w.tax[j]['grade']}, {w.nm[j]} acc), which is {dp} from H. placei")
