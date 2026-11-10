"""ref_top.py <TOP> <prefix>...: rerun nemabiome-its2-refdb's classify_asvs on the DADA2 ASVs with
the candidate cap TOP, and print the gap rule's read-weighted species share."""
import sys, os, csv, collections
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb"
S = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R + "/bin")
def main():
    TOPARG, PREFIXES = sys.argv[1], sys.argv[2:]
    import its2_classify; its2_classify.TOP = int(TOPARG)
    import classify_asvs; classify_asvs.TOP = int(TOPARG)
    for P in PREFIXES:
        out = f"{S}/reftop{TOPARG}.{os.path.basename(P)}.tsv"
        if not os.path.exists(out):
            sys.argv = ["classify_asvs", P + ".asvs.tsv", R + "/work/build_v0.1.0", out, "--workers", "8"]
            classify_asvs.main()
        reads = collections.Counter()
        for r in csv.DictReader(open(P + ".counts.tsv"), delimiter="\t"):
            reads[r["asv"]] += int(r["reads"])
        calls = {r["asv"]: r for r in csv.DictReader(open(out), delimiter="\t")}
        tot = sum(reads.values())
        sp = sum(n for a, n in reads.items() if calls[a]["gap_rank"] == "species")
        print(f"TOP={its2_classify.TOP} {os.path.basename(P)}: gap species {100*sp/tot:.1f}% of {tot:,} reads", flush=True)

if __name__ == '__main__':
    main()
