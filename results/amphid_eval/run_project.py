"""run_project.py <project> <read_length> [extra assign_db_pairs args...]:
amphid's assign_db_pairs on every sample of a test project; 2x150 cut from the raw 2x251."""
import gzip, os, subprocess, sys, glob
from concurrent.futures import ThreadPoolExecutor
S = os.path.dirname(os.path.abspath(__file__))
prj, L = sys.argv[1], int(sys.argv[2])
extra = sys.argv[3:]
tag = os.environ.get("TAG", "")
T = f"/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/testdata/{prj}"
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
A = "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline"
out = f"{S}/acc/{prj}.{L}{tag}"
os.makedirs(out, exist_ok=True)
src = T
if L != 250:
    src = f"{S}/reads{L}/{prj}"
    os.makedirs(src, exist_ok=True)
    def cut(f):
        d = os.path.join(src, os.path.basename(f))
        if os.path.exists(d):
            return
        with gzip.open(f, "rt") as a, gzip.open(d + ".tmp", "wt", compresslevel=1) as b:
            for i, line in enumerate(a):
                b.write(line.rstrip("\n")[:L] + "\n" if i % 2 == 1 else line)
        os.rename(d + ".tmp", d)
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(cut, glob.glob(f"{T}/*.fastq.gz")))
def run(r1):
    s = os.path.basename(r1)[:-len("_1.fastq.gz")]
    with open(f"{out}/{s}.log", "w") as log:
        subprocess.run([sys.executable, f"{A}/bin/assign_db_pairs.py", "--r1", r1,
                        "--r2", r1.replace("_1.fastq.gz", "_2.fastq.gz"), "--database", R,
                        "--primers", f"{S}/nc_primers.tsv", "--sample-id", s,
                        "--read-length", str(L), "--primer-offset", "3", "--outdir", out,
                        *extra], stdout=log, stderr=log, check=True)
with ThreadPoolExecutor(8) as ex:
    list(ex.map(run, sorted(glob.glob(f"{src}/*_1.fastq.gz"))))
print(prj, L, len(glob.glob(f"{out}/*.clusters.tsv")), "samples")
