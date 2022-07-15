import os.path
from icecream import ic
from glob import glob

DEFAULT_ROOT = 'tree_decomp'

def subset_alignment_paths(root):
    i = 0
    while os.path.isdir(os.path.join(root, "root", f"A_0_{i}")):
        alns = glob(os.path.join(root, "root", f"A_0_{i}","*.fasta"))
        i += 1
        assert os.path.isfile(alns[0])
        yield alns[0]

def query_path(root):
    return os.path.join(root, "backbone", "queries.fasta")