import os.path
from icecream import ic
from glob import glob

DEFAULT_ROOT = 'tree_decomp'

def subset_alignment_paths(root):
    i = 0
    while os.path.isfile(os.path.join(root, "backbone", "crucible_test","subsets", f"{i}.afa")):
        aln = os.path.join(root, "backbone", "crucible_test","subsets", f"{i}.afa")
        
        i += 1
        assert os.path.isfile(aln)
        yield aln

def query_path(root):
    return os.path.join(root, "backbone", "queries.fasta")

class SpillingContainer():
    def __init__(self, capacity):
        self.capacity = capacity
        self.data = []
    
    def see(self, item):
        heapq.heappush(self.data, item)
        if len(self.data) > self.capacity:
            heapq.heappop(self.data)