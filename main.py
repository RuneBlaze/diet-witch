from pyhmmer import easel, plan7
import pyhmmer
import os.path
import heapq
from utils import *
from collections import defaultdict

class SpillingContainer():
    def __init__(self, capacity):
        self.capacity = capacity
        self.data = []
    
    def see(self, item):
        heapq.heappush(self.data, item)
        if len(self.data) > self.capacity:
            heapq.heappop(self.data)

ALPH = easel.Alphabet.dna()

def build_ehmms():
    builder = plan7.Builder(ALPH, ere=0.59, symfrac=0.0)
    background = plan7.Background(ALPH)
    ehmms = {}
    opt_ehmms = {}
    for i, subset_aln in enumerate(subset_alignment_paths(DEFAULT_ROOT)):
        with easel.MSAFile(subset_aln) as msafile:
            msa = msafile.read().digitize(ALPH)
            name = f'{i}'.encode()
            msa.name = name
            hmm = builder.build_msa(msa, background)
            ehmms[name] = hmm[0]
            opt_ehmms[name] = hmm[2]
    return ehmms, opt_ehmms

def load_query_sequences():
    sequences = {}
    with easel.SequenceFile(query_path(DEFAULT_ROOT)) as sf:
        for s in sf:
            ds = s.digitize(ALPH)
            sequences[ds.name] = ds
    return sequences

def main():
    ehmms, opt_ehmms = build_ehmms()
    query_sequences = load_query_sequences()
    query_top_matches = defaultdict(lambda: SpillingContainer(7))
    for tophits in pyhmmer.hmmsearch(opt_ehmms.values(), query_sequences.values(), E=99999999,bias_filter=False, F1=1.0, F2=1.0, F3=1.0,alphabet = ALPH):
        hmm_name = tophits.query_name
        # ic(len(tophits.to_msa(ALPH).alignment))
        for hit in tophits:
            query_top_matches[hit.name].see((-hit.score, hmm_name))
    per_hmm_queries = {k:[] for k in ehmms}
    for sequence_name, top_matches_outer in query_top_matches.items():
        top_matches = top_matches_outer.data
        for neg_score, hmm_name in top_matches:
            per_hmm_queries[hmm_name].append((-neg_score, sequence_name))
    # ic(per_hmm_queries)
    hit_results = defaultdict(list) # mapping from sequence name to a list of aligned version of that sequences with weights
    aligner = plan7.TraceAligner()
    for hmm_name, ground_queries in per_hmm_queries.items():
        hmm = ehmms[hmm_name]
        seqs = [query_sequences[qname] for w, qname in ground_queries]
        traces = aligner.compute_traces(hmm, seqs)
        text_msa = aligner.align_traces(hmm, seqs, traces, all_consensus_cols = True)
        for aligned_row in zip(text_msa.sequences):
            print(aligned_row, "...")

if __name__ == '__main__':
    main()