from pyhmmer import easel, plan7
import pyhmmer
import os.path
import heapq
import logging
from utils import *
from dataclasses import dataclass, field
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple
import numpy as np

# Step 1: parallel convert all subsets into hmms and write them
# Step 2: query all hmms a la hmmsearch (from pyhmmer)
# Step 3: convert scores to adjusted bitscores
# Step 4: create payload

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
            with Path(subset_aln).with_suffix('.h3m').open('wb+') as fh:
                hmm[0].write(fh,True)
    return ehmms, opt_ehmms

def load_query_sequences():
    order = []
    sequences = {}
    with easel.SequenceFile(query_path(DEFAULT_ROOT)) as sf:
        for s in sf:
            ds = s.digitize(ALPH)
            sequences[ds.name] = ds
            order.append(ds.name)
    return order, sequences

def load_crucible_metadata():
    with open(os.path.join(DEFAULT_ROOT, "backbone", "crucible_test","melt.json"), "r") as fh:
        return json.load(fh)

@dataclass
class AdjustedBitscoreTracker:
    hmm_ids : List[int] = field(default_factory=list)
    hmm_sizes : List[int] = field(default_factory=list)
    bitscores : List[float] = field(default_factory=list)

    # 98% from Chengze's calculation
    def calc_adjusted_scores(self, limit = 7) -> List[Tuple[int, float]]:
        bitscores = np.array(self.bitscores)
        sizes = np.array(self.hmm_sizes)
        indices = self.hmm_ids
        unsorted_hmm_ids = []
        unsorted_scores = []
        for i in range(len(bitscores)):
            score_i, size_i = bitscores[i], sizes[i]
            exponents = bitscores - score_i \
                    + np.log2(sizes / size_i)
            denominator = np.sum(np.power(2, exponents))
            hmm_id = indices[i]
            unsorted_hmm_ids.append(hmm_id)
            unsorted_scores.append(-1. / denominator) # very ugly to negate things...
        if len(bitscores) > limit:
            # then we need to repartition
            ix_arr = np.argpartition(unsorted_scores, limit - 1)[:limit]
            return [(unsorted_hmm_ids[ix], -unsorted_scores[ix]) for ix in ix_arr]
        else:
            return [(a, -neg_b) for a, neg_b in zip(unsorted_hmm_ids, unsorted_scores)]

def serialize_payload(payload):
    with open(os.path.join(DEFAULT_ROOT, "backbone", "crucible_test","scores.json"), "w+") as fh:
        json.dump(payload, fh)

def main():
    _, opt_ehmms = build_ehmms()
    order, query_sequences = load_query_sequences()
    metadata = load_crucible_metadata()
    start = time.time()
    query_results = defaultdict(AdjustedBitscoreTracker)
    for tophits in pyhmmer.hmmsearch(opt_ehmms.values(), query_sequences.values(), E=99999999,bias_filter=False, F1=1.0, F2=1.0, F3=1.0,alphabet = ALPH):
        hmmid = int(tophits.query_name)
        (seq_lb, seq_ub) = metadata['metadata'][hmmid]['sequence_range']
        hmmsize = seq_ub - seq_lb
        for hit in tophits:
            seq_name = hit.name
            res = query_results[seq_name]
            res.hmm_ids.append(hmmid)
            res.hmm_sizes.append(hmmsize)
            res.bitscores.append(hit.score)
    end = time.time(); logging.info(f"HMMSearch Phase took {end - start} seconds...")
    scores = {k:v.calc_adjusted_scores() for k, v in query_results.items()}
    payload = []
    for taxon in order:
        payload.append(scores[taxon])
    serialize_payload(payload)
    
        # hmm_name = tophits.query_name
    #     # ic(len(tophits.to_msa(ALPH).alignment))
    #     for hit in tophits:
    #         query_top_matches[hit.name].see((-hit.score, hmm_name))
    # per_hmm_queries = {k:[] for k in ehmms}
    # for sequence_name, top_matches_outer in query_top_matches.items():
    #     top_matches = top_matches_outer.data
    #     for neg_score, hmm_name in top_matches:
    #         per_hmm_queries[hmm_name].append((-neg_score, sequence_name))
    # hit_results = defaultdict(list) # mapping from sequence name to a list of aligned version of that sequences with weights
    # aligner = plan7.TraceAligner()
    # for hmm_name, ground_queries in per_hmm_queries.items():
    #     hmm = ehmms[hmm_name]
    #     seqs = [query_sequences[qname] for w, qname in ground_queries]
    #     traces = aligner.compute_traces(hmm, seqs)
    #     text_msa = aligner.align_traces(hmm, seqs, traces, all_consensus_cols = True)
    #     for aligned_row in zip(text_msa.sequences):
    #         print(aligned_row, "...")

if __name__ == '__main__':
    main()