from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from tqdm import tqdm
from swAlign import swAlign
from Bio import pairwise2
import Levenshtein

config = {"adapter" : "CTACACGACGCTCTTCCGATCT", #"AGATCGGAAGAGCGTCGTGTAGGCCGTAATGGCCTTTAGT",
          "adapter_RV" : "AGATCGGAAGAGCGTCGTGTAG" } #"ACTAAAGGCCATTACGGCCTACACGACGCTCTTCCGATCT"}

match = 1
mismatch = -1
gap_open = -1
gap_extend = -1

def aligner(fasta,
           batch_size = 500,
           thread = 30):
    
    
    read_batchs = fasta_batch_generator(fasta, batch_size)
    
    rst = multiprocessing_submit(align_batch,
                    read_batchs, n_process=thread, 
                    pbar_update=batch_size)
    
    fasta_ed = []
    
    for idx, f in enumerate(rst):
        fasta_ed.extend(f.result())
    return fasta_ed


def align_batch(read_batch):
    e_distances = []
    for read in read_batch:
        seq1_align, seq2_align, alignment_score, ed = swAlign(config["adapter"], read)
        seq1_align_rv, seq2_align_rv, alignment_score_rv, ed_rv = swAlign(config["adapter_RV"], read)

        #aligned_read, _, _, start, end = pairwise2.align.localms(read, config["adapter"], match, mismatch, gap_open, gap_extend)[0]
        #ed = Levenshtein.distance(aligned_read[start:end], config["adapter"])
        
        #aligned_read, _, _, start, end = pairwise2.align.localms(read, config["adapter_RV"], match, mismatch, gap_open, gap_extend)[0]
        #ed_rv = Levenshtein.distance(aligned_read[start:end], config["adapter_RV"])
        if ed < ed_rv:
                ED = ed + max(len(adapter) - len(seq1_align_rv.rstrip()), 0)
        else:
            ED = ed_rv + max(len(adapter) - len(seq1_align.rstrip()), 0)
        e_distances.append(ED)
            
    return e_distances


def fasta_batch_generator(fasta, batch_size):
    open_fn = gzip.open if fasta.endswith('.gz') else open

    with open_fn(fasta, 'rt') as fasta:
        batch = []
        line_num = 0
        for line in fasta:
            line_num += 1
            if line_num % 2 == 0:
                batch.append(line.rstrip())
                if len(batch) == batch_size:
                    yield batch
                    batch = []

        if len(batch) > 0:
            yield batch


def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,pbar = True, pbar_update = 500,  *arg, **kwargs):
    executor = ProcessPoolExecutor(n_process)

    max_queue = n_process * 2
    if pbar:
        pbar = tqdm(unit = 'read', desc='Processed')

    futures = {}
    n_job_in_queue = 0
    while True:
        while n_job_in_queue < max_queue:
            i = next(iterator, None)
            if not i:
                break
            futures[executor.submit(func, i, *arg, **kwargs)] = None
            n_job_in_queue += 1

        job = next(as_completed(futures), None)

        # no more job  
        if job is None:
            break
        else:
            n_job_in_queue -= 1
            pbar.update(pbar_update)
            yield job
            del futures[job]


def batch_iterator(iterator, batch_size):
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch