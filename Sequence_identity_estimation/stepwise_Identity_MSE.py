import numpy as np
from scipy.optimize import minimize
import os
import subprocess
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def call_vsearch(tmp_fasta, adapters_fasta="adapter.fasta",  min_adapter_id=0.7):
    """ """
    tmp_vsearch = tmp_fasta.replace(".fasta", ".vsearch.tsv")

    vsearch_cmd = "vsearch --usearch_global {fasta} --db {adapters} \
    --threads 1 --minseqlength 20 --maxaccepts 5 --id {id} --strand plus \
    --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
    'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
    --userout {output}".format(
        fasta=tmp_fasta,
        id=min_adapter_id,
        adapters=adapters_fasta,
        output=tmp_vsearch,
    )
    stdout, stderr = run_subprocess(vsearch_cmd)
    #print(stderr)
    return tmp_vsearch


def vsearsh_align(fasta, adapter):
    vsearch_cmd = call_vsearch(tmp_fasta=fasta, adapters_fasta=adapter,  min_adapter_id=0.7)
    colnames = [
            "query",
            "target",
            "id",
            "alnlen",
            "mism",
            "opens",
            "qilo",
            "qihi",
            "qstrand",
            "tilo",
            "tihi",
            "ql",
            "tl",
        ]
    #print(vsearch_cmd)
    df = pd.read_csv(vsearch_cmd, sep="\t", header=None, names=colnames)
    df = df[df.id != 0]
    df = df[df.alnlen == 22]
    #os.remove(vsearch_cmd)
    return df.id.mean(), list(df.mism)


def run_subprocess(cmd):
    """
    Run OS command and return stdout & stderr
    """
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)


def call_AsaruSim(AsaruSimPath, perfect_reads, tmp_fastq, params, adapters_fasta="adapter.fasta",  min_adapter_id=0.7):
    """ """
    tmp_fasta = tmp_fastq.replace("fastq", "fasta")
    AsaruSim_cmd = "python3 {AsaruSimPath}/bin/AsaruSim.py call_badread --thread 31 \
    -t  {perfect_reads} \
    -o  {tmp_fastq} \
    --badread-identity {mean},{var},{max} \
    --badread-error-model  {AsaruSimPath}/badread/new_error_model \
    --badread-qscore-model {AsaruSimPath}/badread/new_qscore_model".format(
        perfect_reads=perfect_reads,
        tmp_fastq=tmp_fastq,
        mean=round(params[0],2),
        var=round(params[1],2),
        max=round(params[2],2),
        AsaruSimPath=AsaruSimPath
    )
    stdout, stderr = run_subprocess(AsaruSim_cmd)
    
    seqtk_cmd = "seqtk seq -a {tmp_fq} > {tmp_fa}".format(
        tmp_fq=tmp_fastq,
        tmp_fa=tmp_fasta,
    )
    stdout, stderr = run_subprocess(seqtk_cmd)
    return tmp_fasta


def simulation_function(params): 
    AsarusimPath = "/export/home1/Asarusim/"
    perfect = AsarusimPath+"simulated_data/fasta/sub_template_1k.fasta"
    tmp_fastq = AsarusimPath+"simulated_data/error_profiling/fastq/sub_sim_fited_1k.fastq"
    adapter = AsarusimPath+"simulated_data/error_profiling/fasta/adapter.fasta"
    tmp_fasta = call_AsaruSim(AsaruSimPath, perfect, tmp_fastq, params, adapters_fasta=adapter,  min_adapter_id=0.7)
    _, ed_sim = vsearsh_align(tmp_fasta, adapter)
    #print(tmp_fasta)
    dict_sim = dict(Counter(ed_sim))
    dict_sim = dict(sorted(dict_sim.items()))
    val_sim = list(dict_sim.values())
    return val_sim


def mse_calculator(val_real, val_sim):#params, *args):
    """
    Fonction de coût (MSE) à minimiser.
    """
    if len(val_real) > len(val_sim):
        for i in range(abs(len(val_real)-len(val_sim))):
            val_sim = np.append(val_sim, [0])
            
    if len(val_sim) > len(val_real):
        for i in range(abs(len(val_real)-len(val_sim))):
            val_real = np.append(val_real, [0])
   
    mse_value = np.mean((val_real - val_sim)**2)

    return mse_value


def main():
    # Defining the range
    par = np.arange(94, 100, 0.5)
  
    # Opening a file to write the results in real-time
    file_path = 'MSE_params2.csv'
    
    # Start writing to the file
    with open(file_path, 'w') as file:
        counter = 0
        ed_real = np.array([568, 92, 72, 48, 75, 43, 13])
        for i, mean in enumerate(par):
            for max in par[i+1:]:
                for st in np.arange(1, 6, 1):
                    for rep in range(3):
                        ed_sim = np.array(simulation_function([mean, st, max]))
                        val_mse = mse_calculator(ed_real, ed_sim)
                        line = f"{mean},{st},{max},{rep},{round(val_mse, 2)}\n"
                        file.write(line)
                        counter += 1
                        #print(counter)
                        print(round(val_mse, 2))

if __name__ == '__main__':
    main()
