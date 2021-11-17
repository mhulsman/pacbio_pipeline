#!/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/python3

#################################################
# Script to extract subreads composing the hifi #
# reads to investigate mosaicism.               #
#                                               #
# This operation is quite computationally heavy #
# as it needs to go through the whole output of #
# the Pacbio (~500-1000Gb) and extract only a   #
# subset of reads.                              #
# For this reason, this script will attempt to  #
# extract only 1 region (HTT gene) which showed #
# some mosaicism. However, if multiple regions  #
# need to be extracted later on, this is quite  #
# easy to adapt.                                #
#################################################

# Libraries
import os
import sys
import pysam
from multiprocessing import Pool
import pandas as pd
import re
import pickle

# Functions
## Function to save obj to file
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

## Function to load obj from file
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

## Function that given a sample, a path and a set of zmw, look for the reads of interest
def findSubreads(sample):
    ## extract sample id, zmw and path
    sample_id, zmw, path = sample.split(";")
    zmw = sorted([int(x) for x in zmw.split(",")])
    ## define output
    sample_reads = {}
    outname = sample_id + ".subreads_HTT"
    ## start to read bam file
    bamfile = pysam.AlignmentFile(path, "rb", check_sq=False)
    ## start the loop on all reads
    for read in bamfile:
        zmw_id = int(read.qname.split("/")[1])
        if zmw_id in zmw:
            seq = read.query
            if zmw_id in sample_reads.keys():
                print("## Found match for ZMW %s" %(zmw_id))
                sample_reads[zmw_id].append(seq)
            else:
                sample_reads[zmw_id] = [seq]
    ## save object
    save_obj(sample_reads, outname)
    return "Done"

# Main
## subreads data is huge --> better to focus on a specific region across samples --> a good candidate is HTT region
reference = pd.read_csv('/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/known_repeats_diseases.csv', sep = ";")
htt_region = reference[reference['Gene'] == "HTT"]

## take the processed tandem repeats to derive the ZMW id of each read
tandem_repeats = pd.read_csv('/project/holstegelab/Share/nicco/workspaces/20211013_target_approach/tandem_repeats_outputs_onlyMatchesMotif.txt', sep = "\t")
tandem_repeats_htt = tandem_repeats[tandem_repeats["REGION"] == "chr4;3074876;3074940;CAG (CTG);HTT;HD;6-29;29-37;38-180"]

## create a dictionary of each smart cell and the relative ZMW ids for the HTT region of interest
all_read_ids = list(set(tandem_repeats_htt["READ_ID"]))
smart_cells = {}
for x in all_read_ids:
    x = x.split("/")
    smrt_cell, zmw = x[0][1:], x[1]
    if smrt_cell in smart_cells.keys():
        smart_cells[smrt_cell].append(zmw)
    else:
        smart_cells[smrt_cell] = [zmw]

## find all samples for which we have subreads data
all_files_subreads_1 = [x.rstrip() for x in os.popen("find /project/holstegelab/Data/pacbio/ -name '*.subreads.bam'")]
all_files_subreads_2 = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/nicco/additional_data_20201218/ -name '*.subreads.bam'")]
all_files_subreads = all_files_subreads_1 + all_files_subreads_2

## match only the smart cells we are interested in
files_interest = []
for x in all_files_subreads:
    tmp = x.split("/")[-1].split(".subreads")[0]
    if tmp in smart_cells.keys():
        tmp_res = tmp + ";" + ",".join(smart_cells[tmp]) + ";" + x
        smart_cells[tmp].append(x)
        files_interest.append(tmp_res)

## nice now we have the smart cells to look into, the path to the file, and the zmw to look for
with Pool(10) as p:
    print(p.map(findSubreads, files_interest))
## wait until everything is finished

## To submit the job through slurm, type the following:
## sbatch --mail-type=NONE --ntasks=1 --cpus-per-task=10 /project/holstegelab/Software/nicco/bin/Snakemake_pipeline/extract_subreads.py