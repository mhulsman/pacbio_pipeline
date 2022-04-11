#!/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/python3

# calculate read length of all hifi reads and sub reads and also add sample information
########################################################################################

# libraries
import pysam
import subprocess
import os
import statistics as stats
import sys
import re
from matplotlib import rcParams
from tsmoothie.smoother import *
from tsmoothie.utils_func import sim_randomwalk
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
from random import random
from multiprocessing import Pool

# functions
# function to make bins for each chromosome -- 1 bin every 10kb
def makeBins(binsize):
    # read chromosome lengths
    chr_lengths = open("/project/holstegelab/Share/nicco/reference_files/hg38_chromosome_lengths.txt").readlines()
    chr_bins = {}
    for line in chr_lengths:
        line = line.rstrip().split()
        chrom, length = "chr" + str(line[0]), int(line[1])
        for i in range(0, length, binsize):
            if i == 0:
                chr_bins[chrom] = [[]]
            else:
                chr_bins[chrom].append([])
    return chr_bins

# function to assign reads to bins depending on chromosomal position and overlap
def findBins(all_bins, chrom, start, end, binsize, rl):
    # identify bin for start and end position separately
    b_start = int(start/binsize)
    b_end = int(end/binsize)
    # if the bin of the start and end are the same, we need to update the corresponsind bin with the read lenght
    if b_start == b_end:
        all_bins[chrom][b_start].append(rl)
    # otherwise we need to look at the overlap and assign to one of the bins based on larger overlap
    else:
        # calculate overlaps
        b_start_overl = (b_start*binsize + binsize) - start
        b_end_overl = end - (b_end*binsize)
        if b_start_overl >= b_end_overl:
            all_bins[chrom][b_start].append(rl)
        else:
            all_bins[chrom][b_end].append(rl)
    return all_bins

# function to look into the ccs algorithm output specifically
def lookIntoCCS():
    # 5. find all reports
    ccs_reports = os.popen("find /project/holstegelab/Share/nicco/ad_chc_project/ -name '*report*'").read().split("\n")[:-1]

    # 6 produce usable output
    outf = open('/project/holstegelab/Share/nicco/ad_chc_project/stats/summary_ccs_algorithm.txt', "w")
    header = "SAMPLE\ZMW_INPUT\tZMW_PASS_PC\tZMW_FAIL_PC\tZMW_SHORTCUTS_PC\tZMW_WITH_TR_PC\tBELOW_SNR_PC\tLACK_FULL_PASS_PC\tHETERO_INS_PC\tCOV_DROP_PC\tINSUF_COV_PC\tDRAFT_ERROR_PC\tDRAFT_ABOVE_PC\tDRAFT_BELOW_PC\tPOLISH_FAIL_PC\tEMPTY_COV_WIND_PC\tCCS_NO_CONV_PC\tCCS_BELOW_RQ_PC\tUNKNOWN_ERROR_PC\n"
    outf.write(header)

    # 7. main loop across files
    ccs_stats = {}
    for sample in ccs_reports:
        print("## Working on --> %s" %(sample))
        sample_name = sample.split("/")[-1].replace(".txt", "")
        with open(sample) as fname:
            for line in fname:
                if line.startswith("ZMWs input"):
                    line = line.rstrip().split(":")
                    qc_id, qc_num = line[0], line[-1]
                    ccs_stats[sample_name] = [(qc_id, qc_num, "NA")]
                elif line.startswith("Exclusive counts"):
                    pass
                else:
                    line = line.rstrip().split(":")
                    if line == [""]:
                        pass
                    else:
                        qc_id = line[0]
                        pc = float(line[-1].split(" ")[-1].replace("(", "").replace(")", "").replace("%", ""))
                        nm = int(line[-1].split(" ")[1].replace("(", "").replace(")", "").replace("%", ""))
                        ccs_stats[sample_name].append((qc_id, nm, pc))
        # finally write things to file
        outf.write(sample_name)
        outf.write("\t")
        outf.write(ccs_stats[sample_name][0][1])
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][1][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][2][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][3][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][4][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][5][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][7][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][8][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][9][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][10][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][12][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][13][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][14][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][15][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][16][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][17][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][18][-1]))
        outf.write("\t")
        outf.write(str(ccs_stats[sample_name][19][-1]))
        outf.write("\n")
    outf.close()
    return("Done")

# Function to save obj to file
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

# Function to load obj from file
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

# function to calculate bam statistics using multiple processors to go faster
def BAMstats(sample):
    statistics = {}
    print("## Sample of interest is %s" %(sample))
    # 3. loop on each file and store the read length and other info
    read_lengths = []
    number_passes = []
    read_quality = []
    bamfile = pysam.AlignmentFile(sample, "rb")
    for read in bamfile:
        # extract read length and reference positions
        rl = read.query_length
        chrom = str(read.reference_name)
        # also extract number of passes and read quality
        info = read.tags
        np = info[2]
        rq = info[3]
        # assign to lists
        number_passes.append(np[-1])
        read_lengths.append(rl)
        read_quality.append(rq[-1])
    print("## All data is stored. Computing some statistics now..")
    # calculate overall statistics -- mean, median and 95% confidence interval
    if len(read_lengths) >1:
        statistics[sample] = [ stats.mean(read_lengths), stats.median(read_lengths), sms.DescrStatsW(read_lengths).tconfint_mean(), len(read_lengths), stats.mean(number_passes), stats.median(number_passes), sms.DescrStatsW(number_passes).tconfint_mean(), stats.mean(read_quality), stats.median(read_quality), sms.DescrStatsW(read_quality).tconfint_mean() ]
    else:
        statistics[sample] = [ "NA", "NA", "NA", len(read_lengths), "NA", "NA", "NA", "NA", "NA", "NA" ]
    save_obj(statistics, '/project/holstegelab/Share/nicco/ad_chc_project/stats/stats_' + str(random()).replace('.', ''))
    print("## Done with %s. Moving to next sample\n" %(sample))
    bamfile.close()
    return("Done!")

# function to take care of the sample check
def parseSampleCheck(sample_check_filelist, path_sample_checks):
    samples_check_info = {}
    problems = []
    for f in sample_check_filelist:
        sample_check = pd.read_csv(path_sample_checks + f, sep = "\t")          # read file into pandas df
        sample_check.sort_values("ID_GWAS", inplace = True)              # sorting by first name
        sample_check.drop_duplicates(subset ="ID_GWAS", keep = 'first', inplace = True)    # drop all dups
        sample_check.sort_values("PERC_HOMOLOGY", inplace = True, ascending=False, ignore_index=True)              # sorting by percentage again
        if (sample_check['PERC_HOMOLOGY'][0] < 0.90):
            print("!!! Something is weird for sample --> %s: percentage of similarity is low !!!" %(f))
        else:
            if (sample_check['PERC_HOMOLOGY'][0] - sample_check['PERC_HOMOLOGY'][1]) > 0.10:            # check if sample is correctly recognized
                sb = sample_check.iloc[0, ]
                tmp_info = sb.values.tolist()
                fname = f.split(".")[0]
                samples_check_info[fname] = tmp_info
            else:
                print("!!! Something is weird for sample --> %s: multiple IDs matching !!! \n%s" %(f, sample_check.head()))
                problems.append([f, sample_check])
    return(samples_check_info, problems)


# main
# 1. find all aligned bam files of hifi reads
aligned_bam_hifi = os.popen("find /project/holstegelab/Share/nicco/ad_chc_project/ -name '*aligned_*bam'").read().split("\n")[:-1]

# 2. main loop across samples
n_cores = 20
with Pool(n_cores) as p:
    print(p.map(BAMstats, aligned_bam_hifi))

# 3. read all files in and put things together
pkl_filelist = os.popen("find /project/holstegelab/Share/nicco/ad_chc_project/stats/ -name '*pkl'").read().split("\n")[:-1]
all_stats = {}
for f in pkl_filelist:
    tmp = load_obj(f[:-4])
    k = list(tmp.keys())[0]
    all_stats[k] = tmp[k]
    os.system("rm %s" %(f))
save_obj(all_stats, '/project/holstegelab/Share/nicco/ad_chc_project/stats/all_stats_aligned_bam')

# 4. before writing, good to find sample information from the sample check
path_sample_checks = "/project/holstegelab/Software/nicco/bin/Snakemake_pipeline/sample_check_snakemake/checks/"
sample_check_filelist = os.popen("ls /project/holstegelab/Software/nicco/bin/Snakemake_pipeline/sample_check_snakemake/checks/").read().split("\n")[:-1]
samples_check_info, problems = parseSampleCheck(sample_check_filelist, path_sample_checks)
save_obj(samples_check_info, '/project/holstegelab/Share/nicco/ad_chc_project/stats/all_stats_sample_check')

# 4. write to file
outf = open('/project/holstegelab/Share/nicco/ad_chc_project/stats/summary_aligned_bam_files.txt', "w")
header = "SAMPLE\tSAMPLE_TYPE\tSAMPLE_SOURCE\tN_READS\tMEAN_READ_LENGTH\tMEDIAN_READ_LENGTH\t95%_CI_READ_LENGTH\tMEAN_PASSES\tMEDIAN_PASSES\t95%_CI_PASSES\tMEAN_QUALITY\tMEDIAN_QUALITY\t95%_CI_QUALITY\tSAMPLE_MATCH\tperc_homology\tn_snps\tindividual_id\tID_100plus\tI_ID\tID_progressPD\tID_NBB\tID_EMIF_90plus_twins\tStudy\tchip\tdiagnosis\tsex\tage\tdiagnosis_ADC_first\tage_ADC_firstvisit\tPD_phenotype_description\n".upper()
outf.write(header)
for sample in all_stats.keys():
    if '/' in sample:
        sample_splt = sample.split("/")
        sample_name, sample_type = sample_splt[-1].split(".")[0], sample_splt[-1].split(".")[1].split("_")[1]
        sample_source = "radbound" if "e_" in sample_name else "vumc"
        mean_r, median_r, confint_r, n, mean_p, median_p, confint_p, mean_q, median_q, confint_q = all_stats[sample][0], all_stats[sample][1], all_stats[sample][2], all_stats[sample][3], all_stats[sample][4], all_stats[sample][5], all_stats[sample][6], all_stats[sample][7], all_stats[sample][8], all_stats[sample][9]
        if sample_name in samples_check_info.keys():
            sample_match, pc_h, n_snps, individual_id, ID_100plus, I_ID, ID_progressPD, ID_NBB, ID_EMIF_90plus_twins, Study, chip, diagnosis, sex, age, diagnosis_ADC_first, age_ADC_firstvisit, PD_phenotype_description = samples_check_info[sample_name]
        else:
            sample_match, pc_h, n_snps, individual_id, ID_100plus, I_ID, ID_progressPD, ID_NBB, ID_EMIF_90plus_twins, Study, chip, diagnosis, sex, age, diagnosis_ADC_first, age_ADC_firstvisit, PD_phenotype_description = ["NA" for i in range(17)]
        tmp_string = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(str(sample_name), str(sample_type), str(sample_source), str(n), str(mean_r), str(median_r), str(confint_r), str(mean_p), str(median_p), str(confint_p), str(mean_q), str(median_q), str(confint_q), str(sample_match), str(pc_h), str(n_snps), str(individual_id), str(ID_100plus), str(I_ID), str(ID_progressPD), str(ID_NBB), str(ID_EMIF_90plus_twins), str(Study), str(chip), str(diagnosis), str(sex), str(age), str(diagnosis_ADC_first), str(age_ADC_firstvisit), str(PD_phenotype_description))
        outf.write(tmp_string)
outf.close()

# to execute on the cluster, run the following command
#sbatch --mail-type=NONE --ntasks=1 --cpus-per-task=10 /project/holstegelab/Share/nicco/ad_chc_project/stats/stats_aligned_and_ccs_MP_v2.py