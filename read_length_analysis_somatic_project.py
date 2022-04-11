# caluclate read length of all hifi reads
##########################################

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

# main
# 1. first look at parameters
sample_interest = ["C1_blood", "C1_brain", "C1_child", "C2_blood", "C2_brain", "C2_child", "C3_blood", "C3_brain", "C3_child", "C4_blood", "C4_brain", "C4_child", "C5_blood", "C5_brain", "C5_child", "C6_blood", "C6_brain", "C6_child", "C7_blood", "C7_brain", "C7_child", "C9_blood", "C9_brain", "C9_child", "C10_blood", "C10_brain", "C10_child", "C11_blood", "C11_brain", "C11_child"]

# 2. find all aligned bam files of hifi reads
aligned_bam_hifi = os.popen("find /project/holstegelab/Share/nicco/analysis/ -name '*merged_hifi*bam'").read().split("\n")
aligned_bam_hifi = [x for x in aligned_bam_hifi if len(x)>0]

# 3. main loop across samples
statistics = {}
for sample in sample_interest:
    files_interest = list(filter(lambda x:sample in x, aligned_bam_hifi))
    print("## Sample of interest is %s" %(sample))
    print("## Files of interest are %s" %(len(files_interest)))
    # 4. loop on each file and store the read length
    for smrt in files_interest:
        read_lengths = []
        print("## Started to work on %s" %(smrt))
        bamfile = pysam.AlignmentFile(smrt, "rb")
        for read in bamfile:
            # extract read length and reference positions
            rl = read.query_length
            chrom = str(read.reference_name)
            read_lengths.append(rl)
        print("## All data is stored. Computing some statistics now..")
        # calculate overall statistics -- mean, median and 95% confidence interval
        if len(read_lengths) >1:
            statistics[smrt] = [ stats.mean(read_lengths), stats.median(read_lengths), sms.DescrStatsW(read_lengths).tconfint_mean(), len(read_lengths) ]
        else:
            statistics[smrt] = [ "NA", "NA", "NA", len(read_lengths) ]
        print("## Done with %s. Moving to next smrt cell\n" %(smrt))
        bamfile.close()

# 5. write to file
outf = open('/project/holstegelab/Share/nicco/read_length_analysis/20220331_summary_bam_files_Hifi.txt', "w")
header = "SAMPLE\tSMRT_CELLS\tMEAN_READ_LENGTH\tMEDIAN_READ_LENGTH\t95%_CI_READ_LENGTH\tN_HIFI\n"
outf.write(header)
for sample in statistics.keys():
    sample_splt = sample.split("/")
    sample_name, sample_smrt = sample_splt[-2], sample_splt[-1]
    mean_r, median_r, confint_r, n = statistics[sample][0], statistics[sample][1], statistics[sample][2], statistics[sample][3]
    outf.write(str(sample_name))
    outf.write("\t")
    outf.write(str(sample_smrt))
    outf.write("\t")
    outf.write(str(mean_r))
    outf.write("\t")
    outf.write(str(median_r))
    outf.write("\t")
    outf.write(str(confint_r))
    outf.write("\t")
    outf.write(str(n))
    outf.write("\n")
outf.close()

# 6. run r-script for the plot
#os.system("/usr/bin/Rscript /project/holstegelab/Share/nicco/read_length_analysis/plt_summary.R")

# to execute on the cluster, run the following command
#sbatch --mail-type=NONE --ntasks=1 --cpus-per-task=2 /project/holstegelab/Share/nicco/read_length_analysis/read_length_analysis_v2.sh