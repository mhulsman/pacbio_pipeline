################################################################
# Script to quickly check whether the analysis of the samples  #
# went smooth and produced all expected files of something was #
# missed along the way.                                        #
# This script was written after we had problems on Spider for  #
# the submission of jobs through slurm that was not working    #
# properly anymore.                                            #
#                                                              #
# The script prints on screen whether each sample has all the  #
# expected files or not. If not, it also list the files that   #
# are expected but not found.                                  #
#                                                              #
# To run the script, make sure to have loaded the conda py37   #
# environment, then type the following:                        #
#                                                              #
# Rscript check_files.R                                        #
################################################################

# Libraries
library(data.table)
library(stringr)

# Functions

# main
# 1. retrieve all files in the directory of interest
path = "/project/holstegelab/Share/nicco/ad_chc_project/"
path_nijmegen = "/project/holstegelab/Share/nicco/ad_chc_project/nijmegen_samples/"
flist = system(paste0("ls ", path, "*.gz"), intern = T)
flist_nijmegen = system(paste0("ls ", path_nijmegen, "*.gz"), intern = T)
all_files = c(flist, flist_nijmegen)

# 2. loop on files -- there should be 11 files for each of these
for (f in all_files){
	# split name -- we just need the initial run name
	sample_name = strsplit(f, "\\.")[[1]][1]
        cat(paste0("# Checking sample --> ", sample_name, "\n"))
	# list all files (and store them) with this name
	flist_sample = system(paste0("ls ", sample_name, "*"), inter = T)
	# check the length
	if (length(flist_sample) != 11){
		# place here the expected output extensions
		expected = c(".aligned_hifi_reads.bam", ".aligned_hifi_reads.bam.bai", ".aligned_sub_reads.bam", ".aligned_sub_reads.bam.bai", ".hifi_reads.bam", ".sub_reads.bam", ".subreads.ccs.bam", ".subreads.ccs.bam.pbi", ".subreads.ccs.log", ".subreads.ccs.report.txt", ".subreads.ccs.zmw_metrics.json.gz")
		# if something is weird, let's see what it is
		expected_correct = paste0(sample_name, expected)
		# check who's missing there
		missing = expected_correct[which(!(expected_correct %in% flist_sample))]
		cat(paste0("## A file is missing for this sample --> ", missing, "\n\n"))
	} else {
	cat(paste0("## No files missing for this sample.\n\n"))
    }
}