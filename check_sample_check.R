############################################################
# Small script to check if all samples have been correctly #
# sample-checked.                                          #
#                                                          #
# Basically, the script will take all samples for which    #
# aligned hifi data are available and check whether the    #
# sample-check for that specific sample went OK and we     #
# reliable estimates of the real sample name.              #
#                                                          #
# The script prints on screen whether each sample has been #
# successfully checked. If not, it lists the files that    #
# need to be looked up.                                    #
#                                                          #
# To run the script, make sure to load the conda py37      #
# environment, then type the following:                    #
#                                                          #
# Rscript check_sample_check.R                             #
################################################################

# list all files for which we have aligned, hifi bam files
list_bam = system("find /project/holstegelab/Share/nicco/ad_chc_project/ -name '*aligned*hifi*bam'", intern = T)

# list all files for which we have a sample check executed
list_checks = system("find /project/holstegelab/Software/nicco/bin/Snakemake_pipeline/sample_check_snakemake/checks -name '*check*txt'", intern = T)

# do a small quality check of the sample_check data -- see if everything worked fine
quality_check = function(i, list_checks){
    tmp = data.table::fread(list_checks[i], h=T, stringsAsFactors=F, sep="\t")
    if (tmp$PERC_HOMOLOGY[1] > 0.95 & tmp$PERC_HOMOLOGY[2] < 0.90){
        return(cat(paste0("## ", stringr::str_split_fixed(unlist(strsplit(list_checks[i], "/"))[length(unlist(strsplit(list_checks[i], "/")))], "\\.", 2)[, 1], " --> OK\n")))
    } else {
        return(cat(paste0("## ", stringr::str_split_fixed(unlist(strsplit(list_checks[i], "/"))[length(unlist(strsplit(list_checks[i], "/")))], "\\.", 2)[, 1], " --> Have a look!!\n")))
    }
}
for (f in 1:length(list_checks)){ print(quality_check(f, list_checks)) }
all_checked = c()
for (f in list_checks){
    fname = stringr::str_split_fixed(unlist(strsplit(f, "/"))[length(unlist(strsplit(f, "/")))], "\\.sample", 2)[, 1]
    all_checked = c(all_checked, fname)
}

# then check files that failed sample check/have to be run again
correct = c()
non_correct = c()
for (f in list_bam){
    fname = stringr::str_split_fixed(unlist(strsplit(f, "/"))[length(unlist(strsplit(f, "/")))], "\\.aligned", 2)[, 1]
    if (fname %in% all_checked){
        correct = c(correct, fname)
    } else {
        non_correct = c(non_correct, fname)
    }
}

print(paste("You should have a look for these samples --> ", non_correct))