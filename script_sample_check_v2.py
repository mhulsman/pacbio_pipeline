#####################################################
# Script to perform sample-check in python          #
#                                                   #
# The script will extract 20,000 SNPs from Pacbio   #
# data (no correction for technical issues here)    #
# and will compare the genotypes of all these SNPs  #
# with those of all individuals for which we have   #
# genotype data for. This means, at the moment,     #
# N~10,000 samples.                                 #
#                                                   #
# The script will:                                  #
#   1. extract SNPs from Pacbio                     #
#   2. clean them based on coverage (coverage > 3)  #
#       and remove all multiallelic SNPs            #
#   3. check alleles of Pacbio with all samples     #
#       with genotype data                          #
#                                                   #
# The script will output a table with the           #
# percentage of homology between Pacbio genotypes   #
# and the array genotypes, for all samples.         #
# Normally, the samples that matches has percentage #
# of homology exceeding 95%, while the second best  #
# is normally around 70%.                           #
#                                                   #
# To run the script, make sure you loaded py37      #
# environment from conda, then type:                #
#                                                   #      
# python3 script_sample_check_v2.py inp.bam out.txt #
#                                                   #
# This script is already part of the general        #
# pipeline for the Pacbio analysis (last step).     #
#####################################################

# Libraries
import os
import re
import pandas as pd
import sys
from multiprocessing import Pool
import collections

# Functions
# function to extract SNPs from pacbio aligned bam
def extractPacbioSNPs(snp):
    try:
        # extract chromosomes and positions
        chrom = snp.split(":")[0]
        startPos = snp.split(":")[1]
        endPos = int(snp.split(":")[1]) + 1
        tmp = os.popen("pysamstats -c chr%s -s %s -e %s -u -t variation --fasta /project/holstegelab/Share/pacbio/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa %s" %(chrom, startPos, endPos, input_file)).read().split("\n")
        header = tmp[0].split("\t")
        info = tmp[1].split("\t")
        tmp_df = pd.concat([pd.Series(x) for x in info], axis=1)
        tmp_df.set_axis(header, axis=1, inplace=True)
        alleles = []
        if tmp_df['A'][0] != '0':
            alleles.append('A')
        if tmp_df['C'][0] != '0':
            alleles.append('C')
        if tmp_df['G'][0] != '0':
            alleles.append('G')
        if tmp_df['T'][0] != '0':
            alleles.append('T')
        alleles = '/'.join(alleles)
        snp_id = str(tmp_df['chrom'][0] + ":" + tmp_df['pos'][0])
        return [snp_id, alleles, int(tmp_df['reads_all'])]
    except:
        return ["NA", "NA", "NA"]

# function to compare pacbio genotypes with array genotypes
def comparePacbioArray(sample):
    if sample in [x for x in range(1000, 15000, 1000)]:
        print("## Processed %s samples thus far.." %(sample))
    count_matches = 0
    count_all = 0
    # get genotypes of that sample
    sample_name, sample_dosage = list(genotypes_matches.loc[sample])[0], list(genotypes_matches.loc[sample])[1:]
    # also get alleles
    alleles = list(genotypes_matches.columns)[1:]
    # check alleles and dosages
    for snp in range(len(alleles)):
        # get pacbio genotype
        pacbio_allele = pacbio_snps_hq[snp][1]
        if len(pacbio_allele) == 1:
            pacbio_genotype = pacbio_allele + "/" + pacbio_allele
        else:
            pacbio_genotype = pacbio_allele
        # get array genotype
        array_effect_allele = alleles[snp].split("_")[-1]
        try:
            array_dosage = float(sample_dosage[snp])
        except:
            array_dosage = "na"
        if array_dosage == 0:
            array_alleles = [alleles[snp].split(":")[-2], alleles[snp].split(":")[-1].split("_")[-2]]
            array_effect_allele = array_alleles[0] if array_alleles[0] != array_effect_allele else array_alleles[1]
            array_dosage = 2 - array_dosage
            array_genotype = array_effect_allele + "/" + array_effect_allele
        elif array_dosage == 1:
            array_alleles = [alleles[snp].split(":")[-2], alleles[snp].split(":")[-1].split("_")[-2]]
            array_genotype = [array_alleles[0] + "/" + array_alleles[1], array_alleles[1] + "/" + array_alleles[0]]
        elif array_dosage == 2:
            array_genotype = array_effect_allele + "/" + array_effect_allele
        else:
            array_genotype = "na"
        # compare genotypes
        if array_genotype == "na":
            pass
        elif pacbio_genotype in array_genotype:
            count_matches += 1
            count_all += 1
        else:
            if isinstance(array_genotype, list):
                array_genotype = "/".join(array_genotype)
            all_alleles = list(set(sorted(list(set(pacbio_genotype.split("/")))) + sorted(array_genotype.split("/"))))
            if len(all_alleles) > 2:
                pass
            else:
                count_all += 1
    if (count_matches > 0) and (count_all > 0):
        pc_match = count_matches / count_all
    else:
        pc_match = "na"
        count_all = "na"
    return([sample_name, pc_match, count_all])

# Main
# 1. manage arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
MAIN = "/project/holstegelab/Software/snakemake_pipeline/sample_check_data/"

# 2. read snps, genotypes and phenotypes
print("## Reading SNPs, Phenotypes and Genotypes...\n")
snps_info = pd.read_csv(MAIN + "random_set_20K_snps.txt", sep = "\t")
phenotypes = pd.read_csv(MAIN + "phenotypes_20211027.txt", sep = "\t")
genotypes = pd.read_csv(MAIN + "dosages_random_set_snps_all_samples.raw.gz", sep = " ")

# 3. extract SNPs from pacbio
snps_to_serch = snps_info['id']
with Pool(4) as p:
    pacbio_snps = p.map(extractPacbioSNPs, snps_to_serch)

# 4. clean pacbio snps -- use also coverage
thr_coverage = 3
pacbio_snps_hq = []
snps_for_comparison = []
for x in pacbio_snps:
    if (x[0] != "NA") and (x[-1] > thr_coverage) and (len(x[1].split("/")) <= 2):
        pacbio_snps_hq.append(x)
        snps_for_comparison.append(x[0])

# 5. now check alleles with gwas
# match pacbio snps with array genotype file
snp_names_array = list(genotypes.columns)
matches = ["IID"]
for snp in pacbio_snps_hq:
    tmp_match = list(filter(lambda x: re.search(r"^%s:" %(snp[0][3:]), x), snp_names_array))
    if len(matches) in [x for x in range(1000, 20000, 1000)]:
        print("## Processed %s snps thus far.." %(len(matches)))
    if len(tmp_match) == 1:
        matches.append(tmp_match[0])
        snp_names_array.remove(matches[-1])
    elif len(tmp_match) > 1:
        print("!!! Possible duplicate for --> %s" %(snp))

# subset genotype file
genotypes_matches = genotypes[matches]
all_samples = list(range(0, genotypes_matches.shape[0]))
with Pool(4) as p:
    res_check = p.map(comparePacbioArray, all_samples)
res_check_df = pd.DataFrame.from_records(res_check)
res_check_df.columns = ['ID_GWAS', 'PERC_HOMOLOGY', 'SNPS_N']

# 6. merge with phenotype information and order by homology
res_check_df['ID_GWAS'] = res_check_df['ID_GWAS'].apply(str)
final_res = pd.merge(res_check_df, phenotypes, on = 'ID_GWAS')
final_res.sort_values(by = ['PERC_HOMOLOGY'], ascending = False, inplace = True)
final_res.to_csv(output_file, sep = "\t", header = True, index = False, na_rep='NA')
