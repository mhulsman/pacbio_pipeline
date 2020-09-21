import os
import csv

wildcard_constraints:
    sample="[\w\d_\-]+"


SMARTTOOLS="/home/hulsmanm/env/sys_enhance/opt/smartlink9/smrtcmds/bin"
BAMSIEVE=os.path.join(SMARTTOOLS, "bamsieve")
CCS=os.path.join(SMARTTOOLS, "ccs")
PBMM2=os.path.join(SMARTTOOLS, "pbmm2")

H38CCS='/projects/0/qtholstg/marc/pacbio_test/h38_ccs.mmi'
H38SUBREAD='/projects/0/qtholstg/marc/pacbio_test/h38_subread.mmi'

PB_SOURCE="/projects/0/qtholstg/pacbio/"
CACHE=None


PYTHONPATH="/nfs/home2/hulsmanm/env/sys_enhance/lib/python27.zip:/nfs/home2/hulsmanm/env/sys_enhance/bin/lib/python2.7:/nfs/home2/hulsmanm/env/sys_enhance/lib/python2.7/plat-linux2:/nfs/home2/hulsmanm/env/sys_enhance/lib/python2.7/lib-tk:/nfs/home2/hulsmanm/env/sys_enhance/lib/python2.7/lib-old:/nfs/home2/hulsmanm/env/sys_enhance/lib/python2.7/lib-dynload:/nfs/home2/hulsmanm/env/sys_enhance/lib/python2.7/site-packages"
MYPATH= os.path.expanduser('~/projects/pacbio')

def generate_done(wildcards):
    samples = get_sample_listing('.')
    path = os.getcwd()
    return [os.path.join(path, e['sample'], 'processing.done') for e in samples.values()]


rule all:
    input:
        generate_done
    output:
        touch('all.done')



def get_sample_listing(dirname):
    global CACHE
    if CACHE is None:
        import csv
        f = open(os.path.join(dirname,'sample_listing.tsv'),'r')
        r = csv.reader(f,delimiter='\t')
        CACHE = {}
        for line in r:
            CACHE[line[0]] = {'run': line[1], 'sample': line[0], 'sub_run': line[2], 'filename':line[3]}
    return CACHE


rule generate_samplelisting:
    output:
        "{directory}/sample_listing.tsv"
    shell: """
        find {PB_SOURCE} | grep -P "subreads.bam$" | awk -F'/' '{{split($NF,a,"."); print a[1] "\t" $6 "\t" $7 "\t" $0}}' | sort -V > {output}
        """







rule sample_done:
    input:
        "{directory}/{sample}/ccs.1p.stats.txt",
        "{directory}/{sample}/bam_extract.stats.tsv"
    output:
        "{directory}/{sample}/processing.done"
    shell: """
        rclone mkdir RD:100plus/collaborations/pacbio/Pacbio_QC/stats/{wildcards.sample}
        rclone copy {input[0]} RD:100plus/collaborations/pacbio/Pacbio_QC/stats/{wildcards.sample}/
        rclone copy {input[1]} RD:100plus/collaborations/pacbio/Pacbio_QC/stats/{wildcards.sample}/
        touch {output}
           """
    

def get_subreads(wildcards):
    record = get_sample_listing(wildcards['directory'])[wildcards['sample']]
    return record['filename']


rule bamsieve:
    input:
        get_subreads
    output:
        "{directory}/{sample}/subreads.1p.bam"
    shell: """
           {BAMSIEVE} --percentage 1 {input} {output}
           """


rule ccs_1p:
    input:
        "{directory}/{sample}/subreads.1p.bam"
    output:
        "{directory}/{sample}/ccs.1p.bam",
        "{directory}/{sample}/ccs.1p.stats.txt"
    threads: 24
    shell: """
         {CCS} -j 24 {input} {output[0]} --report-file {output[1]}
        """

rule align_ccs_1p:
    input:
        "{directory}/{sample}/ccs.1p.bam"
    output:
        "{directory}/{sample}/ccs.1p.aligned.bam"
    threads: 24
    shell: """
           {PBMM2} align --preset CCS {H38CCS} --unmapped {input} {output}
        """
 
rule align_subreads_1p:
    input:
        "{directory}/{sample}/subreads.1p.bam"
    output:
        "{directory}/{sample}/subreads.1p.aligned.bam"
    threads: 24
    shell: """
           {PBMM2} align --preset SUBREAD {H38SUBREAD} --unmapped {input} {output}
        """


rule bam_extract_ccs_1p:
    input:
        "{directory}/{sample}/ccs.1p.aligned.bam"
    output:
        "{directory}/{sample}/ccs.1p.aligned.bam_extract.tsv"
    threads: 12
    shell: """
            samtools view -@ 24 {input} | awk '{{printf $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5 "\\t" $6 "\\t" length($10) "\\t" $13 "\\t" $14 "\\t" $15 "\\t" $16 "\\t" $17 "\\t" $18 "\\t" $19 "\\n"}}' > {output}
            """

rule bam_extract_subread_1p:
    input:
        "{directory}/{sample}/subreads.1p.aligned.bam"
    output:
        "{directory}/{sample}/subreads.1p.aligned.bam_extract.tsv"
    threads: 12
    shell: """
            samtools view -@ 24 {input}| awk '{{printf $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5 "\\t" $6 "\\t" length($10) "\\t" $12 "\\t" $14 "\\t" $16 "\\t" $17 "\\t" $18 "\\t" $19 "\\t" $20 "\\t" $21 "\\t" $22 "\\n"}}' > {output}
            """
rule analyze_stats:
    input:
        "{directory}/{sample}/subreads.1p.aligned.bam_extract.tsv",
        "{directory}/{sample}/ccs.1p.aligned.bam_extract.tsv"
    output:
        "{directory}/{sample}/bam_extract.stats.tsv"
    shell: """
            PYTHONPATH={PYTHONPATH} python {MYPATH}/analyze.py {input[0]} {input[1]} {output}
          """
        
