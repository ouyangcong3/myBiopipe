configfile: "config.yaml"

rawdir = config["rawdata"]
samp = config["sample"]
LANE, = glob_wildcards(rawdir+"/"+samp+"/"+samp+"_L{lane}_R1.fq")

try:
    refgenome = config["refgenome"]
    ref_idx = config["ref_idx"]
except:
    refgenome = "/DataBase/Human/hg19/refgenome/hg19.fa"
    ref_idx = refgenome + ".fai"


##############################################################
############## rules #########################################

rule all:
    input:
        "variation/"+samp+".anno.vcf"

rule bwa_map:
    input:
        ref = refgenome,
        ref_idx = ref_idx,
        R1 = rawdir+"/"+samp+"/"+samp+"_L{lane}_R1.fq",
        R2 = rawdir+"/"+samp+"/"+samp+"_L{lane}_R2.fq"
    output:
        temp("mapping"+samp+"/"+samp+"_L{lane}.bam")
    message:
        "############ start bwa mapping ####################"
    log:
        "log/"+samp+"_bwa_map.log"
    threads: 5
    shell:
        "bwa mem -t {threads} -M {input.ref} {input.R1} {input.R2} | samtools -b -S -t {input.ref_idx} > {output} 2>{log}"

rule samtools_sort_bam:
    input:
        "mapping"+samp+"/"+samp+"_L{lane}.bam"
    output:
        temp("mapping"+"/"+samp+"/"+samp+"_L{lane}.bam.sort")
    message:
        "############## start sorting bam #####################"
    threads: 4
    log:
        "log/"+samp+"_samtools_sort_bam.log"
    shell:
        "samtools sort -@ {threads} -m 4G {input} -o {output} 1>{log} 2>&1"

rule samtools_merge_bam:
    input:
        expand("mapping/"+samp+"/"+samp+"_L{lane}.bam.sort", lane=LANE)
    output:
        temp("mapping"+samp+"/"+samp+".sort.merge.bam")
    message:
        "################ merge bam from lanes ###############"
    log:
        "log/"+samp+"_samtools_merge_bam.log"
    threads: 8
    params:
        header = "mapping"+samp+"/"+samp+"_L"+str(LANE[0])+".bam.sort",
        lane_num = len(LANE)
    shell:
        """
        if [[ ${{params.lane_num}} -gt 1 ]];
        then 
            samtools merge -l 9 -@ {threads} -h {params.header} {output} {input} 2>{log}"
        else
            mv {input} {output}
        fi
        """

rule samtools_rmdup:
    input:
        "mapping"+samp+"/"+samp+".sort.merge.bam"
    output:
        final_bam = protected("mapping"+samp+"/"+samp+"_final.bam"),
        bai = protected("mapping"+samp+"/"+samp+"_final.bam.bai")
    message:
        "################ remove duplicates ################"
    log:
        "log/"+samp+"_samtools_rmdup.log"
    shell:
        "samtools rmdup -S {input} {output.final_bam} && samtools index {output.final_bam} {output.bai} 2>{log}"

rule snp_indel_calling:
    input:
        ref = refgenome,
        bam = "mapping"+samp+"/"+samp+"_final.bam"
    output:
        protected("variation/"+samp+"_raw_var.vcf")
    message:
        "############## call variation #####################"
    log:
        "log/"+samp+"_call_var.log"
    shell:
        """
        samtools mpileup -t DP,AD -q 1 -ugf {input.ref} {input.bam} \
        | bcftools call -vm -O v -o {output} 2>{log}
        """

rule annotate:
    """annotate refgene database"""
    input:
        "variation/"+samp+"_raw_var.vcf"
    output:
        avinput = temp("variation/"+samp+".avinput"),
        anno_var = protected("variation/"+samp+".anno.vcf")
    message:
        "############## snp_indel annotation ###############"
    log:
        "log/"+samp+"_annotate.log"
    shell:
        """
        perl /opt/annovar/convert2annovar.pl -format vcf4 {input} > {output.avinput} && \
        perl /opt/annovar/annotate_variation.pl {output.avinput} -buildver hg19 -geneanno \
        -dbtype refgene /opt/annovar/humandb/ -out {output.anno_var} 2>{log}
        """
