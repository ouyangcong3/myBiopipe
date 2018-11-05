#!/usr/bin/env python3
# coding: utf-8


import os
import sys
from collections import defaultdict
from functools import wraps
import time
import subprocess
import logging
import argparse


def log_info(mesg):
    logging.info(time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + mesg)


def mkdir(dir_):
    if not os.path.exists(dir_):
        os.makedirs(dir_)


def raise_error(cmd, return_code, output):
    if return_code:
        log_info(cmd + ' ERROR: ' + output)
        raise SystemExit(cmd + output)


##############################################################################
################### work directory for intermediate result ###################
##############################################################################
def commands_config():
    # 两种方式组织目录,第一种每一个项目（样本）一个总目录（即总目录名会与项目名一样保持唯一），下面分别放着各个各步的分析结果；
    # 第二种方式是给出一个总目录（此目录名称不需要与某一个项目名称相关），每一个分析步骤会有一个目录，在每一个分析目录下面的结果以项目名称进行区分
    # 这里采用第一种方式
    global logs
    init_dir = outdir + '/' + sample_name
    Mapping = init_dir + '/mapping'
    cleandata = init_dir + '/cleandata'
    variations = init_dir + '/variations'
    logs = init_dir + '/logs'
    for _dir in [init_dir, Mapping, cleandata, variations, logs]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    # 必须的一些输入文件
    refgenome = '/data/sxuan/gatk_bundle/human_g1k_v37/human_g1k_v37.fasta'
    gatk_bundle_b37 = ''
    gatk_bundle_hg19 = ''
    known_sites_1000G_phase1 = '/data/sxuan/gatk_bundle/human_g1k_v37/1000G_phase1.indels.b37.vcf'
    known_sites_Mills_1000G = '/data/sxuan/gatk_bundle/human_g1k_v37/Mills_and_1000G_gold_standard.indels.b37.vcf'
    known_sites_dbsnp = '/data/sxuan/gatk_bundle/human_g1k_v37/dbsnp_138.b37.vcf'

    ###############################################################################
    ####################### Analysis Command ######################################
    ###############################################################################
    # trimming reads with trimmomatic
    clean_R1 = cleandata + '/' + sample_name + '.clean_R1.fq.gz'
    clean_R2 = cleandata + '/' + sample_name + '.clean_R2.fq.gz'
    adapter = {'nextera': '/data/opt/anaconda3/pkgs/trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/NexteraPE-PE.fa',
               'truseq': '/data/opt/anaconda3/pkgs/trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/TruSeq3-PE.fa'}
    trim = '''trimmomatic PE -phred33 -threads {thread} \
                {raw_R1} {raw_R2} \
                {clean_R1} {clean_R1}.unpaired \
                {clean_R2} {clean_R2}.unpaired \
                ILLUMINACLIP:${adapter}:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:5 MINLEN:36
                '''.format(raw_R1=fq1, raw_R2=fq2, clean_R1=clean_R1, clean_R2=clean_R2, adapter=adapter[adapter_type], thread=thread)

    # bwa mapping
    raw_bam = Mapping + '/' + sample_name + '.bam'
    bwa_mapping = 'bwa mem -t {thread} -M -Y -R "@RG\\tID:{RGID}\\tPL:ILLUMINA\\tSM:{samp}" {ref} {clean_R1} {clean_R2} | \
                   samtools view -bS - > {raw_bam}'.format(RGID=sample_name, samp=sample_name, ref=refgenome, clean_R1=clean_R1, clean_R2=clean_R2, thread=thread, raw_bam=raw_bam)
    # samtools sort
    sort_bam = Mapping + '/' + sample_name + '.sort.bam'
    bam_sort = 'samtools sort -@ {thread} -m 4G -O bam -o {sort_bam} {raw_bam}'.format(thread=thread, sort_bam=sort_bam, raw_bam=raw_bam)

    # gatk picard markduplicates
    markdup_bam = Mapping + '/' + sample_name + '.markdup.sort.bam'
    markdup_metrics = Mapping + '/' + sample_name + '.markdup.metrics'
    bam_markdup = 'gatk MarkDuplicates -I {sort_bam} -M {markdup_metrics} -O {markdup_bam}'.format(sort_bam=sort_bam, markdup_metrics=markdup_metrics, markdup_bam=markdup_bam)

    # indel realign (This step is not necessary for GATK4)
    realign_bam = Mapping + '/' + sample_name + '.realign.markdup.sort.bam'
    indel_realign = 'gatk LeftAlignIndels -R {ref} -I {markdup_bam} -O {realign_bam}'.format(ref=refgenome, realign_bam=realign_bam, markdup_bam=markdup_bam)

    # gatk BQSR
    bqsr_bam = Mapping + '/' + sample_name + '.bqsr.bam'
    recal_table = Mapping + '/' + sample_name + '.recal_data.table'
    BaseRecal_and_ApplyBQSR = '''gatk BaseRecalibrator \
                                    -I {markdup_bam} \
                                    -R {ref} \
                                    --known-sites {known_sites1} \
                                    --known-sites {known_sites2} \
                                    --known-sites {known_sites3} \
                                    -O {recal_table} && \
                                gatk ApplyBQSR \
                                     -R {ref} \
                                     -I {markdup_bam} \
                                     --bqsr-recal-file {recal_table} \
                                     -O {bqsr_bam} '''.format(markdup_bam=markdup_bam, ref=refgenome, known_sites1=known_sites_1000G_phase1, bqsr_bam=bqsr_bam,
                                                              known_sites2=known_sites_Mills_1000G, known_sites3=known_sites_dbsnp, recal_table=recal_table)

    # gatk haplotype_caller
    intervals = '-L ' + region if region else ''
    raw_vcf = variations + '/' + sample_name + '.raw.vcf.gz'
    haplotype_caller = '''gatk HaplotypeCaller  \
                            -R {ref} \
                            -I {bqsr_bam} \
                            {intervals} \
                            -O {raw_vcf}'''.format(ref=refgenome, bqsr_bam=bqsr_bam, raw_vcf=raw_vcf, intervals=intervals)

    # gatk Select Variants
    raw_snp = variations + '/' + sample_name + '.raw_snp.vcf.gz'
    raw_indel = variations + '/' + sample_name + '.raw_indel.vcf.gz'
    gatk_selectVariants = '''gatk SelectVariants \
                               -R {ref} \
                               -V {raw_vcf} \
                               --select-type-to-include SNP \
                               -O {raw_snp} && \
                             gatk SelectVariants \
                               -R {ref} \
                               -V {raw_vcf} \
                               --select-type-to-include INDEL \
                               -O {raw_indel}'''.format(ref=refgenome, raw_vcf=raw_vcf, raw_snp=raw_snp, raw_indel=raw_indel)

    # gatk hard filter
    hf_snp = variations + '/' + sample_name + '.hf_snp.vcf.gz'
    hf_indel = variations + '/' + sample_name + '.hf_indel.vcf.gz'
    gatk_hard_filter = '''gatk VariantFiltration \
                            -R {ref} \
                            -V {raw_snp} \
                            -O {hf_snp} \
                            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                            --filter-name "gatk_hardFilters" && \
                          gatk VariantFiltration \
                            -R {ref} \
                            -V {raw_indel} \
                            -O {hf_indel} \
                            --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
                            --filter-name "gatk_hardFilters"'''.format(ref=refgenome, raw_snp=raw_snp, raw_indel=raw_indel, hf_indel=hf_indel, hf_snp=hf_snp)

    ##########################################################################################
    ############################## To Be Finished ############################################
    strelka2_caller = ''
    bcftools_caller = ''
    ##########################################################################################
    ##########################################################################################
    analysis_step_num = {'1': ['trim', trim],
                         '2': ['bwa_mapping', bwa_mapping],
                         '3': ['bam_sort', bam_sort],
                         '4': ['bam_markdup', bam_markdup],
                         '5': ['indel_realign', indel_realign],
                         '6': ['BaseRecal_and_ApplyBQSR', BaseRecal_and_ApplyBQSR],
                         '7.1': ['haplotype_caller', haplotype_caller],
                         '7.2': ['strelka2_caller', strelka2_caller],
                         '7.3': ['bcftools_caller', bcftools_caller],
                         '8': ['gatk_selectVariants', gatk_selectVariants],
                         '9': ['gatk_hard_filter', gatk_hard_filter]
                         }
    return analysis_step_num


def run_cmd(command_name, command):
    global total, logs
    log_info('>>>>>> Start running ' + '*********************** ' + command_name + ' *********************')
    log_info('CMD: ' + ' '.join(command.split()))
    start = time.time()
    return_code, output = subprocess.getstatusoutput(command)
    raise_error(command_name, return_code, output)
    end = time.time()
    with open(logs+'/'+command_name+'.log', 'wt', encoding='utf-8') as log:
        log.write(output)
    log_info(command_name + ' Done')
    elapsed = end - start
    total += elapsed
    hours = elapsed // 3600
    minutes = (elapsed % 3600) // 60
    seconds = (elapsed % 3600) % 60
    step_analysis_time = '{hours}h {minutes}min {seconds:.2f}sec elapsed'.format(hours=hours, minutes=minutes, seconds=seconds)
    log_info(step_analysis_time)


def main():
    global total
    total = 0
    analysis_step_num = commands_config()
    for num in analysis_pipe:
        command_name, command = analysis_step_num[num]
        run_cmd(command_name, command)
    hours = total // 3600
    minutes = (total % 3600) // 60
    seconds = (total % 3600) % 60
    total_analysis_time = '{hours}h {minutes}min {seconds:.2f}sec elapsed'.format(hours=hours, minutes=minutes, seconds=seconds)
    log_info(total_analysis_time)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="germline variations detection pipeline.")
    # parser.add_argument("--project_name", "-pn", type=str, help="Project name for this analysis. Or maybe you could put sample name here.")
    parser.add_argument("--sample_name", "-sn", type=str, required=True, help="sample name will be regard as the prefix of clean fastq, bam and vcf files, RGID is necesary.")
    parser.add_argument("--outdir", "-o", type=str, required=True, help="work directory for results.")
    parser.add_argument("--fastq", "-fq", nargs="+", type=str, help="Input the paried fastq files, eg. '-fq read1.fastq read2.fastq'.")
    parser.add_argument("--bam", "-b", nargs="+", type=str, help="Input the bam files, multi bam for one sample are allowed.")
    parser.add_argument("--region", "-r", default=None, type=str, help='Input region for calling variants.')
    parser.add_argument("--vcf", "-v", type=str, help="Input vcf file for analysis")
    parser.add_argument("--analysis_pipe", '-p', nargs="+", type=str, choices=['1', '2', '3', '4', '5', '6', '7.1', '7.2', '7.3', '8', '9'], default=['1', '2', '3', '4', '6', '7.1', '8', '9'],
                        help="Each number represent an analysis step, {1: 'trim', 2: 'bwa_mapping', 3: 'bam_sort',4: 'bam_markdup',5: 'indel_realign',6: 'BaseRecal_and_ApplyBQSR ',7.1: 'haplotype_caller',7.2: 'strelka2_caller',7.3: 'bcftools_caller',8: 'gatk_selectVariants',9: 'gatk_hard_filter'}")
    parser.add_argument("--adapter", "-a", type=str, choices=['nextera', 'truseq'], default='truseq', help='trimming adapter type.')
    parser.add_argument("--thread", "-t", default=8, type=int, help='Threads to run program.')

    args = vars(parser.parse_args())

    outdir = args["outdir"]
    # project_name = args["project_name"]
    fq1, fq2 = args['fastq']
    bams = args['bam']
    vcf = args['vcf']
    analysis_pipe = args['analysis_pipe']
    sample_name = args['sample_name']
    adapter_type = args['adapter']
    thread = args['thread']
    region = args['region']

    # 日志文件名
    log_file = outdir + '/' + sample_name + '_' + time.strftime('%Y%m%d') + '.log'
    # 配置日志文件和日志级别
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode='w')
    log_info("===================== {samp} Analysi PIPE Start ===========================".format(samp=sample_name))
    init_cmd = 'python3 ' + ' '.join(sys.argv)
    log_info("Initial CMD: " + init_cmd)
    main()
    log_info('Congratulations, ALL DONE!')
