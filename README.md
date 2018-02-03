# myBiopipe
NGS data analysis pipeline write with snakemake.
The pipe accepts regular fq/fq.gz/fastq/fqstq.gz as input, and finally outputs an annotated snp_indel variation for each sample. Of course, Refgenome should be changed according to you analysis.
### This pipe could be devided into 3 steps:
#### 1. Preparation
All the scripts should be put and excuted in top of the directory, the catalog tree can be seen in __example__, finally the work directory may look like this:
```
.
├── config.yaml
├── samples.txt
├── seperate_sample.py
├── wgs.py
│── work.sh
├── dag.pdf
├── log
│   ├── KPGP-00001_annotate.log
│   ├── KPGP-00001_bwa_map.log
│   ├── KPGP-00001_call_var.log
│   ├── KPGP-00001.log
│   ├── KPGP-00001samtools_index.log
│   ├── KPGP-00001_samtools_merge_bam.log
│   ├── KPGP-00001_samtools_rmdup.log
│   └── KPGP-00001_samtools_sort_bam.log
├── mapping
│   └── KPGP-00001
│       ├── KPGP-00001_final.bam
│       └── KPGP-00001_final.bam.bai
├── rawdata
│   └── KPGP-00001
│       ├── KPGP-00001_L2_R1.fq
│       ├── KPGP-00001_L2_R2.fq
│       ├── KPGP-00001_L3_R1.fq
│       ├── KPGP-00001_L3_R2.fq
│       ├── KPGP-00001_L4_R1.fq
│       └── KPGP-00001_L4_R2.fq
└── variation
    ├── KPGP-00001.anno.vcf.exonic_variant_function
    ├── KPGP-00001.anno.vcf.log
    ├── KPGP-00001.anno.vcf.variant_function
    └── KPGP-00001_raw_var.vcf
```   
#### 2. Seperate all samples
This step will mkdir for each sample, and move all fq files to its directory. This is accomplished by `seperate_sample.py`:

```shell
cd example
python3 seperate_sample.py KPGP-00001
```
if there are over 1 sample, seperated by commas.
#### 3. Data analysis
Modify the config file `config.yaml` according to your analysis.Then start analysis:

`sh work.sh`