#!/usr/bin/python3
# code: utf-8

"""
if all samples's fastq is in a directory, mkdir for each sample with sample name,
and move all fastq to corresponding directory.

then creat a file named "samples.txt" to store the sample name,
example: samples.txt
a.samp
b.samp
c.samp

Excute this file at the parent directory of rawdata.

Usage: python3 seperate_sample.py samples
eg. python3 seperate_sample.py s1,s2,s3
"""

import os
import shutil
import sys


def seperate(samp_list, all_fq):
    for dir_ in samp_list:
        for fq in all_fq:
            if dir_ in fq:
                shutil.move("rawdata/"+fq, "rawdata/"+dir_)


if sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print("Usage:\tpython3 seperate_sample.py samples\neg. python3 seperate_sample.py s1,s2,s3")
else:
    samp_list = sys.argv[1].split(",")
    all_fq = [f for f in os.listdir('rawdata') if f.endswith((".fastq", "fastq.gz", ".fq", "fq.gz"))]
    # test if there are some wrong input sample names
    in_samp = list(map(lambda samp: 1 if samp in " ".join(all_fq) else 0, samp_list))
    if all(in_samp):
        for samp in samp_list:
            if not os.path.exists("rawdata/"+samp):
                os.mkdir("rawdata/"+samp)
        seperate(samp_list, all_fq)
    else:
        raise SystemExit("sample name ERROR!")
    with open("samples.txt", "wt") as samples:
        samples.write("\n".join(samp_list)+"\n")
