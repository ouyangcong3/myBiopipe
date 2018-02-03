#!/bin/bash

# function: assign wgs job for every sample 

if [[ ! -d "log" ]];then
    mkdir log
fi

cat samples.txt | while read samp; do
    nohup snakemake -npr -s wgs.py --config sample=${samp} 1>log/${samp}.log 2>&1 &
done
