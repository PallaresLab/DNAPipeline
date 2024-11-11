#!/bin/bash
fq=$1
header=$(zcat $fq | head -n 1 )
name=$(basename $fq | sed 's/_R.*//')
IFS=':' read -ra fields <<< "$header"
id=$(echo "${fields[0]}" | sed 's/@//')

read_group="@RG\\tID:${id}\\tSM:${name}\\tLB:${name}\\tPL:Illumina"
echo "$read_group"
