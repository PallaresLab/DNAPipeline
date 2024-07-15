#!/bin/bash
fq=$1
header=$(zcat $fq | head -n 1 )
name=$(echo $1 | cut -f1 -d"_")

IFS=':' read -ra fields <<< "$header"

id=$(echo "${fields[0]}" | sed 's/@//')
run=${fields[1]}
flowcell=${fields[2]}
lane=${fields[3]}
tile=${fields[4]}

read_group="@RG\\tID:${id}\\tSM:${name}\\tLB:${name}\\tPL:Illumina"
echo "$read_group"
