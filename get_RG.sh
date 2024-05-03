#!/bin/bash
fq=$1
header=$(zcat $fq | head -n 1 )
name=${fq%%.*}

IFS=':' read -ra fields <<< "$header"

instrument=$(echo "${fields[0]}" | sed 's/@//')
run=${fields[1]}
flowcell=${fields[2]}
lane=${fields[3]}
tile=${fields[4]}

read_group_id="${flowcell}.${lane}"

read_group="@RG\\tID:${read_group_id}\\tSM:${name}\\tLB:${instrument}\\tPL:Illumina"
echo "$read_group"
