#!/bin/bash

dir=../data/fastq

for f_ in A01 A02 A03 A04
do
        printf "Calling peaks for %s\n" "${f_}"

        # finding ATAC fragment size d

        macs3 filterdup -i ${dir}/${f_}*_filtered.bam --keep-dup=1 -o ${dir}/MACS/${f_}_MACS_nodup.bam 2> /dev/null
        macs3 predictd -i ${dir}/MACS/${f_}_MACS_nodup.bam -g hs -m 5 50 2> ${dir}/MACS/frag_infer.log

        # extracting d from log file

        d=$(tail -n3 ${dir}/MACS/frag_infer.log | awk 'NR==1')
        IFS=' '
        read -a d <<< "$d"
        d=${d[12]}
        ss=$((d/4))
        es=$((ss+d+ss))
        printf "d was estimated to be %s\n" "$d"
        #rm ${dir}/MACS/${f_}_MACS_nodup.bam

        # peak calling using estimated d

        macs3 callpeak \
        -t ${dir}/${f_}*_filtered.bam \
        -n ${f_} \
        --outdir ${dir}/MACS/ \
        -f BAM \
        -g hs \
        -q 0.05 \
        --nomodel \
        --shift -${ss} \
        --extsize ${es} \
        --call-summits \
        -B 2> ${dir}/MACS/${f_}.log
done
