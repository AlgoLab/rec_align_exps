#!/bin/bash
 
haplos=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")
rec_haplos=("01" "02" "03")
for h in "${haplos[@]}"; do
    awk '/^>/ {out = "output/sars-cov-2/reads/" substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' output/sars-cov-2/reads/${h}/reads.fasta
    rm -rf output/sars-cov-2/reads/${h}
done

for h in "${rec_haplos[@]}"; do
    awk '/^>/ {out = "output/sars-cov-2/rec_reads/" substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' output/sars-cov-2/rec_reads/${h}/reads.fasta
    rm -rf output/sars-cov-2/rec_reads/${h}
done
python scripts/remove_rev_reads.py