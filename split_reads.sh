#!/bin/bash

# Define the genes and error levels
GENES=$(cat genes_HLA_full.txt)
ERRORS=("0" "3" "5")
RECS=("0" "1" "2")

# Iterate over each gene and error level
for GENE in $GENES; do
    for REC in "${RECS[@]}"; do
        for ERR in "${ERRORS[@]}"; do
            INPUT_FILE="output/HLA/genes/${GENE}/${REC}/reads_${ERR}.fa"
            OUTPUT_DIR="output/HLA/genes/${GENE}/${REC}/reads_${ERR}_split"
            
            # Create the output directory if it doesn't exist
            mkdir -p $OUTPUT_DIR
            
            # Split the input file into individual read files
            awk -v OUT_DIR="$OUTPUT_DIR" '
            /^>/ {
                if (out_file) close(out_file)
                read_count++
                out_file = sprintf("%s/read_%d.fa", OUT_DIR, read_count)
            }
            { print > out_file }
            ' $INPUT_FILE
            rm -f $INPUT_FILE

            INPUT_LOG="output/HLA/genes/${GENE}/${REC}/reads_${ERR}.log"
            counter=0
            while IFS= read -r line; do
                output_file="${OUTPUT_DIR}/cigar_${counter}.log"
                echo "$line" > "$output_file"
                ((counter++))
            done < "$INPUT_LOG"
            rm -f $INPUT_LOG
            rm -f "${OUTPUT_DIR}/cigar_0.log"
        done
    done
done
