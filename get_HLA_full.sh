#!/bin/bash

input_file="genes_HLA_full.txt"
graphs_dir="https://raw.githubusercontent.com/ekg/HLA-zoo/master/graphs/spoa/"
seqs_dir="https://raw.githubusercontent.com/ekg/HLA-zoo/master/seqs/"

if [ ! -f "$input_file" ]; then
    echo "Error: File $input_file not exisiting."
    exit 1
fi

mkdir -p "output/HLA/genes"

while IFS= read -r gene_name

do
    gene_name_trimmed=$(echo "$gene_name" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//')
    graph_url="${graphs_dir}${gene_name_trimmed}.gfa"
    seq_url="${seqs_dir}${gene_name_trimmed}.fa"
    directory="$gene_name_trimmed"
    mkdir -p "output/HLA/genes/$directory"
    mkdir -p "output/HLA/genes/$directory/0"
    curl -o "output/HLA/genes/$directory/graph.gfa" "$graph_url"
    curl -o "output/HLA/genes/$directory/0/reads_0.fa" "$seq_url"
done < "$input_file"