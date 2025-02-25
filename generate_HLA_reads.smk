import os

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("output/HLA/genes/{gene}/{rec}/reads_{err}.fa", 
        gene=genes,
        rec=[0,1,2],
        err=[0,3,5]
        ),
        

rule generate_rec_1:
    input:
        graph = "output/HLA/genes/{gene}/graph.gfa"
    output:
        results_file = "output/HLA/genes/{gene}/1/reads_0.fa"
    shell:
        "python scripts/generate_rec_reads.py {input.graph} 10 1 > {output.results_file}"

rule generate_rec_2:
    input:
        graph = "output/HLA/genes/{gene}/graph.gfa"
    output:
        results_file = "output/HLA/genes/{gene}/2/reads_0.fa"
    shell:
        "python scripts/generate_rec_reads.py {input.graph} 10 2 > {output.results_file}"

rule generate_3:
    input:
        reads_file = "output/HLA/genes/{gene}/{rec}/reads_0.fa"
    output:
        results_file = "output/HLA/genes/{gene}/{rec}/reads_3.fa"
    log:
        cigar = "output/HLA/genes/{gene}/{rec}/reads_3.log"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file} 2> {log.cigar}"
        
rule generate_5:
    input:
        reads_file = "output/HLA/genes/{gene}/{rec}/reads_0.fa"
    output:
        results_file = "output/HLA/genes{gene}/{rec}/reads_5.fa"
    log:
        cigar = "output/HLA/genes/{gene}/{rec}/reads_5.log"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file} 2> {log.cigar}"

