
### Software
All the dependencies for this experiment are available on conda, so it is sufficient to create an environment with:
```
mamba create -c bioconda -n rg-exps --file requirements.lock 
```

### HLA experiment
In this experiment,  we downloaded the pangenome graphs for the HLA genes from [HLA-zoo]{https://github.com/ekg/HLA-zoo}, simulate recombinants haplotypes by traversing the graphs, and then align both recombinant and true haplotypse to their original graph using `RecAlign`, `RecGraph`, `GraphAligner` and `Minichain`.

```bash
# Download the selected genes from HLA-zoo
bash get_HLA_full.sh

# Generate and prepare the reads
snakemake -s generate_HLA_reads.smk -c 1
bash split_reads.sh

# Align the reads to their graphs
snakemake -s align_all_HLA.smk -c 16

# results are in output/HLA/{switches_edit_scores.csv, not_aligned.csv}
```

### Sars-Cov-2 experiment

In this experiment, we build graph using `pggb` using the haplotypes in `data/sars-cov-2/haplotypes`, we simulated reads using `pbsim3` both on the original ones and the recombinants in `data/sars-cov-2/rec_haplos` and then align each read using `RecAlign`, `Minigraph`, `GraphAligner` and `Minichain`.

```bash
# Prepare pbsim3 and simulates the reads
snakemake -s sars_cov2_reads_graphs.smk -c 1

# Align the reads
snakemake -s sars_cov2_align.smk -c 4

# results are in output/sars-cov-2/{switches_edit_scores.csv, not_aligned.csv}
```
