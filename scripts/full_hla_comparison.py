import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib.patches as mpatches

plt.rcParams["axes.grid"] = True 
plt.rcParams["grid.linestyle"] = "--"  
plt.rcParams["grid.alpha"] = 0.5  
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.facecolor"] = "#fafafa"

rec_palette = {
    "Fast_0": sns.color_palette("Oranges", 3)[0],
    "Fast_1": sns.color_palette("Oranges", 3)[1],
    "Fast_2": sns.color_palette("Oranges", 3)[2],
    "Seed_0": sns.color_palette("Blues", 3)[0],
    "Seed_1": sns.color_palette("Blues", 3)[1],
    "Seed_2": sns.color_palette("Blues", 3)[2]
}

version_patches = [
    mpatches.Patch(color=rec_palette["Fast_1"], label="Chain"),
    mpatches.Patch(color=rec_palette["Seed_1"], label="Seed"),
]
rec_patches = [
    mpatches.Patch(color="#f1f1f1", label="0 "),
    mpatches.Patch(color="lightgray", label="1 "),
    mpatches.Patch(color="darkgray", label="2 "),
]


genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

modes = ["0", "3", "5"]
recs = ["0", "1", "2"]

def parse_logs(directory):
    times = {}
    memory = {}
    edit_distance = {}
    for gene in genes:
        for rec in recs:
            for mode in modes:
                dir_path = f"{directory}/{gene}/{rec}/reads_{mode}_split/"
                count = 0
                count_gaf = 0
                for filename in os.listdir(dir_path):
                    if filename.endswith(".log"):
                        count += 1
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f:
                                if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                                    pattern = r"(?:(\d+):)?(\d{1,2}):(\d{2})(?:\.(\d+))?"

                                    # Cerca il tempo nella stringa
                                    match = re.search(pattern, line)

                                    if match:
                                        # Estrai ore, minuti, secondi e frazioni di secondo
                                        ore = int(match.group(1)) if match.group(1) else 0
                                        minuti = int(match.group(2)) if match.group(2) else 0
                                        secondi = int(match.group(3)) if match.group(3) else 0
                                        frazione_di_secondo = int(match.group(4)) if match.group(4) else 0
                                        
                                        # Converte tutto in millisecondi
                                        millisecondi_totali = (
                                            (ore * 3600 * 1000) +  # ore in millisecondi
                                            (minuti * 60 * 1000) +  # minuti in millisecondi
                                            (secondi * 1000) +  # secondi in millisecondi
                                            (frazione_di_secondo * (10 ** (3 - len(match.group(4)))) if match.group(4) else 0)  # frazioni di secondo
                                        )
                                        times[(gene, mode, rec, count)] = millisecondi_totali
                                if "Maximum resident set size (kbytes)" in line:
                                    memory[(gene, mode, rec, count)] = int(line.split()[-1])
                    if filename.endswith(".gaf"):
                        count_gaf += 1
                        if "ga" in directory:
                            with open(os.path.join(dir_path, filename), "r") as f:
                                for line in f.readlines():
                                    cigar = line.strip().split("\t")[-1].split(",")[0].split(":")[-1]
                                    edit_score = edit_distance_from_cigar(cigar)
                                    if "recombination" in line:
                                        edit_score += 4
                                edit_distance[(gene, mode, count_gaf)] = edit_score
                        else:
                            with open(os.path.join(dir_path, filename), "r") as f:
                                for line in f.readlines():
                                    if line.startswith("@CO"):
                                        edit_distance[(gene, mode, rec, count_gaf)] = int(line.split("\t")[1])
                    
    return times, memory, edit_distance

def parse_cigar(cigar_string):
    
    operations = []
    current_length = ""
    for char in cigar_string:
        if char.isdigit():
            current_length += char
        else:
            operations.append((int(current_length), char))
            current_length = ""
    return operations

def edit_distance_from_cigar(cigar_string):
    
    operations = parse_cigar(cigar_string)
    edit_distance = 0
    for length, op in operations:
        if op in 'IDX':
            edit_distance += length
    return edit_distance

old_times, old_memory, old_edit = parse_logs("output/HLA/ra_f")
new_times, new_memory, new_edit = parse_logs("output/HLA/ra_s")

data_time = []
data_memory = []

for gene in genes:
    for mode in modes:
        for rec in recs:
            for count in range(1, max(len(old_times), len(new_times)) + 1):
                if (gene, mode, rec, count) in old_times and (gene, mode, rec, count) in new_times:
                    data_time.append((gene, mode, rec, old_times[(gene, mode, rec, count)], "Fast"))
                    data_time.append((gene, mode, rec, new_times[(gene, mode, rec, count)], "Seed"))

                if (gene, mode,rec, count) in old_memory and (gene, mode, rec,count) in new_memory:
                    data_memory.append((gene, mode, rec, old_memory[(gene, mode, rec, count)], "Fast"))
                    data_memory.append((gene, mode, rec, new_memory[(gene, mode, rec, count)], "Seed"))


data_edit = []

for gene, mode, rec, count in old_edit:
    data_edit.append((gene, mode, rec, old_edit[(gene, mode, rec, count)], "Fast"))
    data_edit.append((gene, mode, rec, new_edit[(gene, mode, rec, count)], "Seed"))

df_edit = pd.DataFrame(data_edit, columns=["Gene", "Error", "Rec", "Edit", "Version"])

import math
                    
genes_part1 = genes[:10]
genes_part2 = genes[10:20]

# Numero di colonne e righe
ncols = 2
nrows = 5

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part1):
    df_edit_gene = df_edit[df_edit["Gene"] == gene]
    df_edit_gene["Version_Rec"] = df_edit_gene["Version"] + "_" + df_edit_gene["Rec"].astype(str)
    sns.boxplot(ax=axes[i], x="Error", y="Edit", hue="Version_Rec", data=df_edit_gene, palette=rec_palette, legend=False, linewidth=0, showfliers=False, boxprops=dict(alpha=.5))
    sns.boxplot(ax=axes[i], x="Error", y="Edit", hue="Version_Rec", data=df_edit_gene, palette=rec_palette, legend=False, fill=False)
    axes[i].set_yscale("symlog") 
    axes[i].set_ylim(-1, math.ceil(df_edit_gene["Edit"].max() / 1000) * 1000 + 1)
    axes[i].set_xlabel("Error", fontsize=16)
    axes[i].set_title(f"Gene: {gene}", fontsize=18)
    axes[i].set_ylabel("Edit distance", fontsize=16)

fig.legend(handles=version_patches, loc="lower center", ncols=2, bbox_to_anchor=(0.3, -0.0), title="Heuristic" , fontsize=18, title_fontsize=18, facecolor="white")
fig.legend(handles=rec_patches, loc="lower center", ncols=3, title="Recombination", bbox_to_anchor=(0.7, 0.0), fontsize=18, title_fontsize=18,facecolor="white")

plt.tight_layout()
plt.subplots_adjust(hspace=0.4, bottom=0.10)

plt.savefig("/home/dcmonti/Scaricati/fast_vs_seed_ed.png")
plt.show()

df_time = pd.DataFrame(data_time, columns=["Gene", "Error", "Rec", "Time", "Version"])
df_time = df_time[["Gene", "Version", "Error", "Rec", "Time"]]
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part1):
    df_time_gene = df_time[df_time["Gene"] == gene].copy()
    df_time_gene["Version_Rec"] = df_time_gene["Version"] + "_" + df_time_gene["Rec"].astype(str)
    
    sns.boxplot(ax=axes[i], x="Error", y="Time", hue="Version_Rec", data=df_time_gene, palette=rec_palette, legend=False, linewidth=0, showfliers=False, boxprops=dict(alpha=.5))
    sns.boxplot(ax=axes[i], x="Error", y="Time", hue="Version_Rec", data=df_time_gene, palette=rec_palette, legend=False, fill=False)
    axes[i].set_title(f"Gene: {gene}", fontsize=18)
    axes[i].set_yscale("symlog")
    y_min = math.floor(df_time_gene["Time"].min() / 1000) * 1000
    print((gene, y_min))
    axes[i].set_ylim(y_min -1, math.ceil(df_time_gene["Time"].max() / 1000) * 1000 + 1)
    axes[i].set_xlabel("Error", fontsize=16)
    axes[i].set_ylabel("Time (ms)", fontsize=16)

fig.legend(handles=version_patches, loc="lower center", ncols=2, bbox_to_anchor=(0.3, -0.0), title="Heuristic" , fontsize=18, title_fontsize=18, facecolor="white")
fig.legend(handles=rec_patches, loc="lower center", ncols=3, title="Recombination", bbox_to_anchor=(0.7, 0.0), fontsize=18, title_fontsize=18, facecolor="white")

plt.tight_layout()
plt.subplots_adjust(hspace=0.4, bottom=0.10)

plt.savefig("/home/dcmonti/Scaricati/fast_vs_seed_time.png")
plt.show()


