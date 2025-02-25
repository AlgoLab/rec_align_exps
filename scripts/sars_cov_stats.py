import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import Levenshtein
from Bio import Align


tools = ["ga", "mc", "ra", "mg"]
types = ["ori", "rec"]

def extract_nodes_path(graph_file):
    paths_in_nodes = {}
    nodes_label = {}
    with open(graph_file, "r") as f:
        for line in f:
            if line.startswith("P"):
                path_id = line.strip().split("\t")[1]
                path_nodes = line.strip().split("\t")[2].split(",")
                path_nodes = [node[:-1] for node in path_nodes]
                for node in path_nodes:
                    if node not in paths_in_nodes:
                        paths_in_nodes[node] = [path_id]
                    paths_in_nodes[node].append(path_id)
            elif line.startswith("S"):
                node_id = line.strip().split("\t")[1]
                node_label = line.strip().split("\t")[2]
                nodes_label[node_id] = node_label
    return paths_in_nodes, nodes_label

def edit_distance_from_path(alignment_path, node_labels, query_sequence, path_start, path_end):
    path_sequence = build_path_label(alignment_path, node_labels)
    sub_path = path_sequence[path_start:path_end]
    return Levenshtein.distance(sub_path, query_sequence)

def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(sequence))

def build_path_label(path, nodes_label):
    path_label = ""
    for node in path:
        if node[0] == ">":
            path_label += nodes_label[node[1:]]
        else:
            path_label += reverse_complement(nodes_label[node[1:]])
    return path_label

def extract_path_from_alignment(alignment_file, nodes_in_paths, query, nodes_label, tool):
    with open(alignment_file, "r") as f:
        alignment = f.readline().strip()
        if alignment == "":
            return 0, len(query), False, []
        path = alignment.split("\t")[5]
        path_start = int(alignment.split("\t")[7])
        if tool == "ra":
            path_end = int(alignment.split("\t")[8].strip())-1
        else:
            path_end = int(alignment.split("\t")[8].strip())
        path = ">".join(path.split("<"))
        alignment_nodes = path.split(">")[1:]
        alignment_paths = nodes_in_paths[alignment_nodes[0]]
        
        swithces = 0
        for node in alignment_nodes[1:]:
            if node in nodes_in_paths:
                common_paths = set(alignment_paths).intersection(set(nodes_in_paths[node]))
                if len(common_paths) == 0:
                    swithces += 1
                    alignment_paths = nodes_in_paths[node]
                else:
                    alignment_paths = common_paths
        nodes_with_directions = re.findall(r'[><]\d+', alignment.split("\t")[5])
        edit_score = edit_distance_from_path(nodes_with_directions, nodes_label, query, path_start, path_end)
        return swithces , edit_score, True, alignment_paths
    
def extract_time_mem(log_file):
    try:
        with open(log_file, "r") as f:
            for line in f:
                if "User time (seconds)" in line:
                    time = float(line.split(":")[1].strip())
                elif "Maximum resident set size (kbytes)" in line:
                    memory = float(line.split(":")[1].strip())
            return True, time, memory
    except FileNotFoundError:
        return False, 0, 0
        
if __name__ == "__main__":
    results = pd.DataFrame(columns=["tool", "type", "read", "switches", "edit_score", "time", "memory", "read_length", "paths"])    
    not_aligned = pd.DataFrame(columns=["tool", "total_reads", "not_aligned_reads"])

    not_aligned = not_aligned._append({"tool": "ga", "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
    not_aligned = not_aligned._append({"tool": "ra", "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
    not_aligned = not_aligned._append({"tool": "mc", "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
    
    graph_file = f"data/sars-cov-2/graph.gfa"
    nodes_in_paths, nodes_label = extract_nodes_path(graph_file)
    reads = os.listdir("output/sars-cov-2/reads")
    reads += os.listdir("output/sars-cov-2/rec_reads")
    for t in types:
        for tool in tools:
            for read in reads:
                read_id = read.split(".")[0]
                if read.endswith(".fasta"):
                    query = ""
                    
                    try:
                        if t == "ori":
                            read_folder = f"output/sars-cov-2/reads" 
                        else:
                            read_folder = f"output/sars-cov-2/rec_reads" 

                        with open(f"{read_folder}/{read_id}.fasta", "r") as f:
                            for line in f:
                                if not line.startswith(">"):
                                    query += line.strip()
                        
                        switches, edit_score, aligned, alignment_paths = extract_path_from_alignment(f"output/sars-cov-2/{t}/{tool}/{read_id}.gaf", nodes_in_paths, query, nodes_label, tool)
                        candidate_paths = ["None"] if switches > 0 else alignment_paths
                        _, time, memory = extract_time_mem(f"output/sars-cov-2/{t}/{tool}/{read_id}.log")
                        
                        not_aligned.loc[(not_aligned["tool"] == tool), "total_reads"] += 1
                        if aligned:
                            results = results._append({"tool": tool,"type": t, "read": read_id, "switches": switches, "edit_score": edit_score + 4*switches, "time": time, "memory": memory, "read_length": len(query), "paths": "\t".join(candidate_paths)}, ignore_index=True)
                        else:
                            not_aligned.loc[(not_aligned["tool"] == tool), "not_aligned_reads"] += 1
                        
                        
                            
                    except FileNotFoundError:
                        not_aligned.loc[(not_aligned["tool"] == tool), "total_reads"] += 1
                        not_aligned.loc[(not_aligned["tool"] == tool), "not_aligned_reads"] += 1
                    
    results.to_csv("output/sars-cov-2/switches_edit_scores.csv", index=False)
    not_aligned.to_csv("output/sars-cov-2/not_aligned.csv", index=False)
    