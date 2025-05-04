import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import Levenshtein

genes = []
with open("genes_HLA_full.txt", "r") as f:
    for line in f:
        genes.append(line.strip())

recs = ["0", "1", "2"]
errs = ["1", "3", "5"]
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

def edit_distance_from_path(alignment_path, node_labels, query_sequence):
    path_sequence = build_path_label(alignment_path, node_labels)
    return Levenshtein.distance(path_sequence, query_sequence)
    
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

def extract_path_from_alignment(alignment_file, nodes_in_paths, query, nodes_label):
    with open(alignment_file, "r") as f:
        alignment = f.readline().strip()
        if alignment == "":
            print(alignment_file)
            return 0, len(query), False
        path = alignment.split("\t")[5]
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
        edit_score = edit_distance_from_path(nodes_with_directions, nodes_label, query)

        return swithces , edit_score, True

def compute_tools_performances(results):
    means_df = pd.DataFrame(columns=["tool", "gene", "mean_edit_score"])
    for gene in genes:
        results_gene = results[results["gene"] == gene]
        for tool in ["GraphAligner", "RecAlign", "Minichain"]:
            results_tool = results_gene[results_gene["tool"] == tool]
            mean_edit_score = round(results_tool["edit_score"].mean(),2)
            means_df = means_df._append({"tool": tool, "gene": gene, "mean_edit_score": mean_edit_score}, ignore_index=True)
    means_df.to_csv("output/HLA/mean_edit_scores.csv", index=False)
    return means_df

def compute_tools_performances2(results):
    means_df = pd.DataFrame(columns=["tool", "gene", "mean_time", "mean_memory"])
    for gene in genes:
        results_gene = results[results["gene"] == gene]
        for tool in ["GraphAligner", "RecAlign", "Minichain"]:
            results_tool = results_gene[results_gene["tool"] == tool]
            mean_time = round(results_tool["time"].sum(),2)
            mean_memory = round(results_tool["memory"].mean(),2)
            means_df = means_df._append({"tool": tool, "gene": gene, "mean_time": mean_time, "mean_memory": mean_memory}, ignore_index=True)
    means_df.to_csv("output/HLA/mean_time_mem.csv", index=False)
    return means_df

def extract_time_mem(log_file):
    time = 0.0
    memory = 0
    with open(log_file, "r") as f:
        for line in f:
            if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                time = line.split("):")[1].strip()
                time_parts = time.split(":")
                if len(time_parts) == 3:  # h:mm:ss
                    time = int(time_parts[0]) * 3600 + int(time_parts[1]) * 60 + float(time_parts[2])
                elif len(time_parts) == 2:  # m:ss
                    time = int(time_parts[0]) * 60 + float(time_parts[1])
            elif "Maximum resident set size (kbytes)" in line:
                memory = float(line.split(":")[1].strip())
    
    if time == 0.0:
        time = 0.001
    
    return time, memory

if __name__ == "__main__":
    results = pd.DataFrame(columns=["tool", "gene", "rec", "err", "read", "switches", "edit_score"])
    performances = pd.DataFrame(columns=["tool", "gene",  "rec", "err", "read","time", "memory"])
    not_aligned = pd.DataFrame(columns=["tool", "gene", "total_reads", "not_aligned_reads"])

    for gene in genes:
        not_aligned = not_aligned._append({"tool": "GraphAligner", "gene": gene, "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
        not_aligned = not_aligned._append({"tool": "RecAlign", "gene": gene, "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
        not_aligned = not_aligned._append({"tool": "Minichain", "gene": gene, "total_reads": 0, "not_aligned_reads": 0}, ignore_index=True)
        
        graph_file = f"output/HLA/genes/{gene}/graph.gfa"
        nodes_in_paths, nodes_label = extract_nodes_path(graph_file)

        for rec in recs:
            for err in errs:
                for read in os.listdir(f"output/HLA/genes/{gene}/{rec}/reads_{err}_split"):
                    if read.endswith(".fa"):
                        query = ""
                        with open(f"output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read}", "r") as f:
                            for line in f:
                                if not line.startswith(">"):
                                    query += line.strip()
                        read_id = re.search(r'read_(\d+)', read).group(1)
                        
                        switches, edit_score, aligned = extract_path_from_alignment(f"output/HLA/ga/{gene}/{rec}/reads_{err}_split/read_{read_id}.gaf", nodes_in_paths, query, nodes_label)
                        time, memory = extract_time_mem(f"output/HLA/ga/{gene}/{rec}/reads_{err}_split/read_{read_id}.log")
                        not_aligned.loc[(not_aligned["tool"] == "GraphAligner") & (not_aligned["gene"] == gene), "total_reads"] += 1
                        if aligned:
                            results = results._append({"tool": "GraphAligner","gene": gene, "rec": rec, "err": err, "read": read_id, "switches": switches, "edit_score": edit_score + 4*switches}, ignore_index=True)
                            performances = performances._append({"tool": "GraphAligner","gene": gene, "rec": rec, "err":err, "read": read_id, "time": time, "memory": memory}, ignore_index=True)
                        else:
                            not_aligned.loc[(not_aligned["tool"] == "GraphAligner") & (not_aligned["gene"] == gene), "not_aligned_reads"] += 1
                        
                        switches, edit_score, aligned = extract_path_from_alignment(f"output/HLA/ra_f/{gene}/{rec}/reads_{err}_split/read_{read_id}.gaf", nodes_in_paths, query, nodes_label)
                        time, memory = extract_time_mem(f"output/HLA/ra_f/{gene}/{rec}/reads_{err}_split/read_{read_id}.log")
                        not_aligned.loc[(not_aligned["tool"] == "RecAlign") & (not_aligned["gene"] == gene), "total_reads"] += 1
                        if aligned:
                            results = results._append({"tool": "RecAlign","gene": gene, "rec": rec, "err": err, "read": read_id, "switches": switches, "edit_score": edit_score + 4*switches}, ignore_index=True)
                            performances = performances._append({"tool": "RecAlign","gene": gene, "rec": rec, "err":err, "read": read_id, "time": time, "memory": memory}, ignore_index=True)
                        else:
                            not_aligned.loc[(not_aligned["tool"] == "RecAlign") & (not_aligned["gene"] == gene), "not_aligned_reads"] += 1
                            
                        switches, edit_score, aligned = extract_path_from_alignment(f"output/HLA/mc/{gene}/{rec}/reads_{err}_split/read_{read_id}.gaf", nodes_in_paths, query, nodes_label)
                        time, memory = extract_time_mem(f"output/HLA/mc/{gene}/{rec}/reads_{err}_split/read_{read_id}.log")
                        not_aligned.loc[(not_aligned["tool"] == "Minichain") & (not_aligned["gene"] == gene), "total_reads"] += 1
                        if aligned:
                            results = results._append({"tool": "Minichain","gene": gene, "rec": rec, "err": err, "read": read_id, "switches": switches, "edit_score": edit_score + 4*switches}, ignore_index=True)
                            performances = performances._append({"tool": "Minichain","gene": gene, "rec": rec, "err":err, "read": read_id, "time": time, "memory": memory}, ignore_index=True)
                        else:
                            not_aligned.loc[(not_aligned["tool"] == "Minichain") & (not_aligned["gene"] == gene), "not_aligned_reads"] += 1
    results.to_csv("output/HLA/switches_edit_scores.csv", index=False)
    performances.to_csv("output/HLA/time_mem.csv", index=False)
    not_aligned.to_csv("output/HLA/not_aligned.csv", index=False)
    means_perf = compute_tools_performances(results)
    means_perf.set_index(["tool", "gene"]).join(not_aligned.set_index(["tool", "gene"])).to_csv("output/HLA/final_perf.csv")
    means_time_mem = compute_tools_performances2(performances)