import gfapy
import random
import sys

def extract_paths(file):
    with open(file) as f:
        lines = f.readlines()
        graph = gfapy.Gfa()
        paths = []
        sequences = {}
        for line in lines:
            line = line.strip()
            graph.append(line)
            if line.startswith("P"):
                paths.append(line)
            elif line.startswith("S"):
                fields = line.split("\t")
                idx = fields[1]
                sequence = fields[2]
                sequences[idx] = sequence
    return paths, graph, sequences

def parse_paths(paths):
    parsed_paths = []
    for path in paths:
        fields = path.split("\t")
        segments = fields[2].split(",")
        parsed_paths.append((fields[1], segments))
    return parsed_paths

def can_be_formed_by_other_paths(path, other_paths):
    for other_path1 in other_paths:
        for other_path2 in other_paths:
            if other_path1 != other_path2:
                for k in range(1, len(other_path1)):
                    for l in range(1, len(other_path2)):
                        if other_path1[:k] + other_path2[l:] == path:
                            return True
    return False

def find_common_nodes(path1, path2):
    return set(path1).intersection(set(path2))

def clean_node(node):
    return node.rstrip("+-")

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def generate_reads_from_paths_2(paths, sequences, max_reads=None):
    reads = []
    generated_reads = 0
    tries = 0
    while True:
        tries += 1
        path1_id, path_name1, segments1 = random.choice(paths)
        path2_id, path_name2, segments2 = random.choice(paths)
        
        if path1_id == path2_id or tries > 1000:
            continue
        
        common_nodes = find_common_nodes(segments1, segments2)
        
        if len(common_nodes) < 3:
            continue
        
        common_node1 = list(common_nodes)[len(common_nodes)//3]
        common_node2 = list(common_nodes)[len(common_nodes)//3*2]
        
        if common_node1 == common_node2:
            continue
        
        idx11 = segments1.index(common_node1)
        idx12 = segments1.index(common_node2)
        idx21 = segments2.index(common_node1)
        idx22 = segments2.index(common_node2)
        if idx21 > idx22:
            continue
        combined_segments = segments1[:idx11] + segments2[idx21:idx22] + segments1[idx12:]
        
        new_read = ''
        for node in combined_segments:
            cleaned_node = clean_node(node)
            if node.endswith('-'):
                new_read += reverse_complement(sequences[cleaned_node])
            else:
                new_read += sequences[cleaned_node]
                
        reads.append((path_name1, path_name2, (common_node1, common_node2), new_read))
        generated_reads += 1
        
        if max_reads is not None and generated_reads >= max_reads:
            break
            
    return reads

def generate_reads_from_paths(paths, sequences, max_reads=None):
    reads = []
    generated_reads = 0
    tries = 0
    while True:
        tries += 1
        path1_id, path_name1, segments1 = random.choice(paths)
        path2_id, path_name2, segments2 = random.choice(paths)
        
        if path1_id == path2_id or tries > 1000:
            continue
        
        common_nodes = find_common_nodes(segments1, segments2)
        
        if len(common_nodes) < 3:
            continue
        
        common_node = list(common_nodes)[len(common_nodes)//2]
        
        idx1 = segments1.index(common_node)
        idx2 = segments2.index(common_node)
        
        combined_segments = segments1[:idx1] + segments2[idx2:]
        
        new_read = ''
        for node in combined_segments:
            cleaned_node = clean_node(node)
            if node.endswith('-'):
                new_read += reverse_complement(sequences[cleaned_node])
            else:
                new_read += sequences[cleaned_node]
                
        reads.append((path_name1, path_name2, [common_node], new_read))
        generated_reads += 1
        if max_reads is not None and generated_reads >= max_reads:
            break
            
    return reads

def print_reads(reads):
    for i, read in enumerate(reads):
        rec_nodes = ",".join([node for node in read[2]])
        print(f">{read[0]},{read[1]},{rec_nodes}\n{read[3]}")

def process_gfa_file(input_file, rec_num, max_reads=None):
    paths, graph, sequences = extract_paths(input_file)
    parsed_paths = parse_paths(paths)
    
    unique_paths = []
    remove_paths = []

    seen = set()
    for i, (path_id, path) in enumerate(parsed_paths):
        path_tuple = tuple(path)
        if path_tuple not in seen:
            seen.add(path_tuple)
            unique_paths.append((i, path_id, path))
        else:
            remove_paths.append((i, path_id, path))
    
    #for i, (path_id, path) in enumerate(unique_paths):
    #    if can_be_formed_by_other_paths(path, [p[1] for p in unique_paths[:i] + unique_paths[i+1:]]):
    #        remove_paths.append((i, path))

    new_graph = []
    graph_lines = graph.lines
    p_count = 0
    for line in graph_lines:
        line = str(line).strip()
        if line.startswith("P"):
            if p_count not in [path[0] for path in remove_paths]:
                new_graph.append(line)
            p_count += 1
        else:
            new_graph.append(line)
    if rec_num == 2:
        reads = generate_reads_from_paths_2(unique_paths, sequences, max_reads=max_reads)
    if rec_num == 1:
        reads = generate_reads_from_paths(unique_paths, sequences, max_reads=max_reads)
    print_reads(reads)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_gfa_file> [max_reads]")
        sys.exit(1)

    input_file = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else None
    rec_num = int(sys.argv[3])
    process_gfa_file(input_file,rec_num, max_reads=max_reads)
