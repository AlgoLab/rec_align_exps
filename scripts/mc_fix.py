import re
import sys

graph_path = sys.argv[1]
counter = 0 
new_graph = []
for line in open(graph_path):
    if line.startswith("P"):
        counter +=1
        line = line.strip()
        _, name, path, _ = line.split("\t")
        nodes =re.split("[+-],", path)
        nodes[-1] = nodes[-1][:-1]
        directions = re.split("\d+", path)
        directions.pop(0)
        directions = ['>' if '+' in dir else '<' for dir in directions]
        updated_path = ''.join([x + y for x, y in zip(directions, nodes)])
        walk = f"W\t{name}\t1\t{counter}\t.\t.\t{updated_path}\n"
        new_graph.append(walk)
    else:
        new_graph.append(line)
new_graph = ''.join(new_graph)
print(new_graph)
#for line in open(graph_path):
#    if line.startswith("W"):
#        line = line.strip()
#        _, name, hapix, seqid, seqstart, seqend, walk = line.split()
#        # seqend = 1
#        # seqstart = 1
#        print("W", f"{name}#{hapix}_W", 1, seqid, seqstart, seqend, walk, sep="\t")
#    else:
#        print(line, end="")