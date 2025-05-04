import argparse
import random
import sys
import re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    
def print_cigar(seq, ed):
    pattern = r"(.)\1*"
    matches = re.finditer(pattern, seq)
    
    cigar_parts = []
    for match in matches:
        res = match.group()
        count = len(res) 
        cigar_parts.append(f"{count}{res[0]}") 
    cigar = ''.join(cigar_parts)
    eprint(f"{cigar}\t{ed}")
    
def get_reads(file):
    reads = {}
    with open(file) as f:
        i = 0
        for line in f:
            if line.startswith(">"):
                i += 1
                reads[i] = ""
                continue
            sequence = line.strip()
            reads[i] += sequence
    return reads

def add_err(sequence, error_rate):
    new_sequence = ""
    cigar = ""
    ed = 0
    for base in sequence:
        if random.random() < error_rate:
            ed += 1
            new_base = random.choice(['A', 'C', 'G', 'T'])
            if new_base == base:
                new_base = random.choice(['', f"{base}{base}"])
                if new_base == '':
                    cigar += 'D'
                else:
                    cigar += 'I'
            else:
                cigar += 'X'
            new_sequence += new_base   
        else:
            new_sequence += base
            cigar += 'M'
    return new_sequence, cigar, ed


if __name__ == "__main__":
    input_file = sys.argv[1]
    error_rate = float(sys.argv[2]) 
    cov = int(sys.argv[3]) 
    reads = get_reads(input_file)
    random.seed(1234)  
    for read_id, sequence in reads.items():
        for _ in range(cov):
            modified_sequence, cigar_line, ed = add_err(sequence, error_rate)
            print(f">{read_id}_{_+1}")
            print(modified_sequence)
            print_cigar(cigar_line, ed)
    
    