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

if __name__ == "__main__":
    import sys
    

    input_file = sys.argv[1]
    total_reads = int(sys.argv[2]) if len(sys.argv) > 2 else None
    reads = get_reads(input_file)
    if total_reads is not None:
        if len(reads) < total_reads:
            # duplicate reads until we reach the desired count
            original_reads = list(reads.items())
            while len(reads) < total_reads:
                for read_id, sequence in original_reads:
                    if len(reads) < total_reads:
                        reads[len(reads) + 1] = sequence
                    else:
                        break
    
    for read_id, sequence in reads.items():
        print(f">{read_id}\n{sequence}")