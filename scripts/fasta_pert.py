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

        
    
def modifica_sequenza(sequenza, probabilita_modifica):
    nuova_sequenza = ""
    cigar = ""
    ed = 0
    for base in sequenza:
        if random.random() < probabilita_modifica:
            ed += 1
            nuova_base = random.choice(['A', 'C', 'G', 'T'])
            if nuova_base == base:
                nuova_base = random.choice(['', f"{base}{base}"])
                if nuova_base == '':
                    cigar += 'D'
                else:
                    cigar += 'I'
            else:
                cigar += 'X'
            nuova_sequenza += nuova_base   
        else:
            nuova_sequenza += base
            cigar += 'M'
    return nuova_sequenza, cigar, ed

def modifica_file_fasta(input_file, probabilita_modifica):
    with open(input_file, 'r') as f_input:
        cigar = ''
        ed = 0
        linea = f_input.readline()
        while linea:
            if linea.startswith('>'):  # Intestazione della sequenza
                print(linea, end='')
                print_cigar(cigar, ed)
                cigar = ''
                ed = 0
                linea = f_input.readline()
            else:  # Sequenza
                sequenza = linea.strip()
                sequenza_modificata,cigar_line, ed_update = modifica_sequenza(sequenza, probabilita_modifica)
                print(sequenza_modificata)
                cigar += cigar_line
                ed += ed_update
                linea = f_input.readline()
        print_cigar(cigar, ed)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modifica un file FASTA con una certa probabilità di mutazione.")
    parser.add_argument("input_file", type=str, help="Percorso del file FASTA di input.")
    parser.add_argument("--p", type=float, default=0.1, help="Probabilità di modifica della base [0.0 - 1.0]. Default: 0.1")
    args = parser.parse_args()

    modifica_file_fasta(args.input_file, args.p)

