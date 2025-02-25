import os

haplos = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17','18','19','20']
for i in haplos:
    with open ('output/sars-cov-2/sim/filtered_sd_' + i + '_0001.maf', 'r') as f:
        for line in f:
            fields = line.strip().split(' ')
            fields = [x for x in fields if x != '']
            if len(fields) > 0 and fields[0] == 's' and fields[1] != 'ref':
                if fields[4] == '-':
                    fasta_name = fields[1]
                    try:
                        os.remove('output/sars-cov-2/reads/' + fasta_name + '.fasta')
                    except: 
                        pass

rec_haplos = ['01', '02', '03']
for i in rec_haplos:
    with open ('output/sars-cov-2/sim_rec/filtered_sd_' + i + '_0001.maf', 'r') as f:
        for line in f:
            fields = line.strip().split(' ')
            fields = [x for x in fields if x != '']
            if len(fields) > 0 and fields[0] == 's' and fields[1] != 'ref':
                if fields[4] == '-':
                    fasta_name = fields[1]
                    try:
                        os.remove('output/sars-cov-2/rec_reads/' + fasta_name + '.fasta')
                    except: 
                        pass