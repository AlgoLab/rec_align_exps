import os

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

error_levels = ["0", "1", "3", "5"]
recs= ["0", "1", "2"]

split_read_files = []
for gene in genes:
    for rec in recs:
        for err in error_levels:
            output_dir = f"output/HLA/genes/{gene}/{rec}/reads_{err}_split"
            if os.path.exists(output_dir):
                for read_file in os.listdir(output_dir):
                    if read_file.startswith("read_") and read_file.endswith(".fa"):
                        path = os.path.join(f"{gene}/{rec}/reads_{err}_split", read_file.split(".")[0])
                        split_read_files.append(path)

genes_short = []
with open("genes_HLA_short.txt", "r") as f:
    genes_short = [gene.strip() for gene in f.readlines()]

split_read_files_short = []
for gene in genes_short:
    for rec in recs:
        for err in error_levels:
            output_dir = f"output/HLA/genes/{gene}/{rec}/reads_{err}_split"
            if os.path.exists(output_dir):
                for read_file in os.listdir(output_dir):
                    if read_file.startswith("read_") and read_file.endswith(".fa"):
                        path = os.path.join(f"{gene}/{rec}/reads_{err}_split", read_file.split(".")[0])
                        split_read_files_short.append(path)
rule all:
    input:
        "bin/recgraph_a_star",
        "bin/minichain",
        "bin/recgraph",
        "bin/minigraph",
        expand("output/HLA/{mode}/{read_path}.gaf",
                mode=["ga", "ra_f","mc", "mg"], 
                read_path=split_read_files
        ),
        
        expand("output/HLA/rg/{read_path}.gaf",
                read_path=split_read_files_short
        )

rule build_minichain:
    output:
        "bin/minichain"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir minichain
        cd minichain
        git clone https://github.com/at-cg/minichain
        cd minichain && make
        cd ../..
        cp minichain/minichain/minichain {output}
        """

rule build_recgraph_a_star:
    output:
        "bin/recgraph_a_star"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir a_star
        cd a_star
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout a_star
        cargo build --release --jobs {threads}
        cd ../..
        cp a_star/RecGraph/target/release/recalign {output}
        """

rule build_recgraph:
    output:
        "bin/recgraph"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout main
        cargo build --release --jobs {threads}
        cd ..
        cp RecGraph/target/release/recgraph {output}
        """ 

rule build_minigraph:
    output:
        "bin/minigraph"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir minigraph
        cd minigraph
        git clone https://github.com/lh3/minigraph
        cd minigraph && make
        cd ../..
        cp minigraph/minigraph/minigraph {output}
        """

rule generate_mc_graphs:
    input:
        gfa="output/HLA/genes/{gene}/graph.gfa",
    output:
        w_gfa="output/HLA/genes/{gene}/w_graph.gfa"
    shell:
        "python scripts/mc_fix.py {input.gfa} > {output.w_gfa} "

rule run_recgraph_a_star_chaining:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/HLA/ra_c/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/ra_c/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 14 -r 4 -k 2 -e chaining> {output} 2> {log}"

rule run_recgraph_a_star_seeding:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/HLA/ra_s/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/ra_s/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 14 -r 4 -k 2 -e seeding> {output} 2> {log}"

rule run_recgraph_a_star_fast:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/HLA/ra_f/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/ra_f/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 14 -r 4 -k 2 -e fast> {output} 2> {log}"

rule run_minichain:
    input:
        minichain="bin/minichain",
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/w_graph.gfa",
    log:
        "output/HLA/mc/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/mc/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.minichain} -cx lr {input.gfa} {input.fa} -R 4 -k 12 -w 8 > {output} 2> {log}"

rule run_graphaligner:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
    log:
        "output/HLA/ga/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/ga/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v GraphAligner -f {input.fa}  -g {input.gfa} -a {output} -x vg 2> {log}"

rule run_minigraph:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        mg = "bin/minigraph"
    log:
        "output/HLA/mg/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/mg/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.mg} -cxsr {input.gfa}  {input.fa} > {output} 2> {log}"

rule run_recgraph:
    input:
        fa="output/HLA/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA/genes/{gene}/graph.gfa",
        rg = "bin/recgraph"
    log:
        "output/HLA/rg/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA/rg/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} {input.fa} {input.gfa} -m 8 > {output} 2> {log}"

