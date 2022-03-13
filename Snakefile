configfile: "config.yaml"

import os
import math

compile_cores = math.ceil(workflow.cores / 2)
data_folder = config["folder"]["data"]
software_folder = config["folder"]["software"]
results_folder = config["folder"]["res"]
input_folder = config["folder"]["input"]

pbwt_folder = os.path.join(software_folder, "pbwt")
htslib_folder = os.path.join(software_folder, "htslib")
rlpbwt_folder = os.path.join(software_folder, "rlpbwt")
rlpbwt_build_folder = os.path.join(rlpbwt_folder, "build")

results_time_folder = os.path.join(results_folder, "time")
results_parsed_time_folder = os.path.join(results_folder, "parsed_time")


query_number = config["params"]["query_number"]
panel_height = config["params"]["panel_height"]
panel_width = config["params"]["panel_width"]

new_panel_height = panel_height - query_number;
new_panel_heightm = panel_height - query_number - 1;

input_name = config["params"]["input_panel"]
input_panel = os.path.join(input_folder, input_name)
input_pbwt_folder = os.path.join(input_folder, "pbwt")
input_rlpbwt_folder = os.path.join(input_folder, "rlpbwt")



rule makeSoftwareFolder:
    output:
        directory(software_folder),
        os.path.join(software_folder, "done.txt")
    shell:
        """
        mkdir -p {software_folder}
        touch {software_folder}/done.txt
        """

rule makeResultsFolder:
    output:
        directory(results_folder),
        os.path.join(results_folder, "done.txt"),
        directory(results_time_folder),
        os.path.join(results_time_folder, "done.txt"),
        directory(results_parsed_time_folder),
        os.path.join(results_parsed_time_folder, "done.txt")
    shell:
        """
        mkdir -p {results_folder}
        touch {results_folder}/done.txt
        mkdir -p {results_time_folder}
        touch {results_time_folder}/done.txt
        mkdir -p {results_parsed_time_folder}
        touch {results_parsed_time_folder}/done.txt
        """

rule downloadHtslib:
    input:
        os.path.join(software_folder, "done.txt")
    output:
        directory(htslib_folder),
        os.path.join(htslib_folder, "htsfile")
    shell:
        """
        cd {software_folder}
        git clone https://github.com/samtools/htslib
        cd htslib
        git submodule update --init --recursive
        autoreconf -i
        ./configure
        make
        """
        
rule downloadPbwt:
    input:
        os.path.join(htslib_folder, "htsfile")
    output:
        directory(pbwt_folder),
        os.path.join(pbwt_folder, "pbwt")
    shell:
        """
        cd {software_folder}
        git clone https://github.com/dlcgold/pbwt
        cd pbwt
        make
        """

rule downloadRlpbwt:
    input:
        os.path.join(software_folder, "done.txt")
    output:
        directory(rlpbwt_folder),
        os.path.join(rlpbwt_folder, "done.txt")
    shell:
        """
        cd {software_folder}
        git clone https://github.com/dlcgold/rlpbwt
        cd rlpbwt
        cmake -S . -B build -D BUILD_TESTS=OFF
        cmake --build build -j{compile_cores}
        touch done.txt
        """

rule makeInputPbwt:
    input:
        os.path.join(pbwt_folder, "pbwt")
    output:
        panel = os.path.join(input_pbwt_folder, "panel.pbwt"),
        query = os.path.join(input_pbwt_folder, "query.pbwt")
    shell:
        """
        ./{pbwt_folder}/pbwt -readMacs {input_panel} -write {input_folder}/tmp.pbwt -writeSites {input_folder}/tmp.sites
        ./{pbwt_folder}/pbwt -read {input_folder}/tmp.pbwt -subsample 0 {new_panel_heightm} -write {output.panel}
        ./{pbwt_folder}/pbwt -read {input_folder}/tmp.pbwt -subsample {new_panel_height} {query_number} -write {output.query}
        rm {input_folder}/tmp.pbwt
        rm {input_folder}/tmp.sites
        """

rule makeInputRlpbwt:
    input:
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
        query = os.path.join(input_rlpbwt_folder, "query.macs"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp")
    shell:
        """
        python script/extract_query.py -i {input_panel} -o {output.panel} -q  {output.query} -n {query_number}
        cd {rlpbwt_folder}/script
        python build_slp.py -i ../../../{output.panel} -o ../../../{output.slp}
        """

rule runPbwtNaive:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        panel = os.path.join(input_pbwt_folder, "panel.pbwt"),
        query = os.path.join(input_pbwt_folder, "query.pbwt")
    output:
        out_time = os.path.join(results_time_folder, "pbwtNaive.txt"),
        out = os.path.join(results_folder, "pbwtNaive.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time} ./{pbwt_folder}/pbwt -read {input.panel} -matchNaive {input.query} > {output.out}
        """
        
rule runPbwtIndexed:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        panel = os.path.join(input_pbwt_folder, "panel.pbwt"),
        query = os.path.join(input_pbwt_folder, "query.pbwt")
    output:
        out_time = os.path.join(results_time_folder, "pbwtIndexed.txt"),
        out = os.path.join(results_folder, "pbwtIndexed.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time} ./{pbwt_folder}/pbwt -read {input.panel} -matchIndexed {input.query} > {output.out}
        """

rule runPbwtDynamic:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        panel = os.path.join(input_pbwt_folder, "panel.pbwt"),
        query = os.path.join(input_pbwt_folder, "query.pbwt")
    output:
        out_time = os.path.join(results_time_folder, "pbwtDynamic.txt"),
        out = os.path.join(results_folder, "pbwtDynamic.txt")
    shell:
       """
        /usr/bin/time --verbose -o {output.out_time} ./{pbwt_folder}/pbwt -read {input.panel} -matchDynamic {input.query} > {output.out}
        """

rule makeRlpbwtNaive:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs")
    output:
        ser = os.path.join(input_rlpbwt_folder, "naive.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -m {output.ser} -N
        """
        
rule runRlpbwtNaive:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "naive.ser"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtNaive.txt"),
        out = os.path.join(results_folder, "rlpbwtNaive.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -q {input.query} -o {output.out} -N
        """

rule makeRlpbwtBitvectors:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs")
    output:
        ser = os.path.join(input_rlpbwt_folder, "bitvector.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -m {output.ser} -B
        """

rule runRlpbwtBitvectors:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "bitvector.ser"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtBitvector.txt"),
        out = os.path.join(results_folder, "rlpbwtBitvector.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -q {input.query} -o {output.out} -B
        """

rule makeRlpbwtPanel:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs")
    output:
        ser = os.path.join(input_rlpbwt_folder, "panel.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -m {output.ser} -P 
        """

rule runRlpbwtPanel:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "panel.ser"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtPanel.txt"),
        out = os.path.join(results_folder, "rlpbwtPanel.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -q {input.query} -o {output.out} -P
        """
        
rule makeRlpbwtPanelExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
    output:
        ser = os.path.join(input_rlpbwt_folder, "panel_ext.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -m {output.ser} -P -e
        """

rule runRlpbwtPanelExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "panel_ext.ser"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtPanelExt.txt"),
        out = os.path.join(results_folder, "rlpbwtPanelExt.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -q {input.query} -o {output.out} -P -e
        """
        
rule makeRlpbwtSlpNoThr:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp")
    output:
        ser = os.path.join(input_rlpbwt_folder, "slp_no_thr.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -s {input.slp} -m {output.ser} -S
        """

rule runRlpbwtSlpNoThr:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "slp_no_thr.ser"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtSlpNoThr.txt"),
        out = os.path.join(results_folder, "rlpbwtSlpNoThr.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output.out} -S
        """
        
rule makeRlpbwtSlpNoThrExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp")
    output:
        ser = os.path.join(input_rlpbwt_folder, "slp_no_thr_ext.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -s {input.slp} -m {output.ser} -S -e
        """

rule runRlpbwtSlpNoThrExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "slp_no_thr_ext.ser"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtSlpNoThrExt.txt"),
        out = os.path.join(results_folder, "rlpbwtSlpNoThrExt.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output.out} -S -e
        """
        

rule makeRlpbwtSlpThr:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp")
    output:
        ser = os.path.join(input_rlpbwt_folder, "slp_thr.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -s {input.slp} -m {output.ser} -S -t
        """
        
rule runRlpbwtSlpThr:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "slp_thr.ser"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtSlpThr.txt"),
        out = os.path.join(results_folder, "rlpbwtSlpThr.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output.out} -S -t
        """
        
rule makeRlpbwtSlpThrExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        panel = os.path.join(input_rlpbwt_folder, "panel.macs"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp")
    output:
        ser = os.path.join(input_rlpbwt_folder, "slp_thr_ext.ser")
    shell:
        """
        ./{rlpbwt_build_folder}/rlpbwt -i {input.panel} -s {input.slp} -m {output.ser} -S -t -e
        """

rule runRlpbwtSlpThrExt:
    input:
        os.path.join(rlpbwt_folder, "done.txt"),
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        rlpbwt = os.path.join(input_rlpbwt_folder, "slp_thr_ext.ser"),
        slp = os.path.join(input_rlpbwt_folder, "panel.slp"),
        query = os.path.join(input_rlpbwt_folder, "query.macs")
    output:
        out_time = os.path.join(results_time_folder, "rlpbwtSlpThrExt.txt"),
        out = os.path.join(results_folder, "rlpbwtSlpThrExt.txt")
    shell:
        """
        /usr/bin/time --verbose -o {output.out_time}  ./{rlpbwt_build_folder}/rlpbwt -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output.out} -S -t -e
        """
        
rule extractVerbose:
    input:
        "results/time/{file}.txt"
    output:
        "results/parsed_time/{file}.txt"
    shell:
        """
        python script/time_verbose_extractor.py -i {input} -o {output}
        """

rule analizeTime:
    input:
        os.path.join(results_parsed_time_folder, "pbwtNaive.txt"),
        os.path.join(results_parsed_time_folder, "pbwtIndexed.txt"),
        os.path.join(results_parsed_time_folder, "pbwtDynamic.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtNaive.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtBitvector.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtPanel.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtPanelExt.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtSlpNoThr.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtSlpNoThrExt.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtSlpThr.txt"),
        os.path.join(results_parsed_time_folder, "rlpbwtSlpThrExt.txt")
    output:
        os.path.join(results_folder, "time.png"),
        os.path.join(results_folder, "mem.png")
    shell:
        """
        python script/analyze_time_verbose.py {results_parsed_time_folder} {panel_height} {panel_width} {query_number}
        """
    
rule run_after_checkpoint:
    input:
        # os.path.join(htslib_folder),
        # os.path.join(pbwt_folder),
        # os.path.join(rlpbwt_folder),
        # os.path.join(input_pbwt_folder, "panel.pbwt"),
        # os.path.join(input_pbwt_folder, "query.pbwt"),
        # os.path.join(input_rlpbwt_folder, "panel.macs"),
        # os.path.join(input_rlpbwt_folder, "query.macs"),
        # os.path.join(input_rlpbwt_folder, "panel.slp"),
        # os.path.join(input_rlpbwt_folder, "naive.ser"),
        # os.path.join(input_rlpbwt_folder, "bitvector.ser"),
        # os.path.join(input_rlpbwt_folder, "panel.ser"),
        # os.path.join(input_rlpbwt_folder, "panel_ext.ser"),
        # os.path.join(input_rlpbwt_folder, "slp_no_thr.ser"),
        # os.path.join(input_rlpbwt_folder, "slp_no_thr_ext.ser"),
        # os.path.join(input_rlpbwt_folder, "slp_thr.ser"),
        # os.path.join(input_rlpbwt_folder, "slp_thr_ext.ser"),
        # os.path.join(results_time_folder, "pbwtNaive.txt"),
        # os.path.join(results_time_folder, "pbwtIndexed.txt"),
        # os.path.join(results_time_folder, "pbwtDynamic.txt"),
        # os.path.join(results_time_folder, "rlpbwtNaive.txt"),
        # os.path.join(results_time_folder, "rlpbwtBitvector.txt"),
        # os.path.join(results_time_folder, "rlpbwtPanel.txt"),
        # os.path.join(results_time_folder, "rlpbwtPanelExt.txt"),
        # os.path.join(results_time_folder, "rlpbwtSlpNoThr.txt"),
        # os.path.join(results_time_folder, "rlpbwtSlpNoThrExt.txt"),
        # os.path.join(results_time_folder, "rlpbwtSlpThr.txt"),
        # os.path.join(results_time_folder, "rlpbwtSlpThrExt.txt"),
        # os.path.join(results_parsed_time_folder, "pbwtNaive.txt"),
        # os.path.join(results_parsed_time_folder, "pbwtIndexed.txt"),
        # os.path.join(results_parsed_time_folder, "pbwtDynamic.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtNaive.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtBitvector.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtPanel.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtPanelExt.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtSlpNoThr.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtSlpNoThrExt.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtSlpThr.txt"),
        # os.path.join(results_parsed_time_folder, "rlpbwtSlpThrExt.txt"),
        os.path.join(results_folder, "time.png"),
        os.path.join(results_folder, "mem.png")


