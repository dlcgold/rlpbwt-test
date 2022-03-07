configfile: "config.yaml"

import os
import math

#compile_cores = math.ceil(workflow.cores / 2)

data_folder = config["folder"]["data"]
software_folder = config["folder"]["software"]
results_folder = config["folder"]["res"]
input_folder = config["folder"]["input"]

pbwt_folder = os.path.join(software_folder, "pbwt")
htslib_folder = os.path.join(software_folder, "htslib")
rlpbwt_folder = os.path.join(software_folder, "rlpbwt")
results_time_folder = os.path.join(results_folder, "time")
results_parsed_time_folder = os.path.join(results_folder, "parsed_time")


query_number = config["params"]["query_number"]
panel_height = config["params"]["panel_height"]
panel_width = config["params"]["panel_width"]

new_panel_height = panel_height - query_number;

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
        cmake --build build -j1
        touch done.txt
        """

rule makeInputPbwt:
    input:
        os.path.join(pbwt_folder, "pbwt")
    output:
        os.path.join(input_pbwt_folder, "panel.pbwt"),
        os.path.join(input_pbwt_folder, "query.pbwt")
    shell:
        """
        cd {pbwt_folder}
        ./pbwt -readMacs ../../{input_panel} -write ../../{input_folder}/tmp.pbwt -writeSites ../../{input_folder}/tmp.sites
        ./pbwt -read ../../{input_folder}/tmp.pbwt -subsample 0 {new_panel_height} -write ../../{input_pbwt_folder}/panel.pbwt
        ./pbwt -read ../../{input_folder}/tmp.pbwt -subsample {new_panel_height} {query_number} -write ../../{input_pbwt_folder}/query.pbwt
        cd ../..
        rm {input_folder}/tmp.pbwt
        rm {input_folder}/tmp.sites
        """

rule makeInputRlpbwt:
    input:
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(input_rlpbwt_folder, "query.macs"),
        os.path.join(input_rlpbwt_folder, "panel.slp")
    shell:
        """
        python script/extract_query.py -i {input_panel} -o {input_rlpbwt_folder}/panel.macs -q  {input_rlpbwt_folder}/query.macs -n {query_number}
        cd {rlpbwt_folder}/script
        python build_slp.py -i ../../../{input_rlpbwt_folder}/panel.macs -o ../../../{input_rlpbwt_folder}/panel.slp
        """

rule runPbwtNaive:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        os.path.join(input_pbwt_folder, "panel.pbwt"),
        os.path.join(input_pbwt_folder, "query.pbwt")
    output:
        os.path.join(results_time_folder, "pbwtNaive.txt"),
        os.path.join(results_folder, "pbwtNaive.txt")
    shell:
        """
        cd {pbwt_folder}
        /usr/bin/time --verbose -o ../../{results_time_folder}/pbwtNaive.txt ./pbwt -read ../../{input_pbwt_folder}/panel.pbwt -matchNaive ../../{input_pbwt_folder}/query.pbwt > ../../{results_folder}/pbwtNaive.txt
        """
        
rule runPbwtIndexed:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        os.path.join(input_pbwt_folder, "panel.pbwt"),
        os.path.join(input_pbwt_folder, "query.pbwt"),
        os.path.join(pbwt_folder, "pbwt")
    output:
        os.path.join(results_time_folder, "pbwtIndexed.txt"),
        os.path.join(results_folder, "pbwtIndexed.txt")
    shell:
        """
        cd {pbwt_folder}
        /usr/bin/time --verbose -o ../../{results_time_folder}/pbwtIndexed.txt ./pbwt -read ../../{input_pbwt_folder}/panel.pbwt -matchIndexed ../../{input_pbwt_folder}/query.pbwt > ../../{results_folder}/pbwtIndexed.txt
        """

rule runPbwtDynamic:
    input:
        os.path.join(results_folder, "done.txt"),
        os.path.join(results_time_folder, "done.txt"),
        os.path.join(input_pbwt_folder, "panel.pbwt"),
        os.path.join(input_pbwt_folder, "query.pbwt"),
        os.path.join(pbwt_folder, "pbwt")
    output:
        os.path.join(results_time_folder, "pbwtDynamic.txt"),
        os.path.join(results_folder, "pbwtDynamic.txt")
    shell:
        """
        cd {pbwt_folder}
        /usr/bin/time --verbose -o ../../{results_time_folder}/pbwtDynamic.txt ./pbwt -read ../../{input_pbwt_folder}/panel.pbwt -matchDynamic ../../{input_pbwt_folder}/query.pbwt > ../../{results_folder}/pbwtDynamic.txt 
        """

rule makeRlpbwtNaive:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "naive.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -m ../../../{input_rlpbwt_folder}/naive.ser -N
        """
        
rule runRlpbwtNaive:
    input:
        os.path.join(input_rlpbwt_folder, "naive.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtNaive.txt"),
        os.path.join(results_folder, "rlpbwtNaive.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtNaive.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/naive.ser -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtNaive.txt -N
        """

rule makeRlpbwtBitvectors:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "bitvector.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -m ../../../{input_rlpbwt_folder}/bitvector.ser -B
        """

rule runRlpbwtBitvectors:
    input:
        os.path.join(input_rlpbwt_folder, "bitvector.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtBitvector.txt"),
        os.path.join(results_folder, "rlpbwtBitvector.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtBitvector.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/bitvector.ser -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtBitvector.txt -B
        """



rule makeRlpbwtPanel:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "panel.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -m ../../../{input_rlpbwt_folder}/panel.ser -P -t
        """

rule runRlpbwtPanel:
    input:
        os.path.join(input_rlpbwt_folder, "panel.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtPanel.txt"),
        os.path.join(results_folder, "rlpbwtPanel.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtPanel.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/panel.ser -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtPanel.txt -P -t
        """
        
rule makeRlpbwtPanelExt:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "panel_ext.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -m ../../../{input_rlpbwt_folder}/panel_ext.ser -P -t -e
        """

rule runRlpbwtPanelExt:
    input:
        os.path.join(input_rlpbwt_folder, "panel_ext.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtPanelExt.txt"),
        os.path.join(results_folder, "rlpbwtPanelExt.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtPanelExt.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/panel_ext.ser -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtPanelExt.txt -P -t -e
        """
        
rule makeRlpbwtSlpNoThr:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(input_rlpbwt_folder, "panel.slp"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "slp_no_thr.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -s ../../../{input_rlpbwt_folder}/panel.slp -m ../../../{input_rlpbwt_folder}/slp_no_thr.ser -S
        """

rule runRlpbwtSlpNoThr:
    input:
        os.path.join(input_rlpbwt_folder, "slp_no_thr.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtSlpNoThr.txt"),
        os.path.join(results_folder, "rlpbwtSlpNoThr.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtSlpNoThr.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/slp_no_thr.ser -s ../../../{input_rlpbwt_folder}/panel.slp  -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtSlpNoThr.txt -S
        """
        
rule makeRlpbwtSlpNoThrExt:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(input_rlpbwt_folder, "panel.slp"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "slp_no_thr_ext.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -s ../../../{input_rlpbwt_folder}/panel.slp -m ../../../{input_rlpbwt_folder}/slp_no_thr_ext.ser -S -e
        """

rule runRlpbwtSlpNoThrExt:
    input:
        os.path.join(input_rlpbwt_folder, "slp_no_thr_ext.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtSlpNoThrExt.txt"),
        os.path.join(results_folder, "rlpbwtSlpNoThrExt.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtSlpNoThrExt.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/slp_no_thr_ext.ser -s ../../../{input_rlpbwt_folder}/panel.slp  -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtSlpNoThrExt.txt -S -e
        """
        

rule makeRlpbwtSlpThr:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(input_rlpbwt_folder, "panel.slp"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "slp_thr.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -s ../../../{input_rlpbwt_folder}/panel.slp -m ../../../{input_rlpbwt_folder}/slp_thr.ser -S -t
        """

rule runRlpbwtSlpThr:
    input:
        os.path.join(input_rlpbwt_folder, "slp_thr.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtSlpThr.txt"),
        os.path.join(results_folder, "rlpbwtSlpThr.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtSlpThr.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/slp_thr.ser -s ../../../{input_rlpbwt_folder}/panel.slp  -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtSlpThr.txt -S -t
        """
        
        
rule makeRlpbwtSlpThrExt:
    input:
        os.path.join(input_rlpbwt_folder, "panel.macs"),
        os.path.join(input_rlpbwt_folder, "panel.slp"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(input_rlpbwt_folder, "slp_thr_ext.ser")
    shell:
        """
        cd {rlpbwt_folder}/build
        ./rlpbwt -i ../../../{input_rlpbwt_folder}/panel.macs -s ../../../{input_rlpbwt_folder}/panel.slp -m ../../../{input_rlpbwt_folder}/slp_thr_ext.ser -S -t -e
        """

rule runRlpbwtSlpThrExt:
    input:
        os.path.join(input_rlpbwt_folder, "slp_thr_ext.ser"),
        os.path.join(rlpbwt_folder, "done.txt")
    output:
        os.path.join(results_time_folder, "rlpbwtSlpThrExt.txt"),
        os.path.join(results_folder, "rlpbwtSlpThrExt.txt")
    shell:
        """
        cd {rlpbwt_folder}/build
        /usr/bin/time --verbose -o ../../../{results_time_folder}/rlpbwtSlpThrExt.txt ./rlpbwt -l ../../../{input_rlpbwt_folder}/slp_thr_ext.ser -s ../../../{input_rlpbwt_folder}/panel.slp  -q  ../../../{input_rlpbwt_folder}/query.macs -o ../../../{results_folder}/rlpbwtSlpThrExt.txt -S -e
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


