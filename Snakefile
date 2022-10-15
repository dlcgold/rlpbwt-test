configfile: "config.yaml"

import os

software_folder = config["folder"]["software"]
results_folder = config["folder"]["res"]
input_folder = config["folder"]["input"]
output_folder = config["folder"]["output"]

pbwt_folder = os.path.join(software_folder, "pbwt")
htslib_folder = os.path.join(software_folder, "htslib")
rlpbwt_folder = os.path.join(software_folder, "rlpbwt")
rlpbwt_build_folder = os.path.join(rlpbwt_folder, "build")

print(software_folder)
print(pbwt_folder)
print(htslib_folder)
print(rlpbwt_folder)
APPROACHES_PBWT = [
    "pbwtIndexed",
    "pbwtDynamic",
];
APPROACHES_RLPBWT = [
    "rlpbwtNaive",
    "rlpbwtBitvector",
    "rlpbwtPanelExt",
    "rlpbwtSlpNoThrExt",
    "rlpbwtSlpThrExt",
    "rlpbwtPanelExtRaw",
    "rlpbwtSlpNoThrExtRaw",
    "rlpbwtSlpThrExtRaw",
];
NAME_RLPBWT = [
    "naive",
    "bitvector",
    "panel_ext",
    "slp_no_thr_ext",
    "slp_thr_ext",
    "panel_ext_raw",
    "slp_no_thr_ext_raw",
    "slp_thr_ext_raw",
];

PANELS = [   
    "samples_30000-sites_4294/input.macs/q_100",
    "samples_30000-sites_4294/input.macs/q_50",
];


rule makeTime:
    input:
        os.path.join(results_folder, "make-time.csv"),
        os.path.join(results_folder, "time.csv"),
        #os.path.join(results_folder, "slp_vs_macs.png"),
        #os.path.join(results_folder, "pbwt_vs_rlpbwt.png"),
        #os.path.join(results_folder, "pbwt_vs_rlpbwt_dyn.png"),
    shell:
       """
       echo ok
       """

#rule compareTime:
#    input:
#        tmp = os.path.join(results_folder, "time_vs_mem.pdf"),
#    output:
#        os.path.join(results_folder, "slp_vs_macs.png"),
#        os.path.join(results_folder, "pbwt_vs_rlpbwt.png"),
#        os.path.join(results_folder, "pbwt_vs_rlpbwt_dyn.png"),
#    params:
#        out_dir = PANELS[0].split('/')[0]
#    shell:
#        """
#        python script/plot.py -f {params.out_dir}
#        """

#rule analyzeTime:
#    input:
#        os.path.join(results_folder, "time.csv"),
#    output:
#        os.path.join(results_folder, "time_vs_mem.pdf"),
#        os.path.join(results_folder, "time_vs_mem-loglog.pdf"),
#    params:
#        out_dir = os.path.join(results_folder)
#    conda: "envs/tidyverse.yml"
#    shell:
#        """
#        Rscript script/analyze_time_verbose.R {input} {params.out_dir}
#        """

rule mergeTime:
    input:
        expand(
            os.path.join(output_folder, "pbwt", "{panel}", "{approach}.time.csv"),
            approach = APPROACHES_PBWT,
            panel = PANELS,
        ),
        expand(
            os.path.join(output_folder, "rlpbwt", "{panel}", "{approach}.time.csv"),
            approach = APPROACHES_RLPBWT,
            panel = PANELS,
        ),
    output:
        os.path.join(results_folder, "time.csv"),
    conda: "envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule mergeTimeMake:
    input:
        expand(
            os.path.join(output_folder, "rlpbwt", "{panel}", "make-{approach}.time.csv"),
            approach = NAME_RLPBWT,
            panel = PANELS,
        ),
    output:
        os.path.join(results_folder, "make-time.csv"),
    conda: "envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule downloadPbwt:
    output:
        d=directory(pbwt_folder),
	h=directory(htslib_folder),
        exe=os.path.join(pbwt_folder, "pbwt")
    conda: "envs/compilation.yml"
    threads: max(workflow.cores * 0.5, 1)
    shell:
        """
        git clone https://github.com/samtools/htslib  {output.h}
        cd {output.h}
        git submodule update --init --recursive
        #autoreconf -i
        #./configure
        make
        cd ../..
        git clone https://github.com/dlcgold/pbwt {output.d}
        cd {output.d}
        #sed -i.bak 's/$(HTSLIB)/-lhts/' {output.d}/Makefile
        #sed -i.bak 's/CPPFLAGS=-I$(HTSDIR)/#CPPFLAGS=-I$(HTSDIR)/' {output.d}/Makefile
        make 
        """

rule downloadRlpbwt:
    output:
        d=directory(rlpbwt_folder),
        cm=os.path.join(rlpbwt_folder, "CMakeLists.txt")
    shell:
        """
        cd {software_folder}
        git clone https://github.com/dlcgold/rlpbwt 
        """

rule compileRlpbwt:
    input:
        d=rlpbwt_folder,
        cm=os.path.join(rlpbwt_folder, "CMakeLists.txt"),
    output:
        exe=os.path.join(rlpbwt_folder, "rlpbwt")
    conda: "envs/compilation.yml"
    threads: max(workflow.cores * 0.5, 1)
    shell:
        """
        cd {input.d}
        git checkout .
        rm -rf build
        cmake -S . -B build -D BUILD_TESTS=OFF
        cmake --build build -j{threads}
        cp build/rlpbwt .
        make -C build/_deps/bigrepair-src -j{threads}
        """

rule makeInputPbwt:
    input:
        exe = os.path.join(pbwt_folder, "pbwt"),
        panel = os.path.join(input_folder, "samples_{samples}-sites_{sites}", "input.macs")
    output:
        panel = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panel.pbwt"),
        query = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "query.pbwt"),
    shadow: "shallow"
    params:
        new_height=lambda wildcards: int(wildcards.samples) - int(wildcards.n_query)
    log:
        read_macs = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "readmacs.log"),
        panel = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panel.log"),
        query = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "query.log"),	
	timer = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "readmacstime.txt"),
	timep = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panelmacstime.txt"),
	timeq = os.path.join(output_folder, "pbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "querymacstime.txt"),
    shell:
        """
        /usr/bin/time --verbose -o {log.timer} {input.exe} -readMacs {input.panel} -write tmp.pbwt -writeSites tmp.sites &> {log.read_macs}
        /usr/bin/time --verbose -o {log.timep} {input.exe} -read tmp.pbwt -subsample 0 {params.new_height} -write {output.panel} &> {log.panel}
        /usr/bin/time --verbose -o {log.timeq} {input.exe} -read tmp.pbwt -subsample {params.new_height} {wildcards.n_query} -write {output.query} &> {log.query}
        """

rule makeInputRlpbwt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(input_folder, "samples_{samples}-sites_{sites}", "{panel}")
    output:
        panel = os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panel.macs"),
        query = os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "query.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panel.slp"),
    params:
        output_panel = lambda wildcards, output: os.path.abspath(output.panel),
        output_slp = lambda wildcards, output: os.path.abspath(output.slp),
    log:
        log = os.path.abspath(
            os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "panel.log")
        ),
        timequery = os.path.abspath(
            os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "query_extract_time.txt")
        ),
        timeslp = os.path.abspath(
            os.path.join(output_folder, "rlpbwt", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "slp_build_time.txt")
        ),
    threads: min(int(workflow.cores * 0.5) + 1, workflow.cores),
    shell:
        """
        /usr/bin/time --verbose -o {log.timequery} python script/extract_query.py -i {input.panel} -o {output.panel} -q  {output.query} -n {wildcards.n_query} &> {log.log}
        cd {rlpbwt_folder}/script
        /usr/bin/time --verbose -o {log.timeslp} python build_slp.py -i {params.output_panel} -o {params.output_slp} &>> {log.log}
        """

rule runPbwtIndexed:
    input:
        exe = os.path.join(pbwt_folder, "pbwt"),
        panel = os.path.join(output_folder, "pbwt", "{folder}", "panel.pbwt"),
        query = os.path.join(output_folder, "pbwt", "{folder}", "query.pbwt"),
    output:
        out = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexed.txt"),
        #outstat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexedStat.txt"),
    log:
        log = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexed.log"),
        time = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexed.time"),
        #logstat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexedStat.log"),
        #timestat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtIndexedStat.time"),
    shell:
        """
        #/usr/bin/time --verbose -o {log.timestat} {input.exe} -read {input.panel} -matchIndexedStat {input.query} > {output.outstat} 2> {log.logstat}
        /usr/bin/time --verbose -o {log.time} {input.exe} -read {input.panel} -matchIndexed {input.query} > {output.out} 2> {log.log}
        """

rule runPbwtDynamic:
    input:
        exe = os.path.join(pbwt_folder, "pbwt"),
        panel = os.path.join(output_folder, "pbwt", "{folder}", "panel.pbwt"),
        query = os.path.join(output_folder, "pbwt", "{folder}", "query.pbwt"),
    output:
        out = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamic.txt"),
        #outstat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamicStat.txt"),
    log:
        log = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamic.log"),
        time = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamic.time"),
        #logstat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamicStat.log"),
        #timestat = os.path.join(output_folder, "pbwt", "{folder}", "pbwtDynamicStat.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -read {input.panel} -matchDynamic {input.query} > {output.out} 2> {log.log}
        #/usr/bin/time --verbose -o {log.timestat} {input.exe} -read {input.panel} -matchDynamicStat {input.query} > {output.outstat} 2> {log.logstat}
        """

rule makeRlpbwtNaive:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "naive.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-naive.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-naive.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -m {output.ser} -N &> {log.log}
        """
        
rule runRlpbwtNaive:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "naive.ser"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtNaive.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtNaive.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtNaive.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -q {input.query} -o {output} -N &> {log.log}
        """

rule makeRlpbwtBitvectors:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "bitvector.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-bitvector.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-bitvector.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -m {output.ser} -B &> {log.log}
        """

rule runRlpbwtBitvectors:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "bitvector.ser"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtBitvector.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtBitvector.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtBitvector.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -q {input.query} -o {output} -B &> {log.log}
        """

rule makeRlpbwtPanel:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -m {output.ser} -P &> {log.log}
        """

rule runRlpbwtPanel:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.ser"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanel.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanel.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanel.time"),
    shell:
        """
        OMP_NUM_THREADS=1  /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -q {input.query} -o {output} -P &> {log.log}
        """

rule makeRlpbwtPanelExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "panel_ext.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel_ext.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel_ext.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -m {output.ser} -P -e &> {log.log}
        """

rule runRlpbwtPanelExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "panel_ext.ser"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExt.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExt.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExt.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -q {input.query} -o {output} -P -e &> {log.log}
        """
rule makeRlpbwtPanelExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "panel_ext_raw.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel_ext_raw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-panel_ext_raw.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -m {output.ser} -P -e -r&> {log.log}
        """

rule runRlpbwtPanelExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "panel_ext_raw.ser"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExtRaw.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExtRaw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtPanelExtRaw.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -q {input.query} -o {output} -P -e -r &> {log.log}
        """

rule makeRlpbwtSlpNoThr:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -v &> {log.log}
        """

rule runRlpbwtSlpNoThr:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThr.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThr.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThr.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S &> {log.log}
        """

rule makeRlpbwtSlpNoThrExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr_ext.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr_ext.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr_ext.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -e &> {log.log}
        """

rule runRlpbwtSlpNoThrExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr_ext.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExt.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExt.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExt.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S -e &> {log.log}
        """
        
rule makeRlpbwtSlpNoThrExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr_ext_raw.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr_ext_raw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_no_thr_ext_raw.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -e -r &> {log.log}
        """

rule runRlpbwtSlpNoThrExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_no_thr_ext_raw.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExtRaw.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExtRaw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpNoThrExtRaw.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S -e -r &> {log.log}
        """
        

rule makeRlpbwtSlpThr:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -t &> {log.log}
        """

rule runRlpbwtSlpThr:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThr.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThr.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThr.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S -t &> {log.log}
        """

rule makeRlpbwtSlpThrExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr_ext.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr_ext.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr_ext.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -t -e &> {log.log}
        """

rule runRlpbwtSlpThrExt:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr_ext.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExt.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExt.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExt.time"),
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S -t -e &> {log.log}
        """


rule makeRlpbwtSlpThrExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        panel = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.macs"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
    output:
        ser = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr_ext_raw.ser"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr_ext_raw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "make-slp_thr_ext_raw.time"),
    shell:
        """
        /usr/bin/time --verbose -o {log.time} {input.exe} -i {input.panel} -s {input.slp} -m {output.ser} -S -t -e -r &> {log.log}
        """

rule runRlpbwtSlpThrExtRaw:
    input:
        exe = os.path.join(rlpbwt_folder, "rlpbwt"),
        rlpbwt = os.path.join(output_folder, "rlpbwt", "{folder}", "slp_thr_ext_raw.ser"),
        slp = os.path.join(output_folder, "rlpbwt", "{folder}", "panel.slp"),
        query = os.path.join(output_folder, "rlpbwt", "{folder}", "query.macs"),
    output:
        os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExtRaw.txt"),
    log:
        log = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExtRaw.log"),
        time = os.path.join(output_folder, "rlpbwt", "{folder}", "rlpbwtSlpThrExtRaw.time"),
    shell:
        """
        OMP_NUM_THREADS=1  /usr/bin/time --verbose -o {log.time} {input.exe} -l {input.rlpbwt} -s {input.slp} -q {input.query} -o {output} -S -t -e -r &> {log.log}
        """

rule extractVerbose:
    input:
        os.path.join(output_folder, "{tool}", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "{file}.time")
    output:
        os.path.join(output_folder, "{tool}", "samples_{samples}-sites_{sites}", "{panel}", "q_{n_query}", "{file}.time.csv")
    shell:
        """
        python script/time_verbose_extractor.py {wildcards.tool} {wildcards.samples} {wildcards.sites} {wildcards.panel} {wildcards.n_query} < {input} > {output}
        """

