from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import os
import statistics

plt.style.use('~/.matplotlib/stylelib/nord-light.mplstyle')
#plt.rc('xtick', labelsize=7.5) 
plt.rc('lines', linewidth=0.5)


width_chr_map = {1055454: 'chr22', 1739315: 'chr20', 2171378: 'chr18', 2596072: 'chr16', 6196151: 'chr1',}
chr_runphi_map = {'chr16': (31187856, 31192761), 'chr18': (24288263, 24293169), 'chr20': (19966504, 19971410), 'chr22': (14772105, 14777009), 'chr1':(69671952, 69676858)}
height = 4908
smems = [118392, 141598, 170789, 228011, 495261]
df_exe = pd.read_csv('../data/time.csv')
df_make = pd.read_csv('../data/make-time.csv')
df_pre = pd.read_csv('../data/script.csv')
df_mean = pd.read_csv('../data/mean.csv')

to_gb = 9.5367431640625e-7

df_exe.insert(0, "chr", '')
for i in range(len(df_exe)):
    df_exe.loc[i, "chr"] = width_chr_map[df_exe.loc[i, "sites"]]
    
df_make.insert(0, "chr", '')
for i in range(len(df_make)):
    df_make.loc[i, "chr"] = width_chr_map[df_make.loc[i, "sites"]]

df_pre.insert(0, "chr", '')
for i in range(len(df_pre)):
    df_pre.loc[i, "chr"] = width_chr_map[df_pre.loc[i, "sites"]]

df_exe.insert(df_exe.shape[1],"run", 0)
for i in range(len(df_exe)):
    df_exe.loc[i, "run"] = chr_runphi_map[df_exe.loc[i, "chr"]][0]

df_make.insert(df_make.shape[1], "run", 0)
for i in range(len(df_make)):
    df_make.loc[i, "run"] = chr_runphi_map[df_make.loc[i, "chr"]][0]

df_pre.insert(df_pre.shape[1], "run", 0)
for i in range(len(df_pre)):
    df_pre.loc[i, "run"] = chr_runphi_map[df_pre.loc[i, "chr"]][0]

df_exe.insert(df_exe.shape[1], "phi", 0)
for i in range(len(df_exe)):
    df_exe.loc[i, "phi"] = chr_runphi_map[df_exe.loc[i, "chr"]][1]
    
df_make.insert(df_make.shape[1], "phi", 0)
for i in range(len(df_make)):
    df_make.loc[i, "phi"] = chr_runphi_map[df_make.loc[i, "chr"]][1]

df_pre.insert(df_pre.shape[1], "phi", 0)
for i in range(len(df_pre)):
    df_pre.loc[i, "phi"] = chr_runphi_map[df_pre.loc[i, "chr"]][1]

df_means = df_mean.sort_values(by='sites').reset_index(drop=True)
df_exes = df_exe.sort_values(by='sites').reset_index(drop=True)
df_makes = df_make.sort_values(by='sites').reset_index(drop=True)
df_pres = df_pre.sort_values(by='sites').reset_index(drop=True)

df_means = pd.read_csv('../data/mean.csv')



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))

#plt.figure(figsize=(9,6))   
#plt.subplot(1, 2, 1)
## make PBWT
tmp_df = df_makes[df_makes["variant"] == "pbwt"].reset_index(drop=True)
x_maked, y_maked = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maked, y_maked)
L1 = ax2.plot(x_maked, y_maked, label="PBWT", linewidth=1.5, marker="o")[0]


## make bitvector panel + threshold
tmp_df = df_makes[df_makes["variant"] == "panel extended"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L2 = ax2.plot(x_makep, y_makep, label="MAP-BV + THR-BV + RA-BV + PERM + PHI",linestyle='dashed',  marker="d")[0]

## make slp+threshold
tmp_df = df_makes[df_makes["variant"] == "slp_thr extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L3 = ax2.plot(x_maket, y_maket, label="MAP-BV + THR-BV + RA-SLP + PERM + PHI", linestyle='dashed', marker="s")[0]

## make slp+lce
tmp_df = df_makes[df_makes["variant"] == "slp_no_thr extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L4 = ax2.plot(x_makes, y_makes, label="MAP-BV + LCE + PERM + PHI",linestyle='dashed',  marker="p")[0]

## make bitvector panel + threshold
tmp_df = df_makes[df_makes["variant"] == "panel extended raw"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L6 = ax2.plot(x_makep, y_makep, label="MAP-INT + THR-INT + RA-BV + PERM + PHI", linewidth=1.5, marker="1")[0]

## make slp+threshold
tmp_df = df_makes[df_makes["variant"] == "slp_thr_raw extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L7 = ax2.plot(x_maket, y_maket, label="MAP-INT + THR-INT + RA-SLP + PERM + PHI", linestyle='dashed', marker="2")[0]

## make slp+lce
tmp_df = df_makes[df_makes["variant"] == "slp_no_thr_raw extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L8 = ax2.plot(x_makes, y_makes, label="MAP-INT + LCE + PERM + PHI",linewidth=1.5, marker="3")[0]

# ## make naive
# tmp_df = df_makes[df_makes["variant"] == "naive"].reset_index(drop=True)
# x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
# #plt.scatter(x_maken, y_maken)
# L9 = ax2.plot(x_maken, y_maken, label="MAP-INT + RLCP", linestyle='dashed', marker="P", color="#4c566a")[0]

# ## make bv
# tmp_df = df_makes[df_makes["variant"] == "bitvectors"].reset_index(drop=True)
# x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
# #plt.scatter(x_makeb, y_makeb)
# L10 = ax2.plot(x_makeb, y_makeb, label="MAP-BV + RLCP", linestyle='dashed', marker="X", color="#2e3440")[0]

ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.tick_params(axis='x', labelrotation = 45)
ax2.set_title(f"(b)\nConstruction time, log scale", fontweight="bold")
ax2.set_xticks(x_makes)
ax2.semilogy()
ax2.set_xlabel("#siti", fontweight="bold")
ax2.set_ylabel("Elapsed (wall clock) time (seconds)", fontweight="bold")
#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.tight_layout()
#plt.savefig(f"../thesis_figures/make_time_paper.pdf", dpi=500)

# MAX MEMORY FOR MAKE RLPBWT

#plt.clf()
#plt.subplot(1, 2, 2)
## make runs

#tmp_df = df_makes[df_makes["variant"] == "pbwt"].reset_index(drop=True)
#x_maker, y_maker = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "run"])) for i in range(len(tmp_df))])
#plt.scatter(x_maker, y_maker)
#plt.plot(x_maker, y_maker, label="O(r)")


## make PBWT
tmp_df = df_makes[df_makes["variant"] == "pbwt"].reset_index(drop=True)
x_maked, y_maked = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maked, y_maked)
L1 = ax1.plot(x_maked, y_maked, label="PBWT",linewidth=1.5, marker="o")[0]


## make bitvector panel + threshold
tmp_df = df_makes[df_makes["variant"] == "panel extended"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L2 = ax1.plot(x_makep, y_makep, label="MAP-BV + THR-BV + RA-BV + PERM + PHI", linestyle='dashed', marker="d")[0]

## make slp+threshold
tmp_df = df_makes[df_makes["variant"] == "slp_thr extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L3 = ax1.plot(x_maket, y_maket, label="MAP-BV + THR-BV + RA-SLP + PERM + PHI",linestyle='dashed',  marker="s")[0]

## make slp+lce
tmp_df = df_makes[df_makes["variant"] == "slp_no_thr extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L4 = ax1.plot(x_makes, y_makes, label="MAP-BV + LCE + PERM + PHI",linestyle='dashed',  marker="p")[0]

## make bitvector panel + threshold
tmp_df = df_makes[df_makes["variant"] == "panel extended raw"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L6 = ax1.plot(x_makep, y_makep, label="MAP-INT + THR-INT + RA-BV + PERM + PHI",linewidth=1.5, marker="1")[0]

## make slp+threshold
tmp_df = df_makes[df_makes["variant"] == "slp_thr_raw extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L7 = ax1.plot(x_maket, y_maket, label="MAP-INT + THR-INT + RA-SLP + PERM + PHI", linestyle='dashed', marker="2")[0]

## make slp+lce
tmp_df = df_makes[df_makes["variant"] == "slp_no_thr_raw extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L8 = ax1.plot(x_makes, y_makes, label="MAP-INT + LCE + PERM + PHI",linewidth=1.5, marker="3")[0]

# ## make naive
# tmp_df = df_makes[df_makes["variant"] == "naive"].reset_index(drop=True)
# x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
# #plt.scatter(x_maken, y_maken)
# L9 = ax1.plot(x_maken, y_maken, label="MAP-INT + RLCP", linestyle='dashed', marker="P", color="#4c566a")[0]

# ## make bv
# tmp_df = df_makes[df_makes["variant"] == "bitvectors"].reset_index(drop=True)
# x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
# #plt.scatter(x_makeb, y_makeb)
# L10 = ax1.plot(x_makeb, y_makeb, label="MAP-BV + RLCP", linestyle='dashed', marker="X", color = "#2e3440")[0]

ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"(a)\nConstruction max memory, log scale", fontweight="bold")
ax1.set_xticks(x_makes)
ax1.semilogy()

ax1.set_xlabel("#siti", fontweight="bold")
ax1.set_ylabel("Maximum resident set size (gigabytes)", fontweight="bold")
#plt.legend(bbox_to_anchor=(0,0), loc='lower center')
#fig.subplots_adjust(bottom=0.1, wspace=0.33)
plt.tight_layout()
line_labels = ["MAP-BV + THR-BV + RA-BV + PERM + PHI", "MAP-BV + THR-BV + RA-SLP + PERM + PHI", "MAP-BV + LCE + PERM + PHI", "MAP-INT + THR-INT + RA-BV + PERM + PHI", "MAP-INT + THR-INT + RA-SLP + PERM + PHI", "MAP-INT + LCE + PERM + PHI", "PBWT"]
fig.legend(handles= [L2, L3, L4, L6, L7, L8,L1],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="lower center",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(0.5, -0.13),
           ncol=3
           #title="Legend Title"  # Title for the legend
           )
plt.savefig(f"../thesis_figures/make_time_mem_paper2.pdf", dpi=500, bbox_inches='tight')


# #############################################################################

plt.clf()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
"""
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=0.5)
ax1 = plt.subplot(gs[0, :2], )
ax2 = plt.subplot(gs[0, 2:])
ax3 = plt.subplot(gs[1, 1:3])
"""
## run indexed
tmp_df = df_exes[df_exes["variant"] == "indexed"].reset_index(drop=True)
x_makei, y_makei = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makei, y_makei)
L1 = ax2.plot(x_makei, y_makei, label="PBWT MatchIndexed",linewidth=1.5, marker="o")[0]


## run bitvector panel + threshold
tmp_df = df_exes[df_exes["variant"] == "panel extended"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L2 = ax2.plot(x_makep, y_makep, label="MAP-BV + THR-BV + RA-BV + PERM + PHI",linestyle='dashed',  marker="d")[0]

## run slp+threshold
tmp_df = df_exes[df_exes["variant"] == "slp_thr extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L3 = ax2.plot(x_maket, y_maket, label="MAP-BV + THR-BV + RA-SLP + PERM + PHI", linestyle='dashed', marker="s")[0]

## run slp+lce
tmp_df = df_exes[df_exes["variant"] == "slp_no_thr extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L4 = ax2.plot(x_makes, y_makes, label="MAP-BV + LCE + PERM + PHI", linestyle='dashed', marker="p")[0]

# ## run dynamic
# tmp_df = df_exes[df_exes["variant"] == "dynamic"].reset_index(drop=True)
# x_maked, y_maked = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
# #plt.scatter(x_maked, y_maked)
# L5 = ax2.plot(x_maked, y_maked, label="PBWT MatchDynamic", marker="v")[0]


## run bitvector panel + threshold
tmp_df = df_exes[df_exes["variant"] == "panel extended raw"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L6 = ax2.plot(x_makep, y_makep, label="MAP-INT + THR-INT + RA-BV + PERM + PHI",linewidth=1.5, marker="1")[0]



## run slp+threshold
tmp_df = df_exes[df_exes["variant"] == "slp_thr_raw extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L7= ax2.plot(x_maket, y_maket, label="MAP-INT + THR-INT + RA-SLP + PERM + PHI",linestyle='dashed',  marker="2")[0]

## run slp+lce
tmp_df = df_exes[df_exes["variant"] == "slp_no_thr_raw extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L8 = ax2.plot(x_makes, y_makes, label="MAP-INT + LCE + PERM + PHI",linewidth=1.5, marker="3")[0]


# ## run naive
# tmp_df = df_exes[df_exes["variant"] == "naive"].reset_index(drop=True)
# x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
# #plt.scatter(x_maken, y_maken)
# L9 = ax2.plot(x_maken, y_maken, label="MAP-INT + RLCP", linestyle='dashed', marker="P", color="#4c566a")[0]

# ## run bv
# tmp_df = df_exes[df_exes["variant"] == "bitvectors"].reset_index(drop=True)
# x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
# #plt.scatter(x_makeb, y_makeb)
# L10 = ax2.plot(x_makeb, y_makeb, label="MAP-BV + RLCP", linestyle='dashed', marker="X", color = "#2e3440")[0]

ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.tick_params(axis='x', labelrotation = 45)
ax2.set_title(f"(b)\nSMEMs computation time \n100 queries, log scale", fontweight="bold")
ax2.set_xticks(x_makes)
ax2.set_xlabel("#siti", fontweight="bold")
ax2.set_ylabel("Elapsed (wall clock) time (seconds)", fontweight="bold")
ax2.semilogy()
#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_time_dyn_log_paper.pdf", dpi=500)



#plt.clf()

## extimation
## x_maked, y_maked = zip(*[(float(tmp_df.loc[i, "sites"]),0.00097656*13*height*float(tmp_df.loc[i, "sites"])) for i in range(len(tmp_df))])
## plt.scatter(x_maked, y_maked)
## plt.plot(x_maked, y_maked, label="13*N*M")

## run indexed
tmp_df = df_exes[df_exes["variant"] == "indexed"].reset_index(drop=True)
x_makei, y_makei = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makei, y_makei)
L1 = ax1.plot(x_makei, y_makei, label="PBWT MatchIndexed",linewidth=1.5, marker="o")[0]


## run bitvector panel + threshold
tmp_df = df_exes[df_exes["variant"] == "panel extended"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L2 = ax1.plot(x_makep, y_makep, label="MAP-BV + THR-BV + RA-BV + PERM + PHI", linestyle='dashed', marker="d")[0]

## run slp+threshold
tmp_df = df_exes[df_exes["variant"] == "slp_thr extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L3 = ax1.plot(x_maket, y_maket, label="MAP-BV + THR-BV + RA-SLP + PERM + PHI",linestyle='dashed',  marker="s")[0]

## run slp+lce
tmp_df = df_exes[df_exes["variant"] == "slp_no_thr extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L4 = ax1.plot(x_makes, y_makes, label="MAP-BV + LCE + PERM + PHI", linestyle='dashed', marker="p")[0]

# ## run dynamic
# tmp_df = df_exes[df_exes["variant"] == "dynamic"].reset_index(drop=True)
# x_maked, y_maked = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
# #plt.scatter(x_maked, y_maked)
# L5 = ax1.plot(x_maked, y_maked, label="PBWT MatchDynamic", marker="v")[0]

## run bitvector panel + threshold
tmp_df = df_exes[df_exes["variant"] == "panel extended raw"].reset_index(drop=True)
x_makep, y_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep)
L6 = ax1.plot(x_makep, y_makep, label="MAP-INT + THR-INT + RA-BV + PERM + PHI",linewidth=1.5, marker="1")[0]




## run slp+threshold
tmp_df = df_exes[df_exes["variant"] == "slp_thr_raw extended"].reset_index(drop=True)
x_maket, y_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket)
L7 = ax1.plot(x_maket, y_maket, label="MAP-INT + THR-INT + RA-SLP + PERM + PHI",linestyle='dashed',  marker="2")[0]

## run slp+lce
tmp_df = df_exes[df_exes["variant"] == "slp_no_thr_raw extended"].reset_index(drop=True)
x_makes, y_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes)
L8 = ax1.plot(x_makes, y_makes, label="MAP-INT + LCE + PERM + PHI",linewidth=1.5, marker="3")[0]

## run naive
# tmp_df = df_exes[df_exes["variant"] == "naive"].reset_index(drop=True)
# x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
# #plt.scatter(x_maken, y_maken)
# L9 = ax1.plot(x_maken, y_maken, label="MAP-INT + RLCP", linestyle='dashed', marker="P", color ="#4c566a" )[0]

# ## run bv
# tmp_df = df_exes[df_exes["variant"] == "bitvectors"].reset_index(drop=True)
# x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
# #plt.scatter(x_makeb, y_makeb)
# L10 = ax1.plot(x_makeb, y_makeb, label="MAP-BV + RLCP", linestyle='dashed', marker="X", color ="#2e3440")[0]

ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"(a)\nSMEMs computation max memory \n100 queries, log scale", fontweight="bold")
ax1.set_xticks(x_makes)
ax1.set_xlabel("#siti", fontweight="bold")
ax1.set_ylabel("Maximum resident set size (gigabytes)", fontweight="bold")
ax1.semilogy()

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_mem_dyn_log_paper.pdf", dpi=500)

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_single_time_dyn_paper.pdf", dpi=500)

line_labels = ["PBWT MatchIndexed", "MAP-BV + THR-BV + RA-BV + PERM + PHI", "MAP-BV + THR-BV + RA-SLP + PERM + PHI", "MAP-BV + LCE + PERM + PHI", "MAP-INT + THR-INT + RA-BV + PERM + PHI", "MAP-INT + THR-INT + RA-SLP + PERM + PHI", "MAP-INT + LCE + PERM + PHI"]
fig.legend(handles= [L1, L2, L3, L4, L6, L7, L8],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="lower center",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(0.5, -0.17),
           ncol=3
           #title="Legend Title"  # Title for the legend
           )
plt.savefig(f"../thesis_figures/exe_time_mem_paper2.pdf", dpi=500, bbox_inches='tight')



fig, (ax3) = plt.subplots(1, 1, figsize=(5,5))


tmp_df = df_means[df_means["variant"] == "indexed"].reset_index(drop=True)
x_makei, y_makei, e_makei = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_makei, y_makei, marker="o", label="PBWT MatchIndexed")
L1 = ax3.errorbar(x_makei, y_makei, yerr=e_makei, label="PBWT MatchIndexed", fmt="o", capsize=2, color='#5e81ac')[0]
L11 = ax3.plot(x_makei, y_makei,linewidth=1.5, color='#5e81ac')[0]


## run bitvector panel + threshold
tmp_df = df_means[df_means["variant"] == "panel extended"].reset_index(drop=True)
x_makep, y_makep, e_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep, marker="d", label="MAP-BV + THR-BV + RA-BV + PERM + PHI")
L2 = ax3.errorbar(x_makep, y_makep, yerr=e_makep, label="MAP-BV + THR-BV + RA-BV + PERM + PHI", fmt="d", capsize=2, color='#88c0d0')[0]
L21 = ax3.plot(x_makep, y_makep,linestyle='dashed',  color='#88c0d0')[0]

## run slp+threshold
tmp_df = df_means[df_means["variant"] == "slp_thr extended"].reset_index(drop=True)
x_maket, y_maket, e_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket, marker="s", label="MAP-BV + THR-BV + RA-SLP + PERM + PHI")
L3 = ax3.errorbar(x_maket, y_maket, yerr=e_maket, label="MAP-BV + THR-BV + RA-SLP + PERM + PHI", fmt="s", capsize=2, color='#a3be8c')[0]
L31 = ax3.plot(x_maket, y_maket,linestyle='dashed',  color='#a3be8c')[0]

## run slp+lce
tmp_df = df_means[df_means["variant"] == "slp_no_thr extended"].reset_index(drop=True)
x_makes, y_makes, e_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes, marker="p", label="MAP-BV + LCE + PERM + PHI")
L4 = ax3.errorbar(x_makes, y_makes, yerr=e_makes, label="MAP-BV + LCE + PERM + PHI", fmt="p", capsize=2, color='#bf616a')[0]
L41 = ax3.plot(x_makes, y_makes,linestyle='dashed',  color='#bf616a')[0]

# ## run dynamic
# tmp_df = df_means[df_means["variant"] == "dynamic"].reset_index(drop=True)
# x_maked, y_maked, e_maked = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
# #plt.scatter(x_maked, y_maked, marker="v", label="PBWT MatchDynamic")
# L5 = ax3.errorbar(x_maked, y_maked, yerr=e_maked, label="PBWT MatchDynamic", fmt="v", capsize=2, color='#ebcb8b')[0]
# L51 = ax3.plot(x_maked, y_maked, color='#ebcb8b')[0]


## run bitvector panel + threshold
tmp_df = df_means[df_means["variant"] == "panel extended raw"].reset_index(drop=True)
x_makep, y_makep, e_makep = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_makep, y_makep, marker="d", label="MAP-BV + THR-BV + RA-BV + PERM + PHI")
L6 = ax3.errorbar(x_makep, y_makep, yerr=e_makep, label="MAP-INT + THR-INT + RA-BV + PERM + PHI", fmt="1", capsize=2, color='#ebcb8b')[0]
L61 = ax3.plot(x_makep, y_makep,linewidth=1.5, color='#ebcb8b')[0]

## run slp+threshold
tmp_df = df_means[df_means["variant"] == "slp_thr_raw extended"].reset_index(drop=True)
x_maket, y_maket, e_maket = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_maket, y_maket, marker="s", label="MAP-BV + THR-BV + RA-SLP + PERM + PHI")
L7 = ax3.errorbar(x_maket, y_maket, yerr=e_maket, label="MAP-INT + THR-INT + RA-SLP + PERM + PHI", fmt="2", capsize=2, color='#b48ead')[0]
L71 = ax3.errorbar(x_maket, y_maket, linestyle='dashed', color='#b48ead',)[0]

## run slp+lce
tmp_df = df_means[df_means["variant"] == "slp_no_thr_raw extended"].reset_index(drop=True)
x_makes, y_makes, e_makes = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "mean"]), float(tmp_df.loc[i, "stdev"])) for i in range(len(tmp_df))])
#plt.scatter(x_makes, y_makes, marker="p", label="MAP-BV + LCE + PERM + PHI")
L8 = ax3.errorbar(x_makes, y_makes, yerr=e_makes, label="MAP-INT + LCE + PERM + PHI", fmt="3", capsize=2, color='#4c566a')[0]
L81 = ax3.errorbar(x_makes, y_makes,linewidth=1.5, color='#4c566a')[0]


ax3.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax3.xaxis.set_minor_formatter(NullFormatter())
ax3.tick_params(axis='x', labelrotation = 45)
ax3.set_title(f"Mean SMEMs computation time \n100 queries one at a time", fontweight="bold")
ax3.set_xticks(x_makes)
ax3.semilogy()
ax3.set_xlabel("#siti", fontweight="bold")
ax3.set_ylabel("Elapsed (wall clock) time (seconds)", fontweight="bold")

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_single_time_dyn_paper.pdf", dpi=500)

line_labels = ["PBWT MatchIndexed", "MAP-BV + THR-BV + RA-BV + PERM + PHI", "MAP-BV + THR-BV + RA-SLP + PERM + PHI", "MAP-BV + LCE + PERM + PHI", "MAP-INT + THR-INT + RA-BV + PERM + PHI", "MAP-INT + THR-INT + RA-SLP + PERM + PHI", "MAP-INT + LCE + PERM + PHI"]
fig.legend(handles= [L1, L2, L3, L4, L6, L7, L8],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="center left",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(1, 0.5),
           ncol=1
           #title="Legend Title"  # Title for the legend
           )
plt.savefig(f"../thesis_figures/exe_time_single_paper2.pdf", dpi=500, bbox_inches='tight')



## preprocess time
plt.clf()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))

## make SLP
tmp_df = df_pres[df_pres["variant"] == "slp"].reset_index(drop=True)
x_makei, y_makei = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makei, y_makei)
L1 = ax2.plot(x_makei, y_makei, label="Build SLP", marker="o")[0]

## extractr queries panel
tmp_df = df_pres[df_pres["variant"] == "query"].reset_index(drop=True)
x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_maken, y_maken)
L2 = ax2.plot(x_maken, y_maken, label="Extract 100 queries", marker="p")[0]

## convert VCF to MACs
tmp_df = df_pres[df_pres["variant"] == "convert"].reset_index(drop=True)
x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "wall_clock"])) for i in range(len(tmp_df))])
#plt.scatter(x_makeb, y_makeb)
L3 = ax2.plot(x_makeb, y_makeb, label="Convert VCF to MACs", marker="d")[0]


ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.tick_params(axis='x', labelrotation = 45)
ax2.set_title(f"(b)\nPreprocessing \ntime {height} samples, log scale", fontweight="bold")
ax2.set_xticks(x_makes)
ax2.set_xlabel("#siti", fontweight="bold")
ax2.set_ylabel("Elapsed (wall clock) time (seconds)", fontweight="bold")
ax2.semilogy()


## make SLP
tmp_df = df_pres[df_pres["variant"] == "slp"].reset_index(drop=True)
x_makei, y_makei = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makei, y_makei)
L1 = ax1.plot(x_makei, y_makei, label="Build SLP", marker="o")[0]

## extractr queries panel
tmp_df = df_pres[df_pres["variant"] == "query"].reset_index(drop=True)
x_maken, y_maken = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_maken, y_maken)
L2 = ax1.plot(x_maken, y_maken, label="Extract 100 queries", marker="p")[0]

## convert VCF to MACs
tmp_df = df_pres[df_pres["variant"] == "convert"].reset_index(drop=True)
x_makeb, y_makeb = zip(*[(float(tmp_df.loc[i, "sites"]),float(tmp_df.loc[i, "max_mem"] * to_gb)) for i in range(len(tmp_df))])
#plt.scatter(x_makeb, y_makeb)
L3 = ax1.plot(x_makeb, y_makeb, label="Convert VCF to MACs", marker="d")[0]



ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"(a)\nPreprocessing\nmax memory {height} samples, log scale", fontweight="bold")
ax1.set_xticks(x_makes)
ax1.set_xlabel("#siti", fontweight="bold")
ax1.set_ylabel("Maximum resident set size (gigabytes)", fontweight="bold")
ax1.semilogy()

line_labels = ["Build SLP", "Extract 100 queries", "Convert VCF to MACs"]
fig.legend(handles= [L1, L2, L3],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="lower center",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(0.5, -0.13),
           ncol=3
           #title="Legend Title"  # Title for the legend
           )
plt.savefig(f"../thesis_figures/prep_mem_time2.pdf", dpi=500, bbox_inches='tight')