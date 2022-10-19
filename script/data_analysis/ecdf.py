from cProfile import label
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import os
import statistics
import statsmodels.api as sm
from statsmodels.distributions.empirical_distribution import ECDF

plt.style.use('~/.matplotlib/stylelib/nord-light.mplstyle')
#plt.rc('xtick', labelsize=7.5) 
plt.rc('lines', linewidth=0.3)

runs = {}
directory = "../runs"
for filename in os.listdir(directory):
    file = os.path.join(directory, filename)
    with open(file) as f:
        tmp = file.split('/')[2].split('.')[0].replace("runs","")
        lines = f.readlines()
        tmp_runs = lines[0].split()
        tmp_runs.pop(-1)
        tmp_runs = [int(i) for i in tmp_runs]
        runs[filename.replace("runs","").split('.')[0]] = np.asarray(tmp_runs)
"""
#for key, value in runs.items():
#    print(key, "->", len(value), value[-1])

fig, (ax1) = plt.subplots(1, 1, figsize=(5,5))

ecdf = ECDF(runs['22'])
L1 = ax1.step(ecdf.x, ecdf.y, label="chr22")[0]

ecdf = ECDF(runs['20'])
L2 = ax1.step(ecdf.x, ecdf.y, label="chr20")[0]

ecdf = ECDF(runs['18'])
L3 = ax1.step(ecdf.x, ecdf.y, label="chr18")[0]

ecdf = ECDF(runs['16'])
L4 = ax1.step(ecdf.x, ecdf.y, label="chr16")[0]

ecdf = ECDF(runs['1'])
L5 = ax1.step(ecdf.x, ecdf.y, label="chr1")[0]

#plt.show()

#ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
#ax1.xaxis.set_minor_formatter(NullFormatter())
#ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"ECDF plot run distribution", fontweight="bold")
#ax1.set_xticks(x_makes)
ax1.set_xlabel("#Run", fontweight="bold")
ax1.set_ylabel("Probability", fontweight="bold")

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_single_time_dyn_paper.pdf", dpi=500)

line_labels = ["chr22", "chr20", "chr18", "chr16", "chr1"]
fig.legend(handles= [L1, L2, L3, L4, L5],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="center left",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(1, 0.5),
           ncol=1
           #title="Legend Title"  # Title for the legend
           )
plt.savefig(f"../thesis_figures/ecdf.pdf", dpi=500, bbox_inches='tight')

"""
plt.clf()
fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(10,5))

data = [runs['22'],runs['20'], runs['18'], runs['16'], runs['1']]
ax1.set_ylim([0, 20])
#for elem in data:
##    print(int(round(statistics.mean(elem),0)), " ", int(round(statistics.median(elem),0)), " ", int(round(statistics.mode(elem),0)))
L1 = ax1.boxplot(data, labels=["Chr22", "Chr20", "Chr18", "Chr16", "Chr1"],patch_artist = True, showmeans=True)
for median in L1['medians']:
    median.set(color='#bf616a',linewidth=1.5)
#ax1.axhline(y=3)
#ax1.axhline(y=13)

#plt.show()

#ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
#ax1.xaxis.set_minor_formatter(NullFormatter())
#ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"(a)\nBoxplot run distribution, max 20 runs", fontweight="bold")
#ax1.set_xticks(x_makes)
#ax1.set_xlabel("CHR", fontweight="bold")
#ax1.semilogy()

ax1.set_ylabel("#Run", fontweight="bold")

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.savefig(f"../thesis_figures/exe_single_time_dyn_paper.pdf", dpi=500)


#for elem in data:
##    print(int(round(statistics.mean(elem),0)), " ", int(round(statistics.median(elem),0)), " ", int(round(statistics.mode(elem),0)))
L1 = ax2.boxplot(data, labels=["Chr22", "Chr20", "Chr18", "Chr16", "Chr1"],patch_artist = True, showmeans=True)
for median in L1['medians']:
    median.set(color='#bf616a',linewidth=1.5)
#ax1.axhline(y=3)
#ax1.axhline(y=13)

#plt.show()

#ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
#ax1.xaxis.set_minor_formatter(NullFormatter())
#ax1.tick_params(axis='x', labelrotation = 45)
ax2.set_title(f"(b)\nBoxplot run distribution, full, log scale", fontweight="bold")
#ax1.set_xticks(x_makes)
#ax1.set_xlabel("CHR", fontweight="bold")
ax2.semilogy()

ax2.set_ylabel("#Run", fontweight="bold")

#plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.tight_layout()
#plt.savefig(f"../thesis_figures/exe_single_time_dyn_paper.pdf", dpi=500)

plt.savefig(f"../thesis_figures/boxplotbi.pdf", dpi=400, bbox_inches='tight')