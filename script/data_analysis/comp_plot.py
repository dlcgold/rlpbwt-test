from math import log
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter

import numpy as np
import os
import statistics

plt.style.use('~/.matplotlib/stylelib/nord-light.mplstyle')
plt.rc('xtick', labelsize=7) 
plt.rc('lines', linewidth=0.5)


width_chr_map = {1055454: 'chr22', 1739315: 'chr20', 2171378: 'chr18', 2596072: 'chr16', 6196151: 'chr1',}
chr_runphi_map = {'chr16': (31187856, 31192761), 'chr18': (24288263, 24293169), 'chr20': (19966504, 19971410), 'chr22': (14772105, 14777009), 'chr1':(69671952, 69676858)}
height = 4908
to_mega = 0.00000095367432
df = pd.read_csv('../data/single_comp.csv')

plt.clf()


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))


x_makem, y_makem = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "p"])+float(df.loc[i, "uv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L1 = ax1.plot(x_makem, y_makem, label="MAP-INT", marker="1")[0]

x_makemb, y_makemb = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "rbv"])+float(df.loc[i, "ubv"])+float(df.loc[i, "vbv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L2 = ax1.plot(x_makemb, y_makemb, label="MAP-BV", marker="o")[0]

x_maket, y_maket = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "t"])) for i in range(len(df))])
L3 = ax1.plot(x_maket, y_maket, label="THR-INT", marker="2")[0]

x_maketb, y_maketb = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "tbv"])) for i in range(len(df))])
L4 = ax1.plot(x_maketb, y_maketb, label="THR-BV", marker="s")[0]

x_maker, y_maker = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "panel"])) for i in range(len(df))])
L5 = ax1.plot(x_maker, y_maker, label="RA-BV", marker="p")[0]

x_makes, y_makes = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "slp"])) for i in range(len(df))])
L6 = ax1.plot(x_makes, y_makes, label="RA-SLP/LCE", marker="v")[0]

x_makep, y_makep = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "pref"])) for i in range(len(df))])
L7 = ax1.plot(x_makep, y_makep, label="PERM", marker="d")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "phis"])) for i in range(len(df))])
L8 = ax1.plot(x_makee, y_makee, label="PHI", marker="X")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "lcp"])) for i in range(len(df))])
L9 = ax1.plot(x_makee, y_makee, label="RLCP", marker="P")[0]

print(x_makem)
print(x_makemb)
print(x_maket)
print(x_maketb)
print(x_maker)
print(x_makes)
print(x_makep)
print(x_makee)
print()

ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"(a)\nComponents size by sites, log scale", fontweight="bold")
ax1.set_xticks(x_makes)
ax1.semilogy()
ax1.set_xlabel("#sites", fontweight="bold")
ax1.set_ylabel("Memory usage (megabytes)", fontweight="bold")


x_makem, y_makem = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "p"])+float(df.loc[i, "uv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L1 = ax2.plot(x_makem, y_makem, label="MAP-INT", marker="1")[0]

x_makemb, y_makemb = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "rbv"])+float(df.loc[i, "ubv"])+float(df.loc[i, "vbv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L2 = ax2.plot(x_makemb, y_makemb, label="MAP-BV", marker="o")[0]

x_maket, y_maket = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "t"])) for i in range(len(df))])
L3 = ax2.plot(x_maket, y_maket, label="THR-INT", marker="2")[0]

x_maketb, y_maketb = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "tbv"])) for i in range(len(df))])
L4 = ax2.plot(x_maketb, y_maketb, label="THR-BV", marker="s")[0]

x_maker, y_maker = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "panel"])) for i in range(len(df))])
L5 = ax2.plot(x_maker, y_maker, label="RA-BV", marker="p")[0]

x_makes, y_makes = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "slp"])) for i in range(len(df))])
L6 = ax2.plot(x_makes, y_makes, label="RA-SLP/LCE", marker="v")[0]

x_makep, y_makep = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "pref"])) for i in range(len(df))])
L7 = ax2.plot(x_makep, y_makep, label="PERM", marker="d")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "phis"])) for i in range(len(df))])
L8 = ax2.plot(x_makee, y_makee, label="PHI", marker="X")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Runs"]),float(df.loc[i, "lcp"])) for i in range(len(df))])
L9 = ax2.plot(x_makee, y_makee, label="RLCP", marker="P")[0]

print(x_makem)
print(x_makemb)
print(x_maket)
print(x_maketb)
print(x_maker)
print(x_makes)
print(x_makep)
print(x_makee)


ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.tick_params(axis='x', labelrotation = 45)
ax2.set_title(f"(b)\nComponents size by runs, log scale", fontweight="bold")
ax2.set_xticks(x_makes)
ax2.semilogy()
ax2.set_xlabel("#runs", fontweight="bold")
ax2.set_ylabel("Memory usage (megabytes)", fontweight="bold")


plt.tight_layout()
line_labels = ["MAP-INT", "MAP-BV", "THR-INT", "THR-BV", "RA-BV", "RA-SLP/LCE", "PERM", "PHI", "RLCP"]
fig.legend(handles= [L1, L2, L3, L4, L5, L6, L7, L8, L9],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="lower center",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(0.5, -0.06),
           ncol=9
           #title="Legend Title"  # Title for the legend
           )
           
plt.savefig(f"../thesis_figures/comp_mem.pdf", dpi=500, bbox_inches='tight')

print()
for idx, row in df.iterrows():
    print(f"{int(row['Sites'])} & {height} & {row['slp']} & {row['panel']} & {round((row['slp']/row['panel'])*100, 2)}\\\\")




plt.clf()


fig, (ax1) = plt.subplots(1, 1, figsize=(5,5))


x_makem, y_makem = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "p"])+float(df.loc[i, "uv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L1 = ax1.plot(x_makem, y_makem, label="MAP-INT", marker="1")[0]

x_makemb, y_makemb = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "rbv"])+float(df.loc[i, "ubv"])+float(df.loc[i, "vbv"])+float(df.loc[i, "Sites"])*4*to_mega+float(df.loc[i, "Sites"])*to_mega) for i in range(len(df))])
L2 = ax1.plot(x_makemb, y_makemb, label="MAP-BV", marker="o")[0]

x_maket, y_maket = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "t"])) for i in range(len(df))])
L3 = ax1.plot(x_maket, y_maket, label="THR-INT", marker="2")[0]

x_maketb, y_maketb = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "tbv"])) for i in range(len(df))])
L4 = ax1.plot(x_maketb, y_maketb, label="THR-BV", marker="s")[0]

x_maker, y_maker = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "panel"])) for i in range(len(df))])
L5 = ax1.plot(x_maker, y_maker, label="RA-BV", marker="p")[0]

x_makes, y_makes = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "slp"])) for i in range(len(df))])
L6 = ax1.plot(x_makes, y_makes, label="RA-SLP/LCE", marker="v")[0]

x_makep, y_makep = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "pref"])) for i in range(len(df))])
L7 = ax1.plot(x_makep, y_makep, label="PERM", marker="d")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "phis"])) for i in range(len(df))])
L8 = ax1.plot(x_makee, y_makee, label="PHI", marker="X")[0]

x_makee, y_makee = zip(*[(float(df.loc[i, "Sites"]),float(df.loc[i, "lcp"])) for i in range(len(df))])
L9 = ax1.plot(x_makee, y_makee, label="RLCP",  marker="P")[0]


ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.tick_params(axis='x', labelrotation = 45)
ax1.set_title(f"Components size, log scale", fontweight="bold")
ax1.set_xticks(x_makes)
ax1.semilogy()
ax1.set_xlabel("#sites", fontweight="bold")
ax1.set_ylabel("Memory usage (megabytes)", fontweight="bold")

plt.tight_layout()
line_labels = ["MAP-INT", "MAP-BV", "THR-INT", "THR-BV", "RA-BV", "RA-SLP/LCE", "PERM", "PHI", "RLCP"]
fig.legend(handles= [L1, L2, L3, L4, L5, L6, L7, L8, L9],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="center left",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           bbox_to_anchor=(1, 0.5),
           ncol=1
           #title="Legend Title"  # Title for the legend
           )
           
plt.savefig(f"../thesis_figures/comp_mem2.pdf", dpi=500, bbox_inches='tight')