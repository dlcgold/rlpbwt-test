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

data_tmp = [0,0,0,0,0,0,0,0,0]
chr_tmp = ""
for index, row in df_makes.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "pbwt":
        data_tmp[0] = int(round(row['wall_clock'],0))
    if row['variant'] == "naive":
        data_tmp[1] = int(round(row['wall_clock'],0))
    if row['variant'] == "bitvectors":
        data_tmp[2] = int(round(row['wall_clock'],0))
    if row['variant'] == "panel extended":
        data_tmp[3] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_thr extended":
        data_tmp[4] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_no_thr extended":
        data_tmp[5] = int(round(row['wall_clock'],0))
    if row['variant'] == "panel extended raw":
        data_tmp[6] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_thr_raw extended":
        data_tmp[7] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_no_thr_raw extended":
        data_tmp[8] = int(round(row['wall_clock'],0))
    if index == len(df_makes)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')

print()
print()
data_tmp = [0,0,0,0,0,0,0,0,0]
chr_tmp = ""
for index, row in df_makes.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "pbwt":
        data_tmp[0] = round(row['max_mem']*to_gb,1)
    if row['variant'] == "naive":
        data_tmp[1] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "bitvectors":
        data_tmp[2] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "panel extended":
        data_tmp[3] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "slp_thr extended":
        data_tmp[4] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "slp_no_thr extended":
        data_tmp[5] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "panel extended raw":
        data_tmp[6] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "slp_thr_raw extended":
        data_tmp[7] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "slp_no_thr_raw extended":
        data_tmp[8] = int(round(row['max_mem']*to_gb,0))
    if index == len(df_makes)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')



print()
print()
data_tmp = [0,0,0,0,0,0,0,0,0,0]
chr_tmp = ""
for index, row in df_exes.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "indexed":
        data_tmp[0] = int(round(row['wall_clock'],0))
    if row['variant'] == "dynamic":
        data_tmp[1] = int(round(row['wall_clock'],0))
    if row['variant'] == "naive":
        data_tmp[2] = int(round(row['wall_clock'],0))
    if row['variant'] == "bitvectors":
        data_tmp[3] = int(round(row['wall_clock'],0))
    if row['variant'] == "panel extended":
        data_tmp[4] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_thr extended":
        data_tmp[5] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_no_thr extended":
        data_tmp[6] = int(round(row['wall_clock'],0))
    if row['variant'] == "panel extended raw":
        data_tmp[7] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_thr_raw extended":
        data_tmp[8] = int(round(row['wall_clock'],0))
    if row['variant'] == "slp_no_thr_raw extended":
        data_tmp[9] = int(round(row['wall_clock'],0))
    if index == len(df_exes)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')

print()
print()
data_tmp = [0,0,0,0,0,0,0,0,0,0]
chr_tmp = ""
for index, row in df_exes.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "indexed":
        data_tmp[0] = round(row['max_mem']*to_gb,1)
    if row['variant'] == "dynamic":
        data_tmp[1] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "naive":
        data_tmp[2] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "bitvectors":
        data_tmp[3] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "panel extended":
        data_tmp[4] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "slp_thr extended":
        data_tmp[5] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "slp_no_thr extended":
        data_tmp[6] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "panel extended raw":
        data_tmp[7] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "slp_thr_raw extended":
        data_tmp[8] = round(row['max_mem']*to_gb,2)
    if row['variant'] == "slp_no_thr_raw extended":
        data_tmp[9] = round(row['max_mem']*to_gb,2)
    if index == len(df_exes)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')


print()
print()
data_tmp = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)]
chr_tmp = ""
for index, row in df_means.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(f"{elem[0]} $\pm$ {elem[1]}", end= ' & ')
                else:
                    print(f"{elem[0]} $\pm$ {elem[1]}", end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "indexed":
        data_tmp[0] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "dynamic":
        data_tmp[1] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "panel extended":
        data_tmp[2] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "slp_thr extended":
        data_tmp[3] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "slp_no_thr extended":
        data_tmp[4] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "panel extended raw":
        data_tmp[5] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "slp_thr_raw extended":
        data_tmp[6] = (round(row['mean'],2),round(row['stdev'],2))
    if row['variant'] == "slp_no_thr_raw extended":
        data_tmp[7] = (round(row['mean'],2),round(row['stdev'],2))
    if index == len(df_means)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(f"{elem[0]} $\pm$ {elem[1]}", end= ' & ')
            else:
                print(f"{elem[0]} $\pm$ {elem[1]}", end= ' ')


print()
print()
data_tmp = [0,0,0]
chr_tmp = ""
for index, row in df_pres.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "slp":
        data_tmp[0] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "convert":
        data_tmp[1] = int(round(row['max_mem']*to_gb,0))
    if row['variant'] == "query":
        data_tmp[2] = round(row['max_mem']*to_gb,5)
   
    if index == len(df_pres)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')


print()
print()
data_tmp = [0,0,0]
chr_tmp = ""
for index, row in df_pres.iterrows():
    if row["sites"] != chr_tmp:
        if chr_tmp != "":
            for elem in data_tmp:
                if data_tmp.index(elem) != len(data_tmp)-1:
                    print(elem, end= ' & ')
                else:
                     print(elem, end= ' ')
            print("\\\\")
        print(width_chr_map[row["sites"]], end=' & ') 
        chr_tmp=row["sites"]
    if row['variant'] == "slp":
        data_tmp[0] = int(round(row['wall_clock'],0))
    if row['variant'] == "convert":
        data_tmp[1] = int(round(row['wall_clock'],0))
    if row['variant'] == "query":
        data_tmp[2] = int(round(row['wall_clock'],0))
   
    if index == len(df_pres)-1:
        for elem in data_tmp:
            if data_tmp.index(elem) != len(data_tmp)-1:
                print(elem, end= ' & ')
            else:
                print(elem, end= ' ')