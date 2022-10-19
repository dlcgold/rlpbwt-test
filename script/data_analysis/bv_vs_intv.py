from re import X
from turtle import color
from click import style
from matplotlib.lines import drawStyles
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import os
import statistics

plt.style.use('~/.matplotlib/stylelib/nord-light.mplstyle')
#plt.rc('xtick', labelsize=7.5) 

to_gb = 9.5367431640625e-7
b_to_gb = 1.1641532182693481e-10

def mystep(x,y, ax=None, where='post', **kwargs):
    assert where in ['post', 'pre']
    x = np.array(x)
    y = np.array(y)
    if where=='post': y_slice = y[:-1]
    if where=='pre': y_slice = y[1:]
    X = np.c_[x[:-1],x[1:],x[1:]]
    Y = np.c_[y_slice, y_slice, np.zeros_like(x[:-1])*np.nan]
    if not ax: ax=plt.gca()
    return ax.plot(X.flatten(), Y.flatten(), **kwargs)

x = np.arange(3,30000)
print(x)
run_prop = 12/4908
print(run_prop)
chr1= 1.1920928955078125e-7 * 1055454
y_bv = np.ceil(np.ceil(run_prop*x)*(2+np.log(np.ceil(x/np.ceil(run_prop*x)))))+128
y_iv = np.ceil(run_prop*x)*(np.ceil(np.log(x-2))+1)
#y_bv = np.ceil((np.ceil(np.ceil(run_prop*x)*(2+np.log(np.ceil(x/np.ceil(run_prop*x)))))+128)*chr1)
#y_iv = np.ceil(np.ceil(run_prop*x)*(np.ceil(np.log(x-2))+1)*chr1)
print(y_iv)
#dis=[]
#for i in range(len(y_iv)-1):
#    if np.abs(y_iv[i]- y_iv[i+1])>0.001:
#        #print(y_iv[i-1], y_iv[i], y_iv[i+1], y_iv[i+1])
#        dis.append((i))
#print(dis)
# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#plt.plot(x, y_bv, label="bitvector sparsi", color="#bf616a")
#plt.plot(x, y_iv, label="intvector compressi")
mystep(x, y_bv, label="bitvector sparsi", color="#bf616a")
mystep(x, y_iv, label="intvector compressi")

#y2 = np.ma.masked_where(((x>=dis[0])&(x<=dis[0]+1)|(x>=dis[1])&(x<=dis[1]+1)|(x>=dis[2])&(x<=dis[2]+1)|(x>=dis[3])&(x<=dis[3]+1)), y_iv) 
#plt.plot(x, y2, label="intvector compressi")
plt.axvline(x=17180, color='#a3be8c', linestyle='--', label="#sample = 17180")
#plt.axvline(x=16800, color='#a3be8c', linestyle='--', label="#aplotipi = 16800")

plt.legend(loc="upper left")
plt.xlabel("#sample", fontweight="bold")
plt.ylabel("Memory usage (bits)", fontweight="bold")
#plt.title("#sites = 1055454", fontweight="bold")
plt.savefig(f"../thesis_figures/bv_vs_iv.pdf", dpi=500, bbox_inches='tight')

# SLP vs MACS

## hardcoded 
macs=[(5201277377 * 0.00097656, "chr22 + 1055454"), (8571026807 * 0.00097656,  "chr20 + 1739315"), (10700222358 * 0.00097656,  "chr18 + 2171378"), (12792976272 * 0.00097656, "chr16 + 2596072"), (30537754694 * 0.00097656, "chr1 + 6196151")]
slps=[(46967317 * 0.00097656, "chr22 + 1055454"), (64899815 * 0.00097656,  "chr20 + 1739315"), (84058835 * 0.00097656,  "chr18 + 2171378"), (103522449 * 0.00097656, "chr16 + 2596072"), (237940983 * 0.00097656, "chr1 + 6196151")]

for i in range(len(macs)):
    print(f"{macs[i][1].split('+')[1].strip()} & 4908 & {round(slps[i][0]*to_gb,2)} & {round(macs[i][0]*to_gb,2)} & {round(float(1-((slps[i][0]*to_gb)/(macs[i][0]*to_gb)))*100,2)}\\\\")

plt.clf()
fig, ax = plt.subplots()
x_slp, y_slp = zip(*[(float(i[1].split('+')[1]),float(i[0])*to_gb ) for i in slps])
#plt.scatter(x_slp, y_slp)
plt.plot(x_slp, y_slp, label="SLP", marker="d")
x_macs, y_macs = zip(*[(float(i[1].split('+')[1]),float(i[0])*to_gb) for i in macs])
#plt.scatter(x_macs, y_macs)
plt.plot(x_macs, y_macs, label="MaCS", marker="p", color="#bf616a")

plt.title(f"SLP vs MaCS, log scale", fontweight="bold")
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_minor_formatter(NullFormatter())
plt.tick_params(axis='x', labelrotation = 45)
plt.xticks(x_slp)
plt.xlabel("#sites", fontweight="bold")
plt.ylabel("Maximum resident set size (gigabytes)", fontweight="bold")
plt.semilogy()
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig(f"../thesis_figures/slp_vs_macs_log.pdf", dpi=500)
