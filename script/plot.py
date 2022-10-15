import sys
import getopt
import os
import math
import matplotlib.pyplot as plt

class data:
    def __init__(self, height, width, slpsize, macssize, naivesize, bvsize, panelsize, msthrsize, mslcesize):
        self.height = height
        self.width = width
        self.slpsize = slpsize
        self.macssize = macssize
        self.naivesize = naivesize
        self.bvsize = bvsize
        self.panelsize = panelsize
        self.msthrsize = msthrsize
        self.mslcesize = mslcesize
        self.dynsize = (self.height + self.width) * 4 * 0.00097656
        self.durbin = self.height * self.width * 13 * 0.00097656

    def __repr__(self):
        return f"{self.height} x {self.width}: {self.slpsize} vs {self.macssize} | {self.naivesize} vs {self.bvsize} vs {self.panelsize} vs {self.msthrsize} vs {self.mslcesize} vs {self.durbin} vs {self.dynsize}"

def main(argv):
    folder = ""
    try:
        opts, args = getopt.getopt(argv, "hf:", ["folder="])
    except getopt.GetoptError:
        print('plot.py -f <folder>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('plot.py -f <folder>') 
            sys.exit()
        elif opt in ("-f", "--folder"):
            folder = arg
    tmp = folder.split('-')
    height = int(tmp[0].split('_')[1])
    width = int(tmp[1].split('_')[1])

    directory = f"output/rlpbwt/{folder}/input.macs/"
    data_list = []
    for dirname in os.listdir(directory):
        tmp_h = height-int(dirname.split('_')[1])
        tmp_w = width
        tmp_slp = os.path.getsize(f"{directory}{dirname}/panel.slp") * 0.00097656
        tmp_macs = os.path.getsize(f"{directory}{dirname}/panel.macs") * 0.00097656
        tmp_naive = os.path.getsize(f"{directory}{dirname}/naive.ser") * 0.00097656
        tmp_bv = os.path.getsize(f"{directory}{dirname}/bitvector.ser") * 0.00097656
        tmp_panel = os.path.getsize(f"{directory}{dirname}/panel_ext.ser") * 0.00097656
        tmp_thr = os.path.getsize(f"{directory}{dirname}/slp_thr_ext.ser") * 0.00097656
        tmp_lce = os.path.getsize(f"{directory}{dirname}/slp_no_thr_ext.ser") * 0.00097656
        data_list.append(data(tmp_h, tmp_w, tmp_slp, tmp_macs, tmp_naive, tmp_bv, tmp_panel, tmp_thr, tmp_lce + tmp_slp)) 
                
    data_list.sort(key=lambda x: x.height, reverse=False)

    for tmp in data_list:
    	print(f"{tmp.height} & {tmp.width} & {round(tmp.slpsize,2)} & {round(tmp.macssize,2)} & {round(tmp.slpsize/tmp.macssize,6)*100}\\\\")
    	
    print("----------")
    for tmp in data_list:
    	print(f"{tmp.height} & {tmp.width} & {round(tmp.naivesize,2)} & {round(tmp.bvsize,2)} & {round(tmp.panelsize,2)} & {round(tmp.msthrsize,2)} & {round(tmp.mslcesize,2)} & {round(tmp.durbin,2)} & {round(tmp.dynsize,2)}\\\\")
    print("----------")
    for tmp in data_list:
    	print(f"{tmp.height} & {tmp.width} & {round(tmp.mslcesize,2)} & {round(tmp.durbin,2)} & {round(tmp.mslcesize/tmp.durbin,6)*100}\\\\")
    x_slp, y_slp = zip(*[(float(i.height),float(i.slpsize)) for i in data_list])
    plt.scatter(x_slp, y_slp)
    plt.plot(x_slp, y_slp, label="slp")
    x_macs, y_macs = zip(*[(float(i.height),float(i.macssize)) for i in data_list])
    plt.scatter(x_macs, y_macs)
    plt.plot(x_macs, y_macs, label="macs")
    plt.title(f"SLP vs MACs (larghezza = {width})")
    plt.xlabel("Altezza del pannello")
    plt.ylabel("kilobytes")
    plt.yscale("log")
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/slp_vs_macs.png", dpi=900)
    
    plt.clf()
    x_dur, y_dur = zip(*[(float(i.height),float(i.durbin)) for i in data_list])
    plt.scatter(x_dur, y_dur)
    plt.plot(x_dur, y_dur, label="Durbin (~13NM)")
    x_n, y_n = zip(*[(float(i.height),float(i.naivesize)) for i in data_list])
    plt.scatter(x_n, y_n)	
    plt.plot(x_n, y_n, label="naive")
    x_bv, y_bv = zip(*[(float(i.height),float(i.bvsize)) for i in data_list])
    plt.scatter(x_bv, y_bv)
    plt.plot(x_bv, y_bv, label="bitvector")
    x_p, y_p = zip(*[(float(i.height),float(i.panelsize)) for i in data_list])
    plt.scatter(x_p, y_p)
    plt.plot(x_p, y_p, label="dense panel")
    x_thr, y_thr = zip(*[(float(i.height),float(i.msthrsize)) for i in data_list])
    plt.scatter(x_thr, y_thr)
    plt.plot(x_thr, y_thr, label="thr")
    x_lce, y_lce = zip(*[(float(i.height),float(i.mslcesize)) for i in data_list])
    plt.scatter(x_lce, y_lce)
    plt.plot(x_lce, y_lce, label="lce", color="yellow")
    plt.title(f"PBWT vs RLPBWT (larghezza = {width})")
    plt.xlabel("Altezza del pannello")
    plt.ylabel("kilobytes")
    plt.yscale("log")
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/pbwt_vs_rlpbwt.png", dpi=900)

    plt.clf()
    x_dur, y_dur = zip(*[(float(i.height),float(i.durbin)) for i in data_list])
    plt.scatter(x_dur, y_dur)
    plt.plot(x_dur, y_dur, label="Durbin (~13NM)")
    x_dyn, y_dyn = zip(*[(float(i.height),float(i.dynsize)) for i in data_list])
    plt.scatter(x_dyn, y_dyn)
    plt.plot(x_dyn, y_dyn, label="Durbin dynamic O(N+M)")
    x_n, y_n = zip(*[(float(i.height),float(i.naivesize)) for i in data_list])
    plt.scatter(x_n, y_n)
    plt.plot(x_n, y_n, label="naive")
    x_bv, y_bv = zip(*[(float(i.height),float(i.bvsize)) for i in data_list])
    plt.scatter(x_bv, y_bv)    
    plt.plot(x_bv, y_bv, label="bitvector")
    x_p, y_p = zip(*[(float(i.height),float(i.panelsize)) for i in data_list])
    plt.scatter(x_p, y_p)
    plt.plot(x_p, y_p, label="dense panel")
    x_thr, y_thr = zip(*[(float(i.height),float(i.msthrsize)) for i in data_list])
    plt.scatter(x_thr, y_thr)
    plt.plot(x_thr, y_thr, label="thr")
    x_lce, y_lce = zip(*[(float(i.height),float(i.mslcesize)) for i in data_list])
    plt.scatter(x_lce, y_lce)
    plt.plot(x_lce, y_lce, label="lce", color="yellow")
    plt.title(f"PBWT vs RLPBWT (larghezza = {width})")
    plt.xlabel("Altezza del pannello")
    plt.ylabel("kilobytes")
    plt.yscale("log")
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/pbwt_vs_rlpbwt_dyn.png", dpi=900)


    macs=[(35884086533 * 0.00097656, "358653"), (9292947170 * 0.00097656,  "100000"), (4555781714 * 0.00097656,  "46538")]
    slps=[(15125543 * 0.00097656, "358653"), (9295770 * 0.00097656,  "100000"), (8209520 * 0.00097656,  "46538")]
	
    for i in range(len(macs)):
    	print(f"100000 & {slps[i][1]} & {round(slps[i][0],2)} & {round(macs[i][0],2)} & {round(slps[i][0]/macs[i][0],6)*100}\\\\")
    
    
    plt.clf()

    x_slp, y_slp = zip(*[(float(i[1]),float(i[0])) for i in slps])
    plt.scatter(x_slp, y_slp)
    plt.plot(x_slp, y_slp, label="slp")
    x_macs, y_macs = zip(*[(float(i[1]),float(i[0])) for i in macs])
    plt.scatter(x_macs, y_macs)
    plt.plot(x_macs, y_macs, label="macs")
    plt.title(f"SLP vs MACs (altezza = 100000)")
    plt.xlabel("Larghezza del pannello")
    plt.ylabel("kilobytes")
    plt.yscale("log")
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/slp_vs_macs2.png", dpi=900)
    
"""
    runs=[150398, 146071, 141692, 137288,132985,128485, 124136, 119666, 115235, 110648]
    phi=[179396*2,174069*2,168690*2, 163286*2, 157983*2,152483*2, 147134*2, 141664*2,136233*2, 130646*2]
    runs.reverse()
    phi.reverse()
    app=[]
    for i in range(len(runs)):
    	a1 = 2*runs[i]*(2+math.log((data_list[i].height*data_list[i].width)/runs[i]))
    	a2 = 2*(phi[i]/2)*(2+ math.log((data_list[i].height*data_list[i].width)/(phi[i]/2)))
    	#b = 2*4*(data_list[i].width+runs[i]+phi[i])
    	b = 2*(data_list[i].width+runs[i]+(phi[i]/2))
    	c = data_list[i].width
    	d = 128 * (2*data_list[i].height+3*data_list[i].width) 
    	app.append((a1+a2+b+c+d)*0.00097656)
    	#app.append(a+b+c+d)
    
    plt.clf()
    x_lce, y_lce = zip(*[(float(runs[i]),float(data_list[i].mslcesize - data_list[i].slpsize)) for i in range(len(data_list))])
    plt.scatter(x_lce, y_lce)
    plt.plot(x_lce, y_lce, label="lce")
    x_app, y_app = zip(*[(float(runs[i]),float(app[i])) for i in range(len(runs))])
    plt.scatter(x_app, y_app)
    plt.plot(x_app, y_app, label="stima")
    plt.title(f"Stima RLPBWT SLP/LCE")
    plt.xlabel("Numero di run")
    plt.ylabel("Memoria")
    plt.xticks(fontsize=5)
    #plt.yscale("log")
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/runs.png", dpi=900)
	"""
if __name__ == "__main__":
    main(sys.argv[1:])
