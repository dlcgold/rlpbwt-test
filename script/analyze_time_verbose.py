import sys
import os
import matplotlib.pyplot as plt

class TimeObj:
    def __init__(self, command=0, version="", user_time=0, sys_time=0, mem=0):
        self.command = command
        self.version = version
        self.user_time = user_time
        self.sys_time = sys_time
        self.mem = mem

    def __repr__(self):
        return f"[{self.command}, {self.version}, " \
               f"{self.user_time}, {self.sys_time}, {self.mem}] "


def main(argv):
    directory = argv[1]
    height = argv[2]
    width = argv[3]
    query_number = argv[4]
    time_obj_list = []
    files = ["pbwtNaive.txt", "pbwtIndexed.txt", "pbwtDynamic.txt",
             "rlpbwtNaive.txt", "rlpbwtBitvector.txt",
             "rlpbwtPanel.txt", "rlpbwtPanelExt.txt",
             "rlpbwtSlpThr.txt", "rlpbwtSlpThrExt.txt",
             "rlpbwtSlpNoThr.txt", "rlpbwtSlpNoThrExt.txt"]
    #for filename in os.listdir(directory):
    for filename in files:
        if "pbwt" in filename or "rlpbwt" in filename:
            file = os.path.join(directory, filename)
            if os.path.isfile(file):
                with open(file) as f:
                    time_obj = TimeObj()
                    lines = f.readlines()
                    if lines[0].strip() == "rlpbwt":
                        time_obj.command = "rlpbwt"
                        if " " in lines[1].strip():
                            tmp = lines[1].strip().replace(" ", "\n")
                            time_obj.version = tmp
                        else:
                            time_obj.version = lines[1].strip()
                        time_obj.user_time = float(lines[2].strip())
                        time_obj.sys_time = float(lines[3].strip())
                        time_obj.mem = float(lines[4].strip())
                    else:
                        time_obj.command = "pbwt"
                        time_obj.version = lines[1].strip()
                        time_obj.user_time = float(lines[2].strip())
                        time_obj.sys_time = float(lines[3].strip())
                        time_obj.mem = float(lines[4].strip())
                    time_obj_list.append(time_obj)
    print(time_obj_list)
    title = f"{query_number} queries on {int(height)-int(query_number)}x{width} panel"
    version = []
    mems = []
    time_user = []
    time_sys = []
    for elem in time_obj_list:
        version.append(elem.version)
        mems.append(elem.mem)
        time_user.append(elem.user_time)
        time_sys.append(elem.sys_time)
    plt.figure()
    plt.yscale('log')
    ax = plt.gca()
    ax.tick_params(axis = 'x', labelsize = 3.5)
    plt.bar(version, mems, color="#5e81ac")
    plt.ylabel("kilobytes")
    plt.xlabel("version")
    plt.title(title)
    plt.savefig("results/mem.png", dpi=200)
    plt.figure()
    plt.yscale('log')
    ax = plt.gca()
    ax.tick_params(axis = 'x', labelsize = 3.5)
    plt.bar(version, time_user, color='#bf616a', label="user time")
    plt.bar(version, time_sys, color='#5e81ac', label="system time") 
    plt.ylabel("seconds")
    plt.xlabel("version")
    plt.title(title)
    plt.legend(fontsize = 7)
    plt.savefig("results/time.png", dpi=200)

if __name__ == "__main__":
    main(sys.argv)
