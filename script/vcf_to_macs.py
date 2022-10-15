import os
import sys
import getopt


def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    
def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('vcf_to_macs.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('vcf_to_macs.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    with open(inputfile) as f:
        count = 0
        with open(outputfile, 'w') as out:
            out.write(f"NA\n")
            out.write(f"SEED:\t0\n")
            sites = ""
            nsites = 0
            samples = 0
            for line in f:
                print(count, end="\r")
                if line[0] == '#' and line[1] == '#':
                    continue
                elif line[0] == '#' and line[1] != '#':
                    samples = len(line.split()[9:])
                else:
                    data = line.split()
                    nsites += 1
                    sites += f"{data[1]}\t"
                    tmp = "".join(data[9:])
                    tmp = tmp.replace('|', '').strip()
                    out.write(f"SITE:\t{data[1]}\t0\t0\t{tmp}\n")
                count += 1
            out.write(f"TOTAL_SAMPLES:\t {samples*2}\n")
            out.write(f"TOTAL_SITES:\t {nsites}\n")
            out.write(f"BEGIN_SELECTED_SITES\n{sites}\nEND_SELECTED_SITES\n")
        replace_line(outputfile,0, f"COMMAND:\t./macs {samples*2} NA -t NA -r NA\n")
if __name__ == "__main__":
    main(sys.argv[1:])
