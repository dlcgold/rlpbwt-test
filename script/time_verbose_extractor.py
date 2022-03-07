import sys
import getopt


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('time_verbose_extractor.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('time_verbose_extractor.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    if inputfile == "done.txt":
        exit
    with open(inputfile, "r") as f:
        with open(outputfile, "w") as out:
            for line in f:
                line = line[1:-1]
                tokens = line.split(sep=":")
                if tokens[0] == "Command being timed":
                    # out.write(tokens[1].lstrip())
                    # out.write("\n")
                    if '/rlpbwt' in tokens[1]:
                        out.write("rlpbwt\n")
                        if '-N' in tokens[1]:
                            out.write("naive\n")
                        if '-B' in tokens[1]:
                            out.write("bitvectors\n")
                        if '-P' in tokens[1]:
                            if '-e' in tokens[1]:
                                out.write("panel extended\n")
                            else:
                                out.write("panel\n")
                        if '-S' in tokens[1]:
                            if '-e' in tokens[1] and '-t' in tokens[1]:
                                out.write("slp_thr extended\n")
                            elif '-e' in tokens[1] and '-t' not in tokens[1]:
                                out.write("slp_no_thr extended \n")
                            elif '-e' not in tokens[1] and '-t' not in tokens[1]:
                                out.write("slp_no_thr \n")
                            else:
                                out.write("slp_thr \n")
                                
                    if '/pbwt' in tokens[1]:
                        out.write("pbwt\n")
                        if '-matchNaive' in tokens[1]:
                             out.write("original\n")
                        if '-matchIndexed' in tokens[1]:
                             out.write("indexed\n")
                        if '-matchDynamic' in tokens[1]:
                             out.write("dynamic\n")
                        
                if tokens[0] == "User time (seconds)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")
                if tokens[0] == "System time (seconds)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")
                if tokens[0] == "Maximum resident set size (kbytes)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")


if __name__ == "__main__":
    main(sys.argv[1:])
