import os
import sys

def main(argv):
    print(argv)
    with open(argv[0]) as f:
        count = 0
        zeros = 0
        ones = 0
        print("read file")
        for line in f:
            c = str(count)
            
            if line.strip().split()[0] == 'TOTAL_SAMPLES:':
                break
            if count > 1:
                # out.write(line.strip().split("\t")[4])
                tmp = line.strip().split()[4]
                zeros += tmp.count('0')
                ones += tmp.count('1')
            count += 1
        print(zeros, ones)
if __name__ == "__main__":
    main(sys.argv[1:])
