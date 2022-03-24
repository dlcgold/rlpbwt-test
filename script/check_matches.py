import sys
import os


class SingleMatchObj:
    def __init__(self, begin, end, length, haplos=None):
        if haplos is None:
            haplos = []
        self.begin = begin
        self.end = end
        self.length = length
        self.haplos = haplos

    def __repr__(self):
        return f"{self.begin}\t{self.end}\t{self.length}\t{self.haplos}\n"

    def __eq__(self, other):
        # check begin, end, length
        if self.begin == other.begin and abs(int(self.end) - int(other.end)) <= 1 and \
                abs(int(self.length) - int(other.length)) <= 1:
            # check for rlpbwtNaive and rlpbwtBitvector
            if self.haplos[0] == "?" or other.haplos[0] == "?":
                if len(self.haplos) == 1 or len(other.haplos) == 1:
                    return True
                elif len(self.haplos) == len(other.haplos):
                    return True
                else:
                    print(f"{len(self.haplos)} vs {len(other.haplos)}")
                    return False
            # check for rlpbwtPanel and rlpbwtSlp not extended
            elif (len(self.haplos) == 1 or len(other.haplos) == 1) and len(self.haplos) != len(other.haplos):
                if len(self.haplos) == 1:
                    if self.haplos[0] in other.haplos:
                        return True
                    else:
                        return False
                else:
                    if other.haplos[0] in self.haplos:
                        return True
                    else:
                        return False
            # check for rlpbwtPanel and rlpbwtSlp extended
            else:
                if self.haplos == other.haplos:
                    return True
                else:
                    print(self.haplos)
                    print(other.haplos)
                    difference = list(set(self.haplos) - set(other.haplos))
                    print(difference)
                    return False
        else:
            print(f"{self.begin} vs {other.begin} \n "
                  f"end: {self.end} vs {other.end}\n"
                  f"length: {self.length} vs {self.length}")
            return False


class QueryMatch:
    def __init__(self, number, matches=None):
        if matches is None:
            matches = []
        self.number = number
        if matches[0].haplos[0] != "?":
            for match in matches:
                match.haplos.sort(key=int)
        self.matches = matches

    def __repr__(self):
        query = f"query {self.number}:\n"
        for matches in self.matches:
            query = query + f"{matches}"
        return query

    def __eq__(self, other):
        if self.number == other.number and self.matches == other.matches:
            return True
        else:
            return False


def match_extractor(file):
    if os.path.isfile(file):
        with open(file) as f:
            prev_query = "0"
            first = True
            query_matches = []
            single_matches = []
            count = 0
            for line in f.readlines():
                c = str(count)
                print(c, end="\r")
                values = line.split("\t")
                query = values[1]
                haplotype = values[2]
                begin = values[3]
                end = values[4]
                length = values[5].strip()
                if query == prev_query and not first:
                    if single_matches[-1].begin == begin and single_matches[-1].end == end:
                        single_matches[-1].haplos.append(haplotype)
                    else:
                        single_matches.append(SingleMatchObj(begin, end, length, [haplotype]))
                    prev_query = query
                else:
                    if first:
                        first = False
                        single_matches.append(SingleMatchObj(begin, end, length, [haplotype]))
                        prev_query = query
                    else:
                        query_matches.append(QueryMatch(prev_query, single_matches))
                        single_matches = [SingleMatchObj(begin, end, length, [haplotype])]
                        prev_query = query
                count = count + 1
    query_matches.append(QueryMatch(prev_query, single_matches))
    return query_matches


def main(argv):
    file1 = os.path.join(argv[1])
    file2 = os.path.join(argv[2])
    # TODO if dynamic not works
    if "Dynamic" in file1 or "Dynamic" in file1:
        print("matches for pbwtDynamic not supported at the moment")
        exit(-1)
    print("extract match first file")
    match1 = match_extractor(file1)
    print("extract match second file")
    match2 = match_extractor(file2)
    correct = True
    print("check matches")
    if len(match1) != len(match2):
        print("DIFFERENT LENGTH")
        correct = False
    else:
        for i in range(0, len(match1)):
            if not (match1[i] == match2[i]):
                print(f"ERROR for match {i}")
                correct = False
            else:
                continue
    if correct:
        print("MATCHES OK")
    else:
        print("ERROR")


if __name__ == "__main__":
    main(sys.argv)
