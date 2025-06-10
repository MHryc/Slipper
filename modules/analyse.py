import io, gzip
from dataclasses import dataclass, astuple

'''
analyse.py, handles the analysis of 'extract' output
'''

'''
What this needs to do:
    * calculate tail length             DONE
    * count nucleotides                 DONE
    * assign strand                     DONE
    * assign tail type:                 DONE
        all G or C -> G/C               DONE
        Empty (len = 0) -> no_tail      DONE
        all T -> polyU                  DONE
        all A -> polyA                  DONE
        A or T -> mixed_AU              DONE
        otherwise -> other              DONE
    * Count terminal Us (Ts)            DONE
'''

@dataclass
class TsvLine:
    '''
    Each field needs to be a str, because they're '\t'.join't before writing to file.
    '''
    QNAME: str
    FLAGS: str
    RNAME: str
    POS: str
    CLIP5: str
    CLIP3: str
    SEQ: str

def get_tail_len(line: list[str, ...]) -> int:
    if line[5] != "NA":
        return len(line[5])
    else:
        return 0

def count_bases(line: list[str, ...]) -> tuple[int, int, int, int]:
    tail_seq = line[5]
    if tail_seq != "NA":
        return tuple(
                [tail_seq.count('A'),
                 tail_seq.count('C'),
                 tail_seq.count('G'),
                 tail_seq.count('T')])
    else:
        return tuple([0, 0, 0, 0])

def get_strand(line: list[str, ...]) -> int:
    flags = bin(int(line[1].split(':')[0]))
    return -1 if flags[-5] else 1

def get_tail_type(base_counts: tuple[int, int, int, int]) -> str:
    '''
    Args:
        base_counts (tuple): (A, C, G, T)
    '''
    total = sum(base_counts)
    if total == 0:
        return "no_tail"

    elif base_counts[3] == total:
        return "polyU"

    elif base_counts[0] == total:
        return "polyA"

    elif base_counts[0] + base_counts[3] == total:
        return "mixed_AU"

    elif base_counts[2] + base_counts[3] == total:
        return "mixed_GC"

    else:
        return "other"

def count_Us(line: list[str]) -> int:
    tail = line[5]
    u_count = 0
    for i in range(1, len(tail) + 1):
        if tail[-i] == 'T':
            u_count += 1
        else:
            break

    return u_count

def main_testing(infile: str) -> None:

    with gzip.open(infile, "rt") as fi, gzip.open("analyse_out.tsv.gz", "wt") as fo:
        header = fi.readline() # throw away the header
        for line in fi:
            if line != '':
                l = line.strip().split('\t')
                print(l)
                print(get_tail_len(l))
                counts = count_bases(l)
                print(get_strand(l))
                print(get_tail_type(counts))
                print(count_Us(l))

                fo.write(f"{line.strip()}\t{get_tail_type(counts)}\n")

#        for line in f:
#            print(line.strip().split('\t'))
        #lines = [line.strip().split('\t') for line in f.readlines()]

#    with gzip.open(out_name, "wt", compresslevel=gzip_level) as f:
#        basic_write(in_bam, f, seq_switch)

    return None
