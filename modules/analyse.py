import io, gzip
from dataclasses import dataclass, astuple

'''
analyse.py, handles the analysis of 'extract' output
'''

'''
What this needs to do:
    * calculate tail length
    * assign strand
    * assign tail type:
        all G or C -> G/C
        Empty (len = 0) -> no_tail
        all T -> polyU
        all A -> polyA
        A or T -> mixed_AU
        otherwise -> other
    * Count terminal Us (Ts)
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

def main_testing(infile: str) -> None:

    with gzip.open(infile, "rt") as f:
        for line in f:
            if line != '':
                print(line)
                print(get_tail_len(line))
                print(count_bases(line))

#        for line in f:
#            print(line.strip().split('\t'))
        #lines = [line.strip().split('\t') for line in f.readlines()]

#    with gzip.open(out_name, "wt", compresslevel=gzip_level) as f:
#        basic_write(in_bam, f, seq_switch)

    return None
