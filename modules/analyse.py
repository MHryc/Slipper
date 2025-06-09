import io, gzip
import pysam
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
