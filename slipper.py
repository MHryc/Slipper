#!/usr/bin/env python3

import argparse, sys, io, gzip
import pysam

parser = argparse.ArgumentParser(
    prog="Slipper",
    description="Extract soft clipped nucleotides from BAM",
    epilog="M.Hryc (2025)"
)
parser.add_argument(
    "-i", "--input",
    required=True,
    type=str,
    help='''
Path to the input BAM file'''
)
parser.add_argument(
    "-o", "--output",
    required=True,
    type=str,
    help='''
Path to the output file (without file extensions)'''
)
parser.add_argument(
    "-g", "--gzip",
    action="store_true",
    help='''
Output to gzip compressed file'''
)
parser.add_argument(
    "-c", "--compress",
    required=False, default=1, type=int,
    help='''
Gzip compression level (default=1)'''
)

# write functions

def default_write(write_file):
    write_file.write(f"readID\tflags\tclipEnd\tisReverse\tsoftClipLen\tclipSeq\n")

    id = 0

    for read in bamfile.fetch():
        cigar = read.cigartuples
        flags = f"{read.flag}:{read.flag:016b}" # read flags as 0 padded 16bit binary
        is_rev = read.is_reverse

        if cigar[-1][0] == 4:
            strand_end = "clip3"
            sfc_idx = cigar[-1][1] # number of terminal soft clipped bases
            seq = read.query_sequence[-sfc_idx:]

            write_file.write(f"{id}\t{flags}\t{strand_end}\t{is_rev}\t{sfc_idx}\t{seq}\n")

        if cigar[0][0] == 4:
            strand_end = "clip5"
            sfc_idx = cigar[0][1] # number of terminal soft clipped bases
            seq = read.query_sequence[sfc_idx:]

            write_file.write(f"{id}\t{flags}\t{strand_end}\t{is_rev}\t{sfc_idx}\t{seq}\n")

        else:
            pass

        id += 1

def gzip_write(out_name: str, c_level: int):
    with gzip.open(out_name, "wt", compresslevel=c_level) as f:
        default_write(f)

def plain_write(out_name: str):
    with io.open(out_name, 'w', buffering=524_288) as f:
        default_write(f)

if __name__ == "__main__":
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.input, "rb")

    if args.gzip:
        out_name = f"{args.output}.tsv.gz"
        gzip_write(out_name, args.compress)

    else:
        out_name = f"{args.output}.tsv"
        plain_write(out_name)
            
    bamfile.close()