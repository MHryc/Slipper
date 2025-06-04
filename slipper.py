#!/usr/bin/env python3

import argparse, sys, io, gzip
import pysam

#
# CLI argument parsing
#

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
Path to the input BAM file (needs an index)'''
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
parser.add_argument(
    "-s", "--with-sequence",
    required=False, action="store_true",
    help='''
Store the whole read sequence (SEQ) in the last column. Default is off (fill
with "NA")'''
)

#
# File writing functions
#

def basic_write(write_file: str, seq_switch: bool) -> None:
    """
    Basic function for writing the output tsv file, it's used for: buffered
    plain text tsv writing, gzip compressed tsv writing.

    Args:
        write_file (str): name of the output file
        seq_switch (bool): controls whether SEQ column should contain whole read
            sequence or NA

    Returns:
        None
    """
    # write column headers to output file
    write_file.write(f"QNAME\tFLAGS\tRNAME\tPOS\tCLIP5\tCLIP3\tSEQ\n")

    for read in bamfile.fetch():

        # initialize dictionary with row values, all values must be strings
        r_dict = {
            "QNAME": read.query_name,
            "FLAGS": f"{read.flag}:{read.flag:012b}",
            "RNAME": read.reference_name,
            "POS": str(read.reference_start),
            "CLIP5": "NA",
            "CLIP3": "NA",
            "SEQ": "NA"
        }

        # get CIGAR string in tuple format
        # (4, 7) is the same as 7S
        cigar = read.cigartuples

        # get whole read sequence if arg switch is on
        if seq_switch:
            r_dict["SEQ"] = read.query_sequence

        # get 5' clip sequence
        if cigar[0][0] == 4:
            r_dict["CLIP5"] = read.query_sequence[:cigar[0][1]]

        # get 3' clip sequence
        if cigar[-1][0] == 4:
            r_dict["CLIP3"] = read.query_sequence[-cigar[-1][1]:]

        else:
            pass
        
        # write row into the output tsv file
        write_file.write('\t'.join(r_dict[k] for k in r_dict) + '\n')

def plain_write(out_name: str, seq_switch: bool) -> None:
    """
    Writes the output into a plain text file

    Args:
        out_name (str): name of the output file

    Returns:
        None
    """
    with io.open(out_name, 'w', buffering=524_288) as f:
        basic_write(f, seq_switch)

def gzip_write(out_name: str, seq_switch: bool, gzip_level: int) -> None:
    """
    Writes the output into a gzip compressed plain text file

    Args:
        out_name (str): name of the output file
        c_level (int): gzip compression level, default argparse value sets it to 1

    Returns:
        None
    """
    with gzip.open(out_name, "wt", compresslevel=gzip_level) as f:
        basic_write(f, seq_switch)

#
# Main logic
#

if __name__ == "__main__":
    args = parser.parse_args()

    # open input BAM for reading
    bamfile = pysam.AlignmentFile(args.input, "rb")

    # write processed output to .tsv or .tsv.gz based on the `-g` switch
    if args.gzip:
        out_name = f"{args.output}.tsv.gz"
        gzip_write(
            out_name=out_name,
            seq_switch=args.with_sequence,
            gzip_level=args.compress
        )

    else:
        out_name = f"{args.output}.tsv"
        plain_write(
            out_name=out_name,
            seq_switch=args.with_sequence
        )
    
    # close input BAM
    bamfile.close()