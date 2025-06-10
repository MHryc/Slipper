import io, gzip
import pysam
from dataclasses import dataclass, astuple

'''
extract.py, handles logic of the 'extract' subcommand
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

def basic_write(
        in_bam: pysam.AlignmentFile,
        write_file: str,
        seq_switch: bool) -> None:
    """
    Basic function for writing the output tsv file, it's used for: buffered
    plain text tsv writing and gzip compressed tsv writing within plain_write()
    and gzip_write().

    Args:
        in_bam (str): path to input BAM file
        write_file (str): name of the output file
        seq_switch (bool): controls whether SEQ column should contain whole read
            sequence or NA

    Returns:
        None
    """

    # write column headers to output file
    write_file.write(f"QNAME\tFLAGS\tRNAME\tPOS\tCLIP5\tCLIP3\tSEQ\n")

    # iterate over bam lines
    for read in in_bam.fetch():

        # get CIGAR string in tuple format
        # eg. (4, 7) is the same as 7S
        cigar = read.cigartuples

        # initialize dataclass with row values, all values must be strings
        line = TsvLine(
            read.query_name,
            f"{read.flag}:{read.flag:012b}",
            read.reference_name,
            str(read.reference_start),
            read.query_sequence[:cigar[0][1]] if cigar[0][0] == 4 else "NA",
            read.query_sequence[:cigar[-1][1]] if cigar[-1][0] == 4 else "NA",
            read.query_sequence if seq_switch else "NA"
        )

        # write row into the output tsv file
        write_file.write('\t'.join(astuple(line)) + '\n')

    return None

def plain_write(
        in_bam: str,
        out_name: str,
        seq_switch: bool=False
        ) -> None:
    """
    Writes the output into a plain text file

    Args:
        in_bam (str): path to input BAM file
        out_name (str): name of the output file
        seq_switch (bool): controls whether the original read SEQ should be
            included

    Returns:
        None
    """
    with io.open(out_name, 'w', buffering=524_288) as f:
        basic_write(in_bam, f, seq_switch)

    return None

def gzip_write(
        in_bam: str,
        out_name: str,
        gzip_level: int,
        seq_switch: bool=False) -> None:
    """
    Writes the output into a gzip compressed plain text file

    Args:
        in_bam (str): path to input BAM file
        out_name (str): name of the output file
        c_level (int): gzip compression level, default argparse value sets it to 1
        seq_switch (bool): controls whether the original read SEQ should be
            included

    Returns:
        None
    """
    with gzip.open(out_name, "wt", compresslevel=gzip_level) as f:
        basic_write(in_bam, f, seq_switch)

    return None

def extractor(in_bam: str,
              out_name: str,
              seq_switch: bool,
              gzip_switch: bool,
              gzip_level: int) -> None:
    '''
    Function for handling the main logic of 'extractor' subcommand

    Args:
        in_bam (str): input BAM path
        out_name (str): name of output TSV
        seq_switch (bool): controls whether to write the original SEQ to the
            TSV
        gzip_switch (bool): controls whether to gzip the output TSV
        gzip_level (int): gzip compression level, default is 2 (set in cli.py
            parsing module)

    Returns:
        None
    '''
    # open input BAM for reading
    bamfile = pysam.AlignmentFile(in_bam, "rb")

    # write to tsv.gz or tsv
    if gzip_switch:
        gzip_write(
            in_bam=bamfile,
            out_name=f"{out_name}.tsv.gz",
            seq_switch=seq_switch,
            gzip_level=gzip_level
        )

    else:
        plain_write(
            in_bam=bamfile,
            out_name=f"{out_name}.tsv",
            seq_switch=seq_switch
        )

        # close input BAM
        bamfile.close()

    return None
