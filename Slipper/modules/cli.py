import argparse

'''
cli.py, handles argument parsing
'''

def build_parser() -> argparse.ArgumentParser:
    '''
    Build CLI arg parser for the main script

    Args:
        None

    Outputs:
        argparse.ArgumentParser: parser object
    '''
    parser = argparse.ArgumentParser(
        prog="Slipper",
        description='''
        Extract soft clipped bases from BAM. For output description see repo's
        README.md
        ''', epilog="M.Hryc (2025)"
    )

    subparsers = parser.add_subparsers(dest='command')

    # === the 'extract' subcommand ===

    extract = subparsers.add_parser(
        "extract", help='''
        Extract soft clipped bases from BAM
        '''
    )
    extract.add_argument(
        "-i", "--input",
        required=True,
        type=str, help='''
    Path to the input BAM file (needs an index)'''
    )
    extract.add_argument(
        "-o", "--output",
        required=True,
        type=str, help='''
    Path to the output file (without file extensions)'''
    )
    extract.add_argument(
        "-g", "--gzip",
        action="store_true",
        help='''
    Output to gzip compressed file'''
    )
    extract.add_argument(
        "-c", "--compress-level", metavar="LEVEL",
        required=False, default=1, type=int,
        help='''
    Gzip compression level (default=2)'''
    )
    extract.add_argument(
        "-s", "--with-sequence",
        required=False, action="store_true",
        help='''
    Store the whole read sequence (SEQ) in the last column. Default is off
    (fill with "NA")'''
    )

    # === the 'analyse' subcommand

    analyse = subparsers.add_parser(
        "analyse", help='''
        Analyse the output of 'Slipper extract'
        '''
    )
    analyse.add_argument(
        "-i", "--input",
        required=True,
        type=str, help='''
        Path to the input TSV produced by 'Slipper extract'
        '''
    )
    analyse.add_argument(
        "--short-output", action="store_true",
        help='''
        Skip the lines with 'no_tail' in the output
        '''
    )
    analyse.add_argument(
        "-o", "--output",
        required=True,
        type=str, help='''
        Path to the output file (without file extensions)
        '''
    )

    # === the 'plot' subcommand ===

    plot = subparsers.add_parser(
        "plot", help='''
        Make plots from 'Slipper analyse' output
        '''
    )
    plot.add_argument(
        "-i", "--input",
        required=True, type=str, 
        nargs='+', help='''
        Path to the input TSVs produced by 'Slipper analyse'
        '''
    )
    plot.add_argument(
        "-o", "--outdir",
        default="slipper_plots",
        type=str, help='''
        Path to directory, where plots will be saved. Defaults to slipper_plots/
        '''
    )

    return parser
