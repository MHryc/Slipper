#!/usr/bin/env python3

from modules.cli import build_parser
from modules.extract import extractor
from modules.analyse import analyser
from modules.plot import ploter


# === Main logic ===

def main() -> None:
    '''
    Main function. Builds arg parser with cli module. Runs a submodule based on
    passed subcommand.
    '''
    # build parser and parse args
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "extract":
        extractor(in_bam=args.input,
                  out_name=args.output,
                  seq_switch=args.with_sequence,
                  gzip_switch=args.gzip,
                  gzip_level=args.compress_level)

    elif args.command == "analyse":
        analyser(args.input,
                 out_name=args.output,
                 briefmode=args.short_output)
    
    elif args.command == "plot":
        ploter(infiles=args.input,
               out_dir=args.outdir)

    return None

if __name__ == "__main__":
    main()
