#!/usr/bin/env python3

from modules.cli import build_parser
from modules.extract import extractor
from modules.analyse import main_testing

#
# Main logic
#

if __name__ == "__main__":

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
        main_testing(args.input)
