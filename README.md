# Slipper

Soft-clipper, script for extracting soft clipped bases from BAM input.

## Dependencies

* pysam

# Installation

Run:

```
git clone https://github.com/MHryc/Slipper.git
cd Slipper/
chmod +x Slipper slipper.py # just in case, should not be necessary
```

`Slipper` is just a symlink to `slipper.py`. You might want to create another
one in your `$PATH`, eg. `~/.local/bin` for easier access.

# Usage

```
usage: Slipper [-h] -i INPUT -o OUTPUT [-g] [-c LEVEL] [-s]

Extract soft clipped bases from BAM. For output description see README.md

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input BAM file (needs an index)
  -o OUTPUT, --output OUTPUT
                        Path to the output file (without file extensions)
  -g, --gzip            Output to gzip compressed file
  -c LEVEL, --compress LEVEL
                        Gzip compression level (default=2)
  -s, --with-sequence   Store the whole read sequence (SEQ) in the last column. Default is off (fill with "NA")

M.Hryc (2025)
```

# Output format

Output tsv has 7 columns:

| name  | description                                                                   |
| ---   | ---                                                                           |
| QNAME | Query template NAME                                                           |
| FLAGS | bitwise FLAG both as int and 12-bit                                           |
| RNAME | Reference sequence NAME                                                       |
| POS   | 0-based leftmost tail POSition, i.e. the 1st base AFTER the last aligned base |
| CLIP5 | 5' soft clipped SEQuence (NA if not present)                                  |
| CLIP3 | 3' soft clipped SEQuence (NA if not present)                                  |
| SEQ   | segment SEQuence (read sequence)                                              |

## Example output

```
QNAME   FLAGS   RNAME   POS     CLIP5   CLIP3   SEQ
A00805:128:H5W25DRX2:1:2242:1190:26882  163:000010100011        1       15161   TTAATAGGAACTCATGTC      NA      TTAATAGGAACTCATGTCAGCCTCAGTCGATCCCTGACCCAGCACCGGGCACTGATGAGACAGCGGCGGTTTGAGGAGCCACCTCCCAGCCACCTCGGGGC
A00805:128:H5W25DRX2:2:2263:27462:7842  419:000110100011        1       15274   TGCATAATGAAAAGCGTCAGG   NA      TGCATAATGAAAAGCGTCAGGATGGGGAAACTGGCCCAGAGAGGTGAGGCAGCTTGCCTGGGGTCACAGAGCAAGGCAAAAGCAGCGCTGGGTACAAGCTC
A00805:128:H5W25DRX2:1:2242:1190:26882  83:000001010011 1       15281   NA      NA      AACTGGCCCAGAGAGGTGAGGCAGCTTGCCTGGGGTCACAGAGCAAGGCAAAAGCAGCGCTGGGTACAAGCTCAAAACCATAGTGCCCAGGGCACTGCCGC
A00805:128:H5W25DRX2:2:2263:27462:7842  339:000101010011        1       15338   NA      NA      CGCTGGGTACAAGCTCAAAACCATAGTGCCCAGGGCACTGCCGCTGCAGGCGCAGGCATCGCATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGT
A00805:128:H5W25DRX2:2:1169:26775:16924 419:000110100011        1       17055   GAGTACGATGTGGTAGTCAGT   NA      GAGTACGATGTGGTAGTCAGTCTGCAAGATTAGGCAGGGACATGTGAGAGGTGACAGGGACCTGCAGGGGCAGCCAACAAGACCTTGTGTGCACCTCCCAT
A00805:128:H5W25DRX2:2:1169:26775:16924 339:000101010011        1       17131   NA      NA      CCATGGGTGGAATAAGGGGCCCAACAGCCTTGACTGGAGAGGAGCTCTGGCAAGGCCCTGGGCCACTGCACCTGTCTCCACCTCTGTCCCACCCCTCCCAC
A00805:128:H5W25DRX2:1:1106:17463:27101 419:000110100011        1       17611   CTTGATACTATATAAGTCAGTGGAAAAAA   NA      CTTGATACTATATAAGTCAGTGGAAAAAATCGGTGGTGTTGAAGAGCAGCAAGGAGCTGACAGAGCTGATGTTGCTGGGAAGACCCCCAAGTCCCTCTTCT
A00805:128:H5W25DRX2:2:2263:18611:10207 419:000110100011        1       22662   GTGACAATGGTAAGGGTCAG    NA      GTGACAATGGTAAGGGTCAGGGCCCTCCCTGGGCTGTGCCAGCAGCTTGGAGAACCCACACTCAATGAACGCAGCACTCCACTACCCAGGAAATGCCTTCC
A00805:128:H5W25DRX2:2:2263:18069:16939 419:000110100011        1       22662   GTGACAATGGTAAGGGTCAG    NA      GTGACAATGGTAAGGGTCAGGGCCCTCCCTGGGCTGTGCCAGCAGCTTGGAGAACCCACACTCAATGAACGCAGCACTCCACTACCCAGGAAATGCCTTCC
```
