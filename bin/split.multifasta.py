#!/usr/bin/env python


import os
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parseArgs():
    parser = ArgumentParser(
        add_help=False,
        description="splits multi-FastA file into individual FastA files",
        epilog="NOTE: for each split file, the output defline will be the "
        "filename excluding the extension",
    )
    req = parser.add_argument_group("Required")
    req.add_argument(
        "-i", "--infile", required=True, metavar="FILE", help="input multi-FastA file"
    )
    opt = parser.add_argument_group("Optional")
    opt.add_argument(
        "-e",
        "--ext",
        metavar="STR",
        type=str,
        default=".fasta",
        help="extension name for each output file ['.fasta']",
    )
    opt.add_argument(
        "-g",
        "--nogaps",
        action="store_true",
        default=False,
        help="toggle on removal of '-' gaps [off]",
    )
    opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    opt.add_argument(
        "-o",
        "--outdir",
        metavar="PATH",
        help="output directory path [dir of input file]",
    )
    opt.add_argument(
        "-p",
        "--prefix",
        metavar="STR",
        type=str,
        default=None,
        help="prefix name for each output file [defline from input file; "
        "first string on whitespace]",
    )
    opt.add_argument(
        "-s",
        "--suffix",
        metavar="STR",
        type=str,
        default=None,
        help="suffix name for each output file ['_<int++>']",
    )
    return parser.parse_args()


def main():
    opt = parseArgs()

    # I/O handling
    infile = os.path.abspath(os.path.expanduser(opt.infile))
    if opt.outdir is None:
        outdir = os.path.dirname(infile)
    else:
        outdir = os.path.abspath(os.path.expanduser(opt.outdir))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    mfasta = SeqIO.parse(infile, "fasta")

    i = 1
    for rec in mfasta:
        # Sequence handling
        seq = str(rec.seq).upper()
        if opt.nogaps:
            seq = seq.replace("-", "")

        # Defline and filename handling
        if opt.prefix is not None:
            pref = opt.prefix
        else:
            pref = (rec.id).split()[0]
        if opt.suffix is not None:
            suff = opt.suffix
        else:
            suff = "_{}".format(i)
            i += 1

        # Output each record as separate file
        new_rec = SeqRecord(Seq(seq), id=pref + suff, description="")
        SeqIO.write(new_rec, os.path.join(outdir, pref + suff + opt.ext), "fasta")


if __name__ == "__main__":
    main()
