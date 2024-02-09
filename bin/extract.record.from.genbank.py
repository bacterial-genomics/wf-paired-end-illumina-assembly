#!/usr/bin/env python


import gzip
import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO


def parseArgs():
    parser = ArgumentParser(
        description="extract record(s) from a GenBank file that contain one "
        "or more query and output in FastA format",
        add_help=False,
        epilog="NOTE: GenBank Feature Table Definition is described at "
        "http://www.insdc.org/files/feature_table.html",
    )
    req = parser.add_argument_group("Required")
    req.add_argument(
        "-i",
        "--infile",
        required=True,
        metavar="FILE",
        help="input GenBank file, optionally gunzip compressed",
    )
    req.add_argument(
        "-q",
        "--query",
        nargs="+",
        metavar="STR",
        help="string(s) to search each defline",
    )
    opt = parser.add_argument_group("Optional")
    opt.add_argument(
        "-f",
        "--query-feature",
        default="CDS",
        metavar="STR",
        help="GenBank feature type to search in, e.g., CDS, gene, rRNA, "
        "source, tRNA, misc_feature [default: CDS]",
    )
    opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    opt.add_argument(
        "-o",
        "--outfile",
        default=None,
        metavar="FILE",
        help="FastA output [default: stdout]",
    )
    opt.add_argument(
        "-u",
        "--query-qualifier",
        default="product",
        metavar="STR",
        help="qualifier term within each genbank feature to "
        "search in, e.g., gene, inference, locus_tag, old_locus_tag, "
        "product, translation [default: product]",
    )
    opt.add_argument(
        "--search-type",
        default="any_q_in_rec",
        choices=["all_q_in_rec", "any_q_in_rec", "any_q_is_rec"],
        help="type of query match [default: any_q_in_rec]",
    )
    return parser.parse_args()


def main():
    opt = parseArgs()
    infile = os.path.abspath(os.path.expanduser(opt.infile))
    inbase = os.path.splitext(os.path.basename(infile))[0]
    qry = opt.query
    query_feat, query_qualif = opt.query_feature, opt.query_qualifier

    if infile.endswith(".gz"):
        infile = gzip.open(infile)

    query_match = []
    for rec in SeqIO.parse(infile, "genbank"):
        for feat in rec.features:
            if query_feat == feat.type and query_qualif in feat.qualifiers:
                qual_vals = feat.qualifiers[query_qualif]
                if opt.search_type == "all_q_in_rec":
                    found = [v for v in qual_vals if all(t in v for t in qry)]
                if opt.search_type == "any_q_in_rec":
                    found = [v for v in qual_vals if any(t in v for t in qry)]
                elif opt.search_type == "any_q_is_rec":
                    found = [v for v in qual_vals if any(t == v for t in qry)]
                if len(found) > 0:
                    hit = " ".join(found)
                    locus_tag = feat.qualifiers["locus_tag"][0]
                    query_match.append(
                        ">{} {} {} {}\n{}".format(
                            inbase, locus_tag, rec.name, hit, feat.extract(rec.seq)
                        )
                    )

    if len(query_match) == 0:
        sys.stderr.write("ERROR: {} absent in {}\n".format(qry, infile))
        # sys.exit(1)
    elif len(query_match) > 1:
        sys.stderr.write("WARNING: found >1 {} in {}\n".format(qry, infile))

    if opt.outfile:
        outfile = os.path.abspath(os.path.expanduser(opt.outfile))
        with open(outfile, "w") as ofh:
            for rec in query_match:
                ofh.write("{}\n".format(rec))
    else:
        for rec in query_match:
            print(rec)


if __name__ == "__main__":
    main()
