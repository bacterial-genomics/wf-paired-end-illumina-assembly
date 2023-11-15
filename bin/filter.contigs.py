#!/usr/bin/env python


import argparse
import gzip
import os
import re
import sys

from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC


def parseArgs():
    parser = argparse.ArgumentParser(
        description="filters contigs (or "
        "scaffolds) based on length, coverage, GC skew, and compositional "
        "complexity",
        add_help=False,
    )
    req = parser.add_argument_group("Required")
    req.add_argument(
        "-i",
        "--infile",
        metavar="FILE",
        required=True,
        help="input FastA file from IDBA, SKESA, SPAdes, "
        "Unicycler, or Velvet, optionally gunzip compressed",
    )
    opt = parser.add_argument_group("Optional")
    opt.add_argument(
        "-b",
        "--baseheader",
        default=None,
        metavar="STR",
        help="contig header prefix (with _<CNTR> suffix) " "[basename infile .ext]",
    )
    opt.add_argument(
        "-c",
        "--cov",
        type=str,
        default="5",
        metavar="INT|FLOAT",
        help="minimum coverage; set to 0 to skip; an integer of at least 1 "
        "applies a raw absolute coverage value filter whereas a float value "
        "filters sequences that meet a minimum percentage relative to the "
        "length-normalized median coverage of sequences that first pass all "
        "other filters [5]",
    )
    opt.add_argument(
        "-d",
        "--discarded",
        metavar="FILE",
        default=None,
        help="output FastA file of discarded sequences which includes "
        "failed filters in the deflines [None]",
    )
    opt.add_argument(
        "-g",
        "--gcskew",
        default=True,
        action="store_false",
        help="keep >88 and <12%% GC contigs",
    )
    opt.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    opt.add_argument(
        "-l",
        "--len",
        type=int,
        default=500,
        metavar="INT",
        help="minimum contig length (in bp) [500]",
    )
    opt.add_argument(
        "-m",
        "--complex",
        default=True,
        action="store_false",
        help="keep low-complexity contigs with only 1 or 2 nucleotides " "represented",
    )
    opt.add_argument(
        "-o",
        "--outfile",
        default=None,
        metavar="FILE",
        help="output filtered FastA file [stdout]",
    )
    opt.add_argument(
        "--deflines",
        choices={"rename", "rename_retain", "retain"},
        default="rename_retain",
        help="rename: fully rename "
        "with -b prefix and _<CNTR> suffix; retain: keep exactly as input; "
        "rename_retain: includes rename style prefix but also adds "
        "OrigDefln=<input defline> suffix [rename_retain]",
    )
    opt.add_argument(
        "--no-sort",
        default=False,
        action="store_true",
        help="skip sorting contigs by descending length",
    )
    opt.add_argument(
        "--quiet-removal",
        default=False,
        action="store_true",
        help="skip explanations for each record removal",
    )
    opt.add_argument(
        "--quiet-coverage",
        default=False,
        action="store_true",
        help="skip reporting coverage statistics (min, Q25, mean, median, " "Q75, max)",
    )
    opt.add_argument(
        "--quiet-stats",
        default=False,
        action="store_true",
        help="skip final reporting of input, discard, and output tallies "
        "(records count, sequence lengths)",
    )
    opt.add_argument(
        "--silent", default=False, action="store_true", help="no screen output"
    )
    return parser.parse_args()


# SPAdes v3.12.0 and Velvet v1.2.10 are the identical
# SKESA v2.2 appends _Circ after cov if suggested circular contig
# IDBA (requires replacing spaces into underscores
#       due to Biopython SeqIO's parsing)
# Unicycler v0.4.7 appends x to relative cov val
ASSEM_REGEXS = {
    "idba": "read_count_([0-9]+)",
    "skesa": "Contig_[0-9]{1,}_[0-9]{1,}.[0-9]{1,}",
    "unicyc": "depth=[0-9]{1,}(\.[0-9]{1,})?",
    "velvet": "cov_([0-9]{1,}.[0-9]{1,})",
}


def string_is_int_or_float(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def string_is_float(s):
    try:
        return float(s) and "." in s
    except ValueError:
        return False


def require_minimum_arg_value(arg_value, arg_name, min_value):
    if arg_value < min_value:
        sys.stderr.write(
            "ERROR: minimum {} must be at least {}\n".format(arg_name, arg_value)
        )
        sys.exit(1)


def detect_assembly_format(record, quiet):
    """
    given a Bio.SeqRecord.SeqRecord object, returns the assembler software
    name used to generate the input FastA file based on the coverage defline
    """
    for k, v in ASSEM_REGEXS.items():
        cov_regex = re.compile(v)
        cov_match = cov_regex.search(record.description)
        if cov_match:
            if not quiet:
                sys.stderr.write("INFO: detected {} defline format\n\n".format(k))
            return k, cov_regex
    else:
        sys.stderr.write(
            "ERROR: unable to parse coverage from first defline\n"
            "Expect raw output file from SKESA, SPAdes, or Velvet assemblers "
            "or an IDBA file with undescores instead of spaces in the deflines. "
            "Coverage cannot be filtered for other formats, so to apply "
            "other filters and skip coverage filtering, set coverage to 0.\n"
        )
        sys.exit(1)


def extract_coverage_value_from_defline(record, assembly_fmt, cov_regex):
    """
    given a Bio.SeqRecord.SeqRecord object and the assembler name,
    extracts and returns the coverage value found in its defline
    """
    cov_val = 0
    cov_match = cov_regex.search(record.description)
    if cov_match:
        if assembly_fmt == "idba":
            cov_val = int((cov_match.group(0)).lstrip("read_count_"))
        elif assembly_fmt == "skesa":
            cov_val = float((cov_match.group(0)).split("_")[-1])
        elif assembly_fmt == "unicyc":
            cov_val = float((cov_match.group(0)).lstrip("depth="))
        elif assembly_fmt == "velvet":
            cov_val = float((cov_match.group(0)).lstrip("cov_"))
    return cov_val


def filter_record_length(record, min_len, quiet):
    """
    evaluates a Bio.SeqRecord.SeqRecord object for sequence length
    """
    seq_len = len(record.seq)
    if seq_len >= min_len:
        return seq_len
    else:
        if not quiet:
            sys.stderr.write("INFO: {} too short\n".format(record.id))
        return False


def filter_record_coverage(record, min_cov, assembly_fmt, cov_regex, quiet):
    """
    evaluates a Bio.SeqRecord.SeqRecord object for read coverage depth
    """
    cov = extract_coverage_value_from_defline(record, assembly_fmt, cov_regex)
    if cov >= min_cov:
        return "{:.2f}".format(cov)
    else:
        if not quiet:
            sys.stderr.write("INFO: {} low coverage\n".format(record.id))
        return False


def filter_record_gc(record, bool_gc, quiet):
    """
    evaluates a Bio.SeqRecord.SeqRecord object for GC composition; note
    that biopython represents GC content as a percentage of the total
    sequence length and not just of ATCG
    """
    if bool_gc:
        gc_content = float(GC(record.seq))
        if 12 <= gc_content <= 88:
            return "{:.2f}".format(gc_content)
        else:
            if not quiet:
                sys.stderr.write("INFO: {} skewed GC\n".format(record.id))
            return False
    else:
        return "skipped filter"


def filter_record_complexity(record, bool_complex, quiet):
    """
    evaluates a Bio.SeqRecord.SeqRecord object for at least 3 unique
    unambiguous nucleotides present
    """
    if bool_complex:
        i = 0
        for nucleotide in ("A", "T", "C", "G"):
            if nucleotide in record.seq:
                i += 1
        if i > 2:
            return i
        else:
            if not quiet:
                sys.stderr.write(
                    "INFO: {} low compositional complexity\n".format(record.id)
                )
            return False
    else:
        return "skipped filter"


def calc_median(srt_floats_list):
    """
    takes a sorted list of floats and returns a tuple of median value and a
    list of the index(es) of the median value(s)
    """
    idxs = []
    len_list = len(srt_floats_list)
    if len_list % 2 == 0:
        idxs.append(int(len_list / 2) - 1)
        idxs.append(int(len_list / 2))
        median = (srt_floats_list[idxs[0]] + srt_floats_list[idxs[1]]) / 2.0
    else:
        idxs.append(int(len_list / 2))
        median = srt_floats_list[idxs[0]]
    return median, idxs


def calc_normalized_median_cov(
    records, assembly_fmt, cov_regex, coverage_filter_type, quiet
):
    """
    calculates the length-normalized median coverage value
    for a list of SeqIO records
    """
    seqcoverages = []
    for rec in records:
        seq_len = len(rec.seq)
        cov = extract_coverage_value_from_defline(rec, assembly_fmt, cov_regex)
        seqcoverages.extend(seq_len * [cov])
    srt_covs = sorted(seqcoverages)
    Q50, median_idxs = calc_median(srt_covs)
    if not quiet:
        d = {"abs": "in", "rel": "ex"}
        Q25, _ = calc_median(srt_covs[: median_idxs[0]])
        Q75, _ = calc_median(srt_covs[median_idxs[-1] + 1 :])
        sys.stderr.write(
            "\nINFO: coverage quartiles of filtered records\n"
            "    ({}cluding cov filt)\n"
            "    25th={:.2f}\n"
            "    50th={:.2f}\n"
            "    75th={:.2f}\n".format(d[coverage_filter_type], Q25, Q50, Q75)
        )
    return Q50


def rename_record_defline(record, baseheader, i, rename_type):
    """
    renames sequence record object according to specified type
    """
    if rename_type == "rename_retain":
        return SeqRecord.SeqRecord(
            id="{}_{} OrigDefln={}".format(baseheader, i, record.description),
            seq=record.seq,
            description="",
        )
    elif rename_type == "rename":
        return SeqRecord.SeqRecord(
            id="{}_{}".format(baseheader, i), seq=record.seq, description=""
        )
    elif rename_type == "retain":
        return SeqRecord.SeqRecord(
            id=record.description, seq=record.seq, description=""
        )
    elif rename_type == "discarded":
        return SeqRecord.SeqRecord(
            id="{}_{} {} OrigDefln={}".format(
                baseheader, i, record.description, record.id
            ),
            seq=record.seq,
            description="",
        )


def write_records(records, baseheader, defline_rename, outfile, bool_nosortbylen):
    if not bool_nosortbylen:
        records = sorted(records, key=lambda x: x[1], reverse=True)
    out_recs = []
    i = 1
    for r, _ in records:
        out_recs.append(rename_record_defline(r, baseheader, i, defline_rename))
        i += 1
    SeqIO.write(out_recs, outfile, "fasta")


def main():
    args = parseArgs()
    if not string_is_int_or_float(args.cov):
        sys.stderr.write(
            "ERROR: --cov {} must be a float or integer\n".format(args.cov)
        )
        sys.exit(1)
    min_cov, min_len = float(args.cov), args.len
    require_minimum_arg_value(min_cov, "coverage", 0)
    require_minimum_arg_value(min_len, "length", 1)
    infile, outfile = args.infile, args.outfile
    defline_rename = args.deflines
    bool_gc, bool_complex = (
        args.gcskew,
        args.complex,
    )
    bool_nosortbylen = args.no_sort
    silent, quiet_coverage = args.silent, args.quiet_coverage
    quiet_removal, quiet_stats = args.quiet_removal, args.quiet_stats
    if silent:
        quiet_removal = quiet_coverage = quiet_stats = True

    # Output handling
    if args.baseheader is None:
        baseheader = os.path.splitext(os.path.basename(infile.strip(".gz")))[0]
    else:
        baseheader = args.baseheader
    if outfile is None:
        out = sys.stdout
    else:
        out = os.path.abspath(os.path.expanduser(outfile))

    # Auto-detections
    if string_is_float(args.cov):
        coverage_filter_type = "rel"
        min_cov = 0
    else:
        coverage_filter_type = "abs"
    if infile.endswith(".gz"):
        infile = gzip.open(infile)

    # Enable skiping coverage filter
    if float(args.cov) > 0:
        first_record = next(SeqIO.parse(infile, "fasta"))
        assembly_fmt, cov_regex = detect_assembly_format(first_record, silent)
    else:
        assembly_fmt, cov_regex = "post-processed", re.compile("xX0Xx")
        quiet_coverage = True

    # Filter each sequence record
    records, discarded = [], []  # items are tuples of seqrec obj and seqlen int
    cnts_input_reclengths = []
    for record in SeqIO.parse(infile, "fasta"):
        cnts_input_reclengths.append(len(record.seq))
        d = {}
        d["length"] = filter_record_length(record, min_len, quiet_removal)
        d["coverage"] = filter_record_coverage(
            record, min_cov, assembly_fmt, cov_regex, quiet_removal
        )
        d["gc_content"] = filter_record_gc(record, bool_gc, quiet_removal)
        d["complexity"] = filter_record_complexity(record, bool_complex, quiet_removal)
        if all(d.values()):
            records.append((record, d["length"]))
        else:
            record.description = ""
            for k, v in d.items():
                if not v:
                    record.description += "Failed={}".format(k)
            discarded.append((record, len(record.seq)))

    # Calculate statistics for read depth of contig coverage
    if len(records) == 0:
        sys.stderr.write(
            "\n\nERROR: no records passed filters\n"
            "To diagnose, check which filter(s) responsible for this in "
            "the stderr info. Coverage is most suspect, especially if input "
            "has relative depths such as Unicycler, which can be quickly "
            "confirmed if the number of 'low coverage' in stderr equals "
            "the number of input records.\n"
        )
        sys.exit(1)
    unfilt_records = [r[0] for r in records]
    median_cov = calc_normalized_median_cov(
        unfilt_records, assembly_fmt, cov_regex, coverage_filter_type, quiet_coverage
    )
    unfilt_coverages = [
        extract_coverage_value_from_defline(r, assembly_fmt, cov_regex)
        for r in unfilt_records
    ]
    if not quiet_coverage:
        min_val, max_val = min(unfilt_coverages), max(unfilt_coverages)
        avg_cov = sum(unfilt_coverages) / float(len(unfilt_coverages))
        d = {"abs": "in", "rel": "ex"}
        sys.stderr.write(
            "\nINFO: coverage of filtered records\n"
            "    ({}cluding cov filt)\n"
            "    minimum={:.2f}\n"
            "    average={:.2f}\n"
            "    maximum={:.2f}\n\n".format(
                d[coverage_filter_type], min_val, avg_cov, max_val
            )
        )

    # Relative coverage filtering of filtered recs only
    if coverage_filter_type == "rel" and len(records) > 0:
        min_cov = median_cov * float(args.cov)
        records = []
        for rec in unfilt_records:
            cov_val = filter_record_coverage(
                rec, min_cov, assembly_fmt, cov_regex, quiet_removal
            )
            if cov_val:
                records.append((rec, len(rec.seq)))
            else:
                rec.description = "Failed=coverage"
                discarded.append((rec, len(rec.seq)))
        if not silent:
            sys.stderr.write(
                "\nINFO: {} minimum coverage applied from "
                "{} input relative coverage value * "
                "{:.2f} calculated median coverage\n\n".format(
                    min_cov, args.cov, median_cov
                )
            )

    # Optionally report tallies
    if not quiet_stats:
        sys.stderr.write("INFO: {} records input\n".format(len(cnts_input_reclengths)))
        sys.stderr.write("INFO: {} bp input\n".format(sum(cnts_input_reclengths)))
        sys.stderr.write("INFO: {} records discarded\n".format(len(discarded)))
        sys.stderr.write("INFO: {} bp discarded\n".format(sum(n for _, n in discarded)))
        sys.stderr.write("INFO: {} records output\n".format(len(records)))
        sys.stderr.write("INFO: {} bp output\n\n".format(sum(n for _, n in records)))

    # Output records
    if len(records) > 0:
        write_records(records, baseheader, defline_rename, out, bool_nosortbylen)
    if args.discarded is not None and len(discarded) > 0:
        outfile = os.path.abspath(os.path.expanduser(args.discarded))
        write_records(discarded, baseheader, "discarded", outfile, bool_nosortbylen)


if __name__ == "__main__":
    main()
