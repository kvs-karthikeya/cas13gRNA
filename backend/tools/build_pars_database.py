#!/usr/bin/env python3
"""
build_pars_database.py
Cas13 gRNA Designer - build the experimental PARS table for NONCODE lncRNAs.

What it does (all on YOUR computer, nothing is uploaded anywhere):
  1. Reads the two GEO PARS files (S1 and V1, GM12878) from GSE50676.
  2. Keeps only transcripts with enough reads (default: average >= 1 read per base,
     the same rule used in the original PARS paper).
  3. Finds each transcript's sequence in reference FASTA files (Ensembl ncRNA,
     optionally RefSeq or your own) and checks the length equals the PARS length.
  4. Looks for that exact sequence in the NONCODE human FASTA and records the
     NONCODE ID. Only EXACT sequence matches are accepted.
  5. Writes the real PARS score for every position:  log2((V1+5)/(S1+5))
     Positions with zero reads in both enzymes are written as NA (no data).

Usage (cmd):
  python build_pars_database.py ^
      --s1 GSM1226157_norm.GM12878_S1.tab.gz ^
      --v1 GSM1226158_norm.GM12878_V1.tab.gz ^
      --noncode NONCODEv6_human.fa.gz ^
      --ref Homo_sapiens_GRCh38_ncrna_fa.gz ^
      --extra-ref xist_extra_reference.fa ^
      --out pars_real_GM12878.tsv.gz --report match_report.txt

Add more reference FASTA files with extra --ref options (e.g. RefSeq NR_ sequences).
Headers must start with the transcript ID (version number is ignored).
"""
import argparse
import collections
import gzip
import hashlib
import math
import re
import sys
import time


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_fasta(path):
    """Yield (id, full header, sequence) from a (gzipped) FASTA file."""
    with open_text(path) as f:
        tid, hdr, parts = None, "", []
        for line in f:
            if line.startswith(">"):
                if tid is not None:
                    yield tid, hdr, "".join(parts)
                hdr = line[1:].strip()
                tid = hdr.split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if tid is not None:
            yield tid, hdr, "".join(parts)


def gene_from_header(hdr, label):
    """Ensembl headers carry 'gene_symbol:XIST'; for your own FASTA the 2nd word is used."""
    m = re.search(r"gene_symbol:(\S+)", hdr)
    if m:
        return m.group(1)
    words = hdr.split()
    if label == "extra" and len(words) > 1:
        return words[1]
    return ""


def norm_seq(seq):
    return seq.upper().replace("T", "U")


def seq_key(seq):
    return hashlib.md5(norm_seq(seq).encode()).digest()


def base_id(tid):
    """Drop a trailing .version (ENST00000429829.6 -> ENST00000429829)."""
    return tid.split(".")[0]


def pars_stats(path):
    """First pass: length and total reads per transcript (no big arrays kept)."""
    stats = {}
    with open_text(path) as f:
        for line in f:
            tid, arr = line.rstrip("\n").split("\t", 1)
            vals = arr.split(";")
            total = 0.0
            for x in vals:
                if x != "0":
                    total += float(x)
            stats[tid] = (len(vals), total)
    return stats


def pars_arrays(path, wanted):
    """Second pass: full per-position counts for the wanted transcripts only."""
    out = {}
    with open_text(path) as f:
        for line in f:
            tid, arr = line.split("\t", 1)
            if tid in wanted:
                out[tid] = [float(x) for x in arr.rstrip("\n").split(";")]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s1", required=True, help="GEO S1 file (.tab.gz)")
    ap.add_argument("--v1", required=True, help="GEO V1 file (.tab.gz)")
    ap.add_argument("--noncode", required=True, help="NONCODE human FASTA (.fa.gz)")
    ap.add_argument("--ref", action="append", default=[], help="reference FASTA with transcript IDs (repeatable)")
    ap.add_argument("--extra-ref", action="append", default=[],
                    help="your own FASTA, e.g. XIST under its RefSeq ID (repeatable; flagged in the report)")
    ap.add_argument("--min-avg-reads", type=float, default=1.0,
                    help="minimum average reads per base, S1+V1 (default 1.0)")
    ap.add_argument("--pseudocount", type=float, default=5.0)
    ap.add_argument("--skip-prefix", action="append", default=None,
                    help="ignore PARS transcripts whose ID starts with this (default: NM_, because "
                         "NM_ are protein-coding mRNAs and NONCODE only holds lncRNAs)")
    ap.add_argument("--out", default="pars_real_GM12878.tsv.gz")
    ap.add_argument("--report", default="match_report.txt")
    args = ap.parse_args()

    t0 = time.time()

    def log(msg):
        print("[%5.0fs] %s" % (time.time() - t0, msg), flush=True)

    # 1) which PARS transcripts have enough reads?
    log("Reading PARS S1 statistics ...")
    s_stats = pars_stats(args.s1)
    log("Reading PARS V1 statistics ...")
    v_stats = pars_stats(args.v1)
    if set(s_stats) != set(v_stats):
        sys.exit("ERROR: S1 and V1 files do not contain the same transcripts.")
    skip = tuple(args.skip_prefix) if args.skip_prefix else ("NM_",)
    candidates = {}
    for tid, (n, s_tot) in s_stats.items():
        if tid.startswith(skip):
            continue
        v_n, v_tot = v_stats[tid]
        if v_n != n:
            sys.exit("ERROR: length mismatch between S1 and V1 for " + tid)
        avg = (s_tot + v_tot) / n if n else 0.0
        if avg >= args.min_avg_reads:
            candidates[tid] = (n, avg)
    log("Transcripts in PARS data: %d; ignoring prefixes %s; with average >= %.1f reads/base: %d"
        % (len(s_stats), skip, args.min_avg_reads, len(candidates)))

    # 2) reference sequences for the candidates only
    wanted = {}
    for tid in candidates:
        wanted.setdefault(base_id(tid), []).append(tid)
    ref_seq = {}  # pars id -> (sequence, source label)
    for label, paths in (("reference", args.ref), ("extra", args.extra_ref)):
        for p in paths:
            log("Reading %s FASTA: %s" % (label, p))
            for rid, hdr, seq in read_fasta(p):
                for pars_id in wanted.get(base_id(rid), ()):
                    if pars_id not in ref_seq or label == "extra":
                        ref_seq[pars_id] = (seq, label + ":" + p.split("/")[-1].split("\\")[-1],
                                            gene_from_header(hdr, label))
    log("Candidates with a reference sequence: %d" % len(ref_seq))

    # 3) NONCODE sequences -> IDs (hash of the sequence keeps memory small)
    need_keys = {}
    status = {}
    for tid, (seq, _src, _gene) in ref_seq.items():
        n = candidates[tid][0]
        if len(seq) != n:
            status[tid] = "length_mismatch (reference %d vs PARS %d)" % (len(seq), n)
        else:
            need_keys[tid] = seq_key(seq)
    key_to_tids = collections.defaultdict(list)
    for tid, k in need_keys.items():
        key_to_tids[k].append(tid)
    log("Reading NONCODE FASTA (this takes a minute) ...")
    found = {}
    nc_seq = {}
    n_nc = 0
    for nid, _hdr, seq in read_fasta(args.noncode):
        n_nc += 1
        k = seq_key(seq)
        if k in key_to_tids:
            for tid in key_to_tids[k]:
                found.setdefault(tid, []).append(nid)
                nc_seq[tid] = norm_seq(seq)
    log("NONCODE sequences read: %d" % n_nc)

    # 4) final status and output
    matched = {}
    for tid in candidates:
        if tid in status:
            continue
        if tid not in ref_seq:
            status[tid] = "no_reference_sequence"
        elif tid not in found:
            status[tid] = "sequence_not_in_NONCODE"
        else:
            matched[tid] = sorted(found[tid])
            status[tid] = "matched"

    log("Reading PARS counts for %d matched transcripts ..." % len(matched))
    s_arr = pars_arrays(args.s1, set(matched))
    v_arr = pars_arrays(args.v1, set(matched))
    pc = args.pseudocount
    with gzip.open(args.out, "wt") as out:
        out.write("noncode_id\tsource_id\tgene_name\tlength\tavg_reads_per_base\tpct_positions_with_reads"
                  "\treference_source\tsequence\tpars_scores\ttotal_reads\n")
        for tid in sorted(matched):
            s, v = s_arr[tid], v_arr[tid]
            scores, reads, covered = [], [], 0
            for a, b in zip(s, v):
                tot = a + b
                reads.append("%g" % tot)
                if tot > 0:
                    covered += 1
                    scores.append("%.3f" % math.log2((b + pc) / (a + pc)))
                else:
                    scores.append("NA")
            n = len(s)
            out.write("%s\t%s\t%s\t%d\t%.2f\t%.1f\t%s\t%s\t%s\t%s\n" % (
                ",".join(matched[tid]), tid, ref_seq[tid][2], n, candidates[tid][1],
                100.0 * covered / n, ref_seq[tid][1], nc_seq[tid], ";".join(scores), ";".join(reads)))
    log("Wrote %s" % args.out)

    # 5) report
    def kind(t):
        return "NM_ (mRNA)" if t.startswith("NM_") else "NR_ (RefSeq non-coding)" if t.startswith("NR_") \
            else "ENST (Ensembl)" if t.startswith("ENST") else "other"

    by = collections.defaultdict(collections.Counter)
    for tid, st in status.items():
        by[kind(tid)][st.split(" ")[0]] += 1
    with open(args.report, "w") as r:
        r.write("PARS -> NONCODE matching report\n")
        r.write("Transcripts in PARS data: %d\n" % len(s_stats))
        r.write("Ignored ID prefixes: %s\n" % (skip,))
        r.write("With average >= %.1f reads/base: %d\n" % (args.min_avg_reads, len(candidates)))
        r.write("Matched to a NONCODE ID by exact sequence: %d\n\n" % len(matched))
        for k in sorted(by):
            r.write("%s: %s\n" % (k, dict(by[k])))
        r.write("\nMatched transcripts that used an --extra-ref sequence (identity NOT independently verified):\n")
        for tid in sorted(matched):
            if ref_seq[tid][1].startswith("extra"):
                r.write("  %s -> %s (length %d)\n" % (tid, ",".join(matched[tid]), candidates[tid][0]))
        multi = [t for t in matched if len(matched[t]) > 1]
        r.write("\nTranscripts whose sequence is shared by several NONCODE IDs: %d\n" % len(multi))
        for tid in multi[:20]:
            r.write("  %s -> %s\n" % (tid, ",".join(matched[tid])))
    log("Wrote %s" % args.report)
    print("\nDONE. Matched %d transcripts. See %s" % (len(matched), args.report))


if __name__ == "__main__":
    main()
