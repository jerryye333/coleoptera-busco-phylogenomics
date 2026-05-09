# filter_by_coverage.py

# Filter gappy sequences (< mincov coverage) from protein and corresponding nucleotide alignments. Writes filtered FASTA files (unless --report-only) and logs removed sequences (Locus,SeqID) to stdout in CSV format.



import os
import sys
import argparse
from Bio import SeqIO

def filter_alignment(in_prot, in_nuc, out_prot, out_nuc, mincov, report_only):
    # Read the protein alignment
    prot_records = list(SeqIO.parse(in_prot, "fasta"))
    if not prot_records:
        sys.stderr.write(f"Warning: no records found in {in_prot}\n")
        return

    # Determine alignment length from first record
    L = len(prot_records[0].seq)

    kept_prot = []
    removed_ids = []

    for rec in prot_records:
        seq_str = str(rec.seq)
        gap_count = seq_str.count('-') + seq_str.count('.')
        non_gap = L - gap_count
        cov = non_gap / L if L > 0 else 0.0

        if cov >= mincov:
            kept_prot.append(rec)
        else:
            removed_ids.append(rec.id)

    # Write filtered protein alignment (unless report-only)
    if not report_only:
        os.makedirs(os.path.dirname(out_prot), exist_ok=True)
        SeqIO.write(kept_prot, out_prot, "fasta")

    # Filter corresponding nucleotide alignment (unless report-only)
    if in_nuc and out_nuc:
        if os.path.isfile(in_nuc):
            nuc_records = list(SeqIO.parse(in_nuc, "fasta"))
            kept_nuc = [r for r in nuc_records if r.id not in removed_ids]
            if not report_only:
                os.makedirs(os.path.dirname(out_nuc), exist_ok=True)
                SeqIO.write(kept_nuc, out_nuc, "fasta")
        else:
            sys.stderr.write(f"Warning: nucleotide file not found: {in_nuc}\n")

    # Log removed sequences
    locus = os.path.basename(in_prot)
    if locus.endswith(".aa.fas"):
        locus = locus[:-len(".aa.fas")]
    elif locus.endswith(".fas"):
        locus = locus[:-4]

    for seqid in removed_ids:
        print(f"{locus},{seqid}")

def main():
    p = argparse.ArgumentParser(
        description="Filter sequences by coverage; remove those < mincov."
    )
    p.add_argument(
        "--in_prot", required=True,
        help="Input protein alignment (.aa.fas)"
    )
    p.add_argument(
        "--in_nuc", required=False,
        help="Input nucleotide alignment (.nuc.fas) (optional)"
    )
    p.add_argument(
        "--out_prot", required=True,
        help="Output filtered protein alignment (use '/dev/null' if report-only)"
    )
    p.add_argument(
        "--out_nuc", required=False,
        help="Output filtered nucleotide alignment (if --in_nuc used; use '/dev/null' if report-only)"
    )
    p.add_argument(
        "--mincov", type=float, default=0.7,
        help="Minimum coverage fraction to keep a sequence (default: 0.7)"
    )
    p.add_argument(
        "--report-only", action="store_true",
        help="Do not write any FASTA output; only print Locus,SeqID lines"
    )
    args = p.parse_args()

    if (args.in_nuc and not args.out_nuc) or (args.out_nuc and not args.in_nuc):
        p.error("Both --in_nuc and --out_nuc must be provided together.")

    filter_alignment(
        in_prot=args.in_prot,
        in_nuc=args.in_nuc,
        out_prot=args.out_prot,
        out_nuc=args.out_nuc,
        mincov=args.mincov,
        report_only=args.report_only
    )

if __name__ == "__main__":
    main()