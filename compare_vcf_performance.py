#!/usr/bin/env python3


"""
compare_vcf_performance.py

Compares the variants produced by a bioinformatics pipeline against an 
external "truth" VCF, computing accuracy metrics.


HOW TO USE

Basic usage with required arguments:

    python compare_vcf_performance.py \
        --truth true_variants.vcf \
        --pred pipeline_variants.vcf

With a BED file restricting evaluation to high-confidence regions:

    python compare_vcf_performance.py \
        --truth true_variants.vcf \
        --pred pipeline_variants.vcf \
        --bed highconf_regions.bed

Saving global results to JSON and CSV:

    python compare_vcf_performance.py \
        --truth true_variants.vcf \
        --pred pipeline_variants.vcf \
        --save-json output.json \
        --save-csv output.csv

Saving per-chromosome results:

    python compare_vcf_performance.py \
        --truth true_variants.vcf \
        --pred pipeline_variants.vcf \
        --per-chrom \
        --save-perchrom-csv output_per_chrom.csv

Combining all options (truth VCF, predicted VCF, BED restriction, global 
output files, and per-chromosome statistics):

    python compare_vcf_performance.py \
        --truth true_variants.vcf \
        --pred pipeline_variants.vcf \
        --bed highconf_regions.bed \
        --save-json output.json \
        --save-csv output.csv \
        --per-chrom \
        --save-perchrom-csv output_per_chrom.csv

"""


import argparse
import csv
import gzip
import json
import os
import shutil
import tempfile
import vcfpy
from typing import Set, Tuple


def main():
    parser = argparse.ArgumentParser(description="Compare pipeline VCF against truth set.")
    parser.add_argument("--truth", required=True, help="Truth VCF (expected variants)")
    parser.add_argument("--pred", required=True, help="Pipeline-produced VCF")
    parser.add_argument("--bed", help="Optional BED file with valid regions (e.g. high-confidence)")
    parser.add_argument("--per-chrom", action="store_true", help="Report per-chromosome stats")
    parser.add_argument("--save-json", help="Path to save global results as JSON")
    parser.add_argument("--save-csv", help="Path to save global results as CSV")
    parser.add_argument("--save-perchrom-csv", help="Path to save per-chromosome results as CSV")

    args = parser.parse_args()

    # Load variants
    truth = load_vcf_as_set(args.truth, restrict_regions=args.bed)
    pred = load_vcf_as_set(args.pred, restrict_regions=args.bed)

    print(f"Loaded truth variants: {len(truth)}")
    print(f"Loaded predicted variants: {len(pred)}")
    print("Running comparison...\n")

    # Global metrics
    tp, fp, fn, precision, recall, f1 = compute_metrics(truth, pred)

    print("==== GLOBAL PERFORMANCE ====")
    print(f"True Positives (TP): {tp}")
    print(f"False Positives (FP): {fp}")
    print(f"False Negatives (FN): {fn}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1 Score:  {f1:.4f}")

    # Save global results to JSON
    if args.save_json:
        global_dict = {
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "truth_variants": len(truth),
            "predicted_variants": len(pred)
        }
        os.makedirs(os.path.dirname(args.save_json), exist_ok=True)
        with open(args.save_json, "w") as f:
            json.dump(global_dict, f, indent=4)

    # Save global results to CSV
    if args.save_csv:
        os.makedirs(os.path.dirname(args.save_csv), exist_ok=True)
        with open(args.save_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["metric", "value"])
            writer.writerow(["TP", tp])
            writer.writerow(["FP", fp])
            writer.writerow(["FN", fn])
            writer.writerow(["precision", precision])
            writer.writerow(["recall", recall])
            writer.writerow(["f1", f1])
            writer.writerow(["truth_variants", len(truth)])
            writer.writerow(["predicted_variants", len(pred)])

    # Per-chromosome metrics
    if args.per_chrom:
        print("\n==== PER-CHROMOSOME ====")
        chroms = set([v[0] for v in truth] + [v[0] for v in pred])
        for chrom in sorted(chroms):
            truth_c = {v for v in truth if v[0] == chrom}
            pred_c  = {v for v in pred if v[0] == chrom}
            tp, fp, fn, precision, recall, f1 = compute_metrics(truth_c, pred_c)

            print(f"\n--- {chrom} ---")
            print(f"TP={tp}, FP={fp}, FN={fn}, Precision={precision:.3f}, Recall={recall:.3f}, F1={f1:.3f}")
            # Save per-chromosome results as CSV
            if args.save_perchrom_csv:
                # Append mode so all chromosomes go into one file
                write_header = not os.path.exists(args.save_perchrom_csv)
                with open(args.save_perchrom_csv, "a", newline="") as f:
                    writer = csv.writer(f)
                    if write_header:
                        writer.writerow(["chrom", "TP", "FP", "FN", "precision", "recall", "f1"])
                    writer.writerow([chrom, tp, fp, fn, precision, recall, f1])
    
    return None


def compute_metrics(truth, pred):
    tp = len(truth & pred)
    fp = len(pred - truth)
    fn = len(truth - pred)

    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0

    return tp, fp, fn, precision, recall, f1


def load_vcf_as_set(path: str, restrict_regions=None) -> Set[Tuple[str, int, str, str]]:
    """
    Load variants from a VCF into a set of (CHROM, POS, REF, ALT).
    Supports optional region restriction using a BED file.
    """
    variants = set()
    allowed = None

    # If BED provided, load valid regions
    if restrict_regions:
        allowed = []
        with open(restrict_regions) as bed:
            for line in bed:
                if line.strip():
                    chrom, start, end = line.split()[:3]
                    allowed.append((chrom, int(start), int(end)))

    clean_path = sanitize_vcf(path)
    reader = vcfpy.Reader.from_path(clean_path)

    for rec in reader:
        pos = rec.POS
        chrom = rec.CHROM

        if allowed:
            in_region = any(chrom == r[0] and r[1] <= pos <= r[2] for r in allowed)
            if not in_region:
                continue

        ref = rec.REF.upper()

        for alt in rec.ALT:
            if alt.type != "SNV" and alt.type != "INDEL":
                continue
            alt_string = alt.value.upper()
            variants.add((chrom, pos, ref, alt_string))

    return variants


def sanitize_vcf(path):
    # Determine if file is gzipped
    gzipped = path.endswith(".gz")

    # Open input accordingly
    opener = gzip.open if gzipped else open
    mode = "rt" if gzipped else "r"

    tmp = tempfile.NamedTemporaryFile(delete=False, mode="wt")
    tmp_path = tmp.name

    with opener(path, mode) as f_in, tmp:
        for line in f_in:
            # Fix malformed ##fileDate lines
            if line.startswith("##fileDate") and "=" not in line:
                tmp.write("##fileDate=unknown\n")
            else:
                tmp.write(line)

    # If original was gzip, re-gzip the sanitized file
    if gzipped:
        gz_path = tmp_path + ".gz"
        with open(tmp_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(tmp_path)
        return gz_path

    return tmp_path


if __name__ == "__main__":
    main()
