#!/usr/bin/env python3
# Purpose: Print variants where nhemi_joint is not 0
# Mireia Marin i Ginestar

import hail as hl
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python print_hemizygous_variants.py <vcf_file_path>")
        sys.exit(1)

    # ------------------ INIT ------------------
    hl.init(log="/tmp/hail_print_hemizygous.log")

    # ------------------ READ ------------------
    vcf_path = sys.argv[1]
    print(f"Reading VCF: {vcf_path}")
    mt = hl.import_vcf(vcf_path, reference_genome='GRCh38')

    # ------------------ FILTER FOR HEMIZYGOUS VARIANTS ------------------
    # Check if nhemi_joint exists
    if 'nhemi_joint' not in mt.info:
        print("Error: 'nhemi_joint' field not found in INFO. Please run add_hemizygotes.py first.")
        sys.exit(1)

    # Filter for variants where nhemi_joint is not 0
    mt_hemi = mt.filter_rows(
        hl.is_defined(mt.info.nhemi_joint) & (mt.info.nhemi_joint != 0)
    )

    # Count how many variants have hemizygous calls
    n_hemi_variants = mt_hemi.count_rows()
    print(f"\nTotal variants with nhemi_joint != 0: {n_hemi_variants}")

    if n_hemi_variants == 0:
        print("No variants found with non-zero nhemi_joint.")
        hl.stop()
        return

    # ------------------ COLLECT AND PRINT 10 VARIANTS ------------------
    # Take first 10 variants
    n_to_show = min(10, n_hemi_variants)
    print(f"\nShowing first {n_to_show} variants:\n")

    # Select relevant fields for display
    rows = mt_hemi.rows()
    rows = rows.select(
        'locus',
        'alleles',
        nhemi_joint=rows.info.nhemi_joint,
        AC_joint_XY=rows.info.AC_joint_XY[0] if 'AC_joint_XY' in rows.info else hl.missing('int32'),
        AN_joint_XY=rows.info.AN_joint_XY if 'AN_joint_XY' in rows.info else hl.missing('int32')
    )

    # Collect first 10 rows
    sample_rows = rows.take(n_to_show)

    # Print header
    print(f"{'Locus':<20} {'Ref':<10} {'Alt':<10} {'nhemi_joint':<12} {'AC_XY':<10} {'AN_XY':<10}")
    print("-" * 80)

    # Print each variant
    for row in sample_rows:
        locus = f"{row.locus.contig}:{row.locus.position}"
        ref = row.alleles[0]
        alt = row.alleles[1] if len(row.alleles) > 1 else "."
        nhemi = row.nhemi_joint
        ac_xy = row.AC_joint_XY if hl.is_defined(row.AC_joint_XY) else "NA"
        an_xy = row.AN_joint_XY if hl.is_defined(row.AN_joint_XY) else "NA"
        
        print(f"{locus:<20} {ref:<10} {alt:<10} {nhemi:<12} {ac_xy:<10} {an_xy:<10}")

    # ------------------ ADDITIONAL STATS ------------------
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    # Calculate some statistics
    stats = rows.aggregate(hl.agg.stats(rows.nhemi_joint))
    print(f"\nnhemi_joint statistics:")
    print(f"  Mean:   {stats.mean:.2f}")
    print(f"  Median: {stats.median:.2f}" if hasattr(stats, 'median') else "  Median: N/A")
    print(f"  Min:    {stats.min}")
    print(f"  Max:    {stats.max}")
    print(f"  StdDev: {stats.stdev:.2f}")

    # ------------------ STOP ------------------
    hl.stop()

if __name__ == "__main__":
    main()