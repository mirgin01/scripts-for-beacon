# Purpose: Import a VCF and annotate hemizygous counts into INFO using gnomAD-style per-group AC fields.
# Mireia Marin i Ginestar

import hail as hl
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python add_hemizygotes.py <file_path>")
        sys.exit(1)

    # ------------------ INIT ------------------
    hl.init(log="/tmp/hail_add_hemizygotes.log")

    # ------------------ READ ------------------
    vcf_path = sys.argv[1]
    mt = hl.import_vcf(vcf_path, reference_genome='GRCh38')

    # ------------------ FILTER VARIANTS ------------------
    filters_to_remove = ['BOTH_FILTERED', 'EXOMES_FILTERED', 'GENOMES_FILTERED']

    n_variants_before = mt.count_rows()
    print(f"Variants before filtering: {n_variants_before}")

    # Keep rows where none of the unwanted FILTERs are present
    mt = mt.filter_rows(~hl.any(lambda f: hl.literal(filters_to_remove).contains(f), mt.filters))

    n_variants_after = mt.count_rows()
    print(f"Variants after filtering:  {n_variants_after}")
    print(f"Variants removed:          {n_variants_before - n_variants_after}")

    # ------------------ DEFINE GROUPS ------------------
    # Ancestries present in file
    ancestries = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'remaining', 'sas']

    # ------------------ HEMIZYGOUS LOGIC ------------------
    # Hemi counts only apply to sex chromosomes outside PAR
    on_sex_chrom_nonpar = ~mt.locus.in_autosome_or_par()

    he_annotations = {}

    # Add a safe nhemi = AC_XY[0] when non-PAR, else 0, if field exists
    def nhemi_from_ac(field_name):
        """Return an expression: AC_field[0] if defined and non-PAR sex chrom; else 0."""
        return hl.if_else(
            on_sex_chrom_nonpar & hl.is_defined(mt.info[field_name]),
            mt.info[field_name][0],  # Index the array here to get int32
            hl.int32(0)
        )

    # -------- Overall / Joint --------
    # nhemi_joint (combined sexes): equal to AC_joint_XY (male AC) on non-PAR sex chrom
    if hasattr(mt.info, 'AC_joint_XY'):
        he_annotations['nhemi_joint'] = nhemi_from_ac('AC_joint_XY')

    # Sex-stratified:
    # nhemi_joint_XX = 0 always; nhemi_joint_XY = AC_joint_XY (on non-PAR sex chrom)
    if hasattr(mt.info, 'AC_joint_XY'):
        he_annotations['nhemi_joint_XY'] = nhemi_from_ac('AC_joint_XY')
        he_annotations['nhemi_joint_XX'] = hl.int32(0)

    # -------- Ancestry-stratified (combined by sex inside ancestry) --------
    # For each ancestry, if AC_joint_{ancestry}_XY exists, that's the hemizygote count on non-PAR sex chrom
    for anc in ancestries:
        ac_xy_field = f'AC_joint_{anc}_XY'
        if hasattr(mt.info, ac_xy_field):
            he_annotations[f'nhemi_{anc}'] = nhemi_from_ac(ac_xy_field)

    # Sex-split hemizygotes within ancestries
    for anc in ancestries:
        ac_xy_field = f'AC_joint_{anc}_XY'
        if hasattr(mt.info, ac_xy_field):
            he_annotations[f'nhemi_{anc}_XY'] = nhemi_from_ac(ac_xy_field)
            he_annotations[f'nhemi_{anc}_XX'] = hl.int32(0)

    # ------------------ APPLY ANNOTATIONS ------------------
    if he_annotations:
        mt = mt.annotate_rows(info=mt.info.annotate(**he_annotations))
    else:
        print("Warning: No hemizygote annotations were added (no matching AC_*_XY fields found).")

    # ------------------ WRITE ------------------
    out_path = f"{vcf_path}.nhemi.vcf.bgz"
    hl.export_vcf(mt, out_path)
    print(f"Wrote: {out_path}")

    print(f"Successfully added {len(he_annotations)} hemizygote count fields")
    print("Fields added:", sorted(he_annotations.keys()))

    # ------------------ STOP ------------------
    hl.stop()

if __name__ == "__main__":
    main()
