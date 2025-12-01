# Purpose: Import a VCF and annotate hemizygous and heterozygous counts into INFO using gnomAD-style per-group AC fields.
# Mireia Marin i Ginestar

import hail as hl
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python add_hemizygotes_and_heterozygotes.py <file_path>")
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

    mt = mt.filter_rows(~hl.any(lambda f: hl.literal(filters_to_remove).contains(f), mt.filters))

    n_variants_after = mt.count_rows()
    print(f"Variants after filtering:  {n_variants_after}")
    print(f"Variants removed:          {n_variants_before - n_variants_after}")

    # ------------------ DEFINE GROUPS ------------------
    ancestries = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'remaining', 'sas']

    # ------------------ REGION LOGIC ------------------
    on_sex_chrom_nonpar = ~mt.locus.in_autosome_or_par()
    in_par = mt.locus.in_x_par() | mt.locus.in_y_par()

    annotations = {}

    # safe access to INFO field with default 0
    def get_info_int(field_name, index=None):
        """Return mt.info[field_name][index] if defined, else 0."""
        if not hasattr(mt.info, field_name):
            return hl.int32(0)
        field = mt.info[field_name]
        if index is not None:
            return hl.if_else(hl.is_defined(field), field[index], hl.int32(0))
        return hl.if_else(hl.is_defined(field), field, hl.int32(0))

    # calculate nhet = AC - 2*nhomalt
    def calc_nhet(ac_field, nhomalt_field):
        """Calculate heterozygous count: AC - 2*nhomalt"""
        ac = get_info_int(ac_field, index=0)
        nhomalt = get_info_int(nhomalt_field, index=0)
        return ac - 2 * nhomalt

    # -------- Overall / Joint --------
    # Non-PAR: nhemi for XY, nhet for XX
    # PAR: nhet for both
    
    if hasattr(mt.info, 'AC_joint_XY'):
        annotations['nhemi_joint_XY'] = hl.if_else(
            on_sex_chrom_nonpar,
            get_info_int('AC_joint_XY', 0),
            hl.int32(0)
        )
    
    if hasattr(mt.info, 'AC_joint_XX') and hasattr(mt.info, 'nhomalt_joint_XX'):
        # XX heterozygotes (non-PAR only)
        annotations['nhet_joint_XX'] = hl.if_else(
            on_sex_chrom_nonpar,
            calc_nhet('AC_joint_XX', 'nhomalt_joint_XX'),
            hl.int32(0)
        )
    
    # PAR heterozygotes (combined)
    if hasattr(mt.info, 'AC_joint') and hasattr(mt.info, 'nhomalt_joint'):
        annotations['nhet_joint'] = hl.if_else(
            in_par,
            calc_nhet('AC_joint', 'nhomalt_joint'),
            hl.int32(0)
        )
    
    # -------- TOTALS (across entire chromosome) --------
    # nhemi_total: hemizygotes in non-PAR (XY only)
    if hasattr(mt.info, 'AC_joint_XY'):
        annotations['nhemi_joint'] = hl.if_else(
            on_sex_chrom_nonpar,
            get_info_int('AC_joint_XY', 0),
            hl.int32(0)
        )
    
    # nhet_total: heterozygotes everywhere
    # Non-PAR: only XX heterozygotes (females)
    # PAR: combined heterozygotes (both sexes)
    if hasattr(mt.info, 'AC_joint') and hasattr(mt.info, 'nhomalt_joint'):
        if hasattr(mt.info, 'AC_joint_XX') and hasattr(mt.info, 'nhomalt_joint_XX'):
            annotations['nhet_joint'] = hl.case()\
                .when(on_sex_chrom_nonpar, calc_nhet('AC_joint_XX', 'nhomalt_joint_XX'))\
                .when(in_par, calc_nhet('AC_joint', 'nhomalt_joint'))\
                .default(hl.int32(0))
        else:
            # Fallback if XX-specific fields don't exist
            annotations['nhet_joint'] = hl.if_else(
                in_par,
                calc_nhet('AC_joint', 'nhomalt_joint'),
                hl.int32(0)
            )

    # -------- Ancestry-stratified --------
    for anc in ancestries:
        prefix = f'joint_{anc}'
        
        # Non-PAR: hemizygotes for XY
        ac_xy = f'AC_{prefix}_XY'
        if hasattr(mt.info, ac_xy):
            annotations[f'nhemi_{anc}_XY'] = hl.if_else(
                on_sex_chrom_nonpar,
                get_info_int(ac_xy, 0),
                hl.int32(0)
            )
        
        # Non-PAR: heterozygotes for XX
        ac_xx = f'AC_{prefix}_XX'
        nhomalt_xx = f'nhomalt_{prefix}_XX'
        if hasattr(mt.info, ac_xx) and hasattr(mt.info, nhomalt_xx):
            annotations[f'nhet_{anc}_XX'] = hl.if_else(
                on_sex_chrom_nonpar,
                calc_nhet(ac_xx, nhomalt_xx),
                hl.int32(0)
            )
        
        # PAR: heterozygotes (combined)
        ac_combined = f'AC_{prefix}'
        nhomalt_combined = f'nhomalt_{prefix}'
        if hasattr(mt.info, ac_combined) and hasattr(mt.info, nhomalt_combined):
            annotations[f'nhet_{anc}'] = hl.if_else(
                in_par,
                calc_nhet(ac_combined, nhomalt_combined),
                hl.int32(0)
            )
        
        # -------- Ancestry totals --------
        # nhemi_total for this ancestry
        if hasattr(mt.info, ac_xy):
            annotations[f'nhemi_{anc}_joint'] = hl.if_else(
                on_sex_chrom_nonpar,
                get_info_int(ac_xy, 0),
                hl.int32(0)
            )
        
        # nhet_total for this ancestry
        if hasattr(mt.info, ac_combined) and hasattr(mt.info, nhomalt_combined):
            if hasattr(mt.info, ac_xx) and hasattr(mt.info, nhomalt_xx):
                annotations[f'nhet_{anc}_joint'] = hl.case()\
                    .when(on_sex_chrom_nonpar, calc_nhet(ac_xx, nhomalt_xx))\
                    .when(in_par, calc_nhet(ac_combined, nhomalt_combined))\
                    .default(hl.int32(0))
            else:
                annotations[f'nhet_{anc}_joint'] = hl.if_else(
                    in_par,
                    calc_nhet(ac_combined, nhomalt_combined),
                    hl.int32(0)
                )

    # ------------------ APPLY ANNOTATIONS ------------------
    if annotations:
        mt = mt.annotate_rows(info=mt.info.annotate(**annotations))
        print(f"Successfully added {len(annotations)} annotations")
        print("Fields added:", sorted(annotations.keys()))
    else:
        print("Warning: No annotations were added (no matching fields found).")

    # ------------------ WRITE ------------------
    out_path = f"{vcf_path}.annotated-both-het-hemi.vcf.bgz"
    hl.export_vcf(mt, out_path)
    print(f"Wrote: {out_path}")

    # ------------------ STOP ------------------
    hl.stop()

if __name__ == "__main__":
    main()
