import hail as hl
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python unified_annotator.py <vcf_path>")
        sys.exit(1)

    # ------------------ INIT ------------------
    hl.init(log="/tmp/hail_unified_annotator.log")

    # ------------------ READ ------------------
    vcf_path = sys.argv[1]
    mt = hl.import_vcf(vcf_path, reference_genome='GRCh37')

    # ------------------ FILTER VARIANTS ------------------
    filters_to_remove = ['AC0', 'InbreedingCoeff', 'RF']
    
    n_variants_before = mt.count_rows()
    print(f"Variants before filtering: {n_variants_before}")
    
    mt = mt.filter_rows(
        ~hl.any(lambda f: hl.literal(filters_to_remove).contains(f), mt.filters)
    )
    
    n_variants_after = mt.count_rows()
    print(f"Variants after filtering:  {n_variants_after}")
    print(f"Variants removed:          {n_variants_before - n_variants_after}")

    # ------------------ REGION LOGIC ------------------
    is_autosome = mt.locus.in_autosome()
    is_chrx = mt.locus.contig == 'X'
    is_chry = mt.locus.contig == 'Y'
    in_par = mt.locus.in_x_par() | mt.locus.in_y_par()
    on_sex_chrom_nonpar = (is_chrx | is_chry) & ~in_par

    # ------------------ FIND ALL AC FIELDS ------------------
    ac_fields = [f for f in mt.info if f.startswith('AC')]
    print(f"\nFound {len(ac_fields)} AC fields to process:")
    print(ac_fields)

    # ------------------ HELPER FUNCTIONS ------------------
    def get_zero_value(field_value):
        """Create a zero value matching the field type"""
        field_type = field_value.dtype
        if isinstance(field_type, hl.tarray):
            return hl.map(lambda x: hl.int32(0), field_value)
        return hl.int32(0)

    def calc_nhet_base(ac_value, nhomalt_value):
        """Calculate base heterozygous count: AC - 2*nhomalt"""
        field_type = ac_value.dtype
        if isinstance(field_type, hl.tarray):
            return ac_value - (nhomalt_value * 2)
        return ac_value - (2 * nhomalt_value)

    # ------------------ STEP 1: CREATE NHEMI ANNOTATIONS FIRST ------------------
    nhemi_annotations = {}

    for ac_field in ac_fields:
        nhemi_field = ac_field.replace('AC', 'nhemi')
        
        if ac_field not in mt.info:
            continue
            
        ac_value = mt.info[ac_field]
        zero_value = get_zero_value(ac_value)
        
        if '_male' in ac_field:
            # male-specific fields: nhemi for chrX/chrY non-PAR
            nhemi_annotations[nhemi_field] = hl.if_else(
                on_sex_chrom_nonpar,
                ac_value,
                zero_value
            )
            print(f"Creating {nhemi_field} from {ac_field}")

        elif '_female' in ac_field:
            # female-specific fields: no hemizygotes
            nhemi_annotations[nhemi_field] = zero_value
            print(f"Creating {nhemi_field} from {ac_field}")

        else:
            # Combined/general fields: use AC_male for sex chromosomes
            ac_male_field = ac_field + '_male'
            if ac_male_field in mt.info:
                nhemi_annotations[nhemi_field] = hl.case()\
                    .when(is_chry, mt.info[ac_male_field])\
                    .when(is_chrx & ~in_par, mt.info[ac_male_field])\
                    .default(zero_value)
                print(f"Creating {nhemi_field} from {ac_field} (using {ac_male_field})")
            else:
                # Fallback if AC_male doesn't exist
                nhemi_annotations[nhemi_field] = zero_value
                print(f"Warning: {nhemi_field} set to 0 (missing {ac_male_field})")

    # Apply nhemi annotations first
    if nhemi_annotations:
        mt = mt.annotate_rows(info=mt.info.annotate(**nhemi_annotations))
        print(f"\nApplied {len(nhemi_annotations)} nhemi annotations")
        
        # Recalculate region logic for the new MatrixTable
        is_autosome = mt.locus.in_autosome()
        is_chrx = mt.locus.contig == 'X'
        is_chry = mt.locus.contig == 'Y'
        in_par = mt.locus.in_x_par() | mt.locus.in_y_par()
        on_sex_chrom_nonpar = (is_chrx | is_chry) & ~in_par

    # ------------------ STEP 2: CREATE NHET ANNOTATIONS (using nhemi) ------------------
    nhet_annotations = {}

    for ac_field in ac_fields:
        # Derive field names
        nhomalt_field = ac_field.replace('AC', 'nhomalt')
        nhet_field = ac_field.replace('AC', 'nhet')
        nhemi_field = ac_field.replace('AC', 'nhemi')
        
        # Check if required fields exist
        has_ac = ac_field in mt.info
        has_nhomalt = nhomalt_field in mt.info
        has_nhemi = nhemi_field in mt.info
        
        if not has_ac:
            continue
            
        ac_value = mt.info[ac_field]
        zero_value = get_zero_value(ac_value)
        
        # ========== NHET ANNOTATIONS ==========
        if has_nhomalt:
            nhomalt_value = mt.info[nhomalt_field]
            nhet_base = calc_nhet_base(ac_value, nhomalt_value)
            
            # Get nhemi value for correction
            nhemi_value = mt.info[nhemi_field] if has_nhemi else zero_value
            
            # Determine nhet based on chromosome and region
            if '_female' in ac_field:
                # female-specific fields: nhet for autosomes and chrX non-PAR
                nhet_annotations[nhet_field] = hl.case()\
                    .when(is_autosome, nhet_base)\
                    .when(is_chrx & ~in_par, nhet_base)\
                    .default(zero_value)
            
            elif '_male' in ac_field:
                # male-specific fields: no heterozygotes (all are hemizygous on sex chroms)
                # But can have hets on autosomes and PAR
                nhet_annotations[nhet_field] = hl.case()\
                    .when(is_autosome, nhet_base)\
                    .when(in_par, nhet_base)\
                    .default(zero_value)
            
            else:
                # Combined/general fields
                # On sex chromosomes non-PAR, subtract hemizygous counts
                nhet_annotations[nhet_field] = hl.case()\
                    .when(is_autosome, nhet_base)\
                    .when(in_par, nhet_base)\
                    .when(on_sex_chrom_nonpar, nhet_base - nhemi_value)\
                    .default(zero_value)
            
            print(f"Creating {nhet_field} from {ac_field} and {nhomalt_field}")
        else:
            # No nhomalt field available - set nhet to zero
            nhet_annotations[nhet_field] = zero_value
            print(f"Warning: {nhet_field} set to 0 (missing {nhomalt_field})")

    # Apply nhet annotations
    if nhet_annotations:
        mt = mt.annotate_rows(info=mt.info.annotate(**nhet_annotations))
        print(f"\nApplied {len(nhet_annotations)} nhet annotations")

    # ------------------ SUMMARY ------------------
    total_annotations = len(nhet_annotations) + len(nhemi_annotations)
    if total_annotations > 0:
        print(f"\n{'='*60}")
        print(f"Successfully added {total_annotations} annotations")
        print(f"  - {len(nhet_annotations)} nhet fields")
        print(f"  - {len(nhemi_annotations)} nhemi fields")
        print(f"{'='*60}")
        print("\nnhet fields added:", sorted(nhet_annotations.keys()))
        print("\nnhemi fields added:", sorted(nhemi_annotations.keys()))
    else:
        print("\nWarning: No annotations were added.")

    # ------------------ WRITE ------------------
    out_path = f"{vcf_path}.annotated_nhet_nhemi.vcf.bgz"
    hl.export_vcf(mt, out_path)
    print(f"\nWrote: {out_path}")

    # ------------------ WRITE MATRIXTABLE ------------------
    out_mt_path = f"{vcf_path}_nhet_nhemi.mt"
    mt.write(out_mt_path, overwrite=True)
    print(f"Wrote MatrixTable: {out_mt_path}")

    # ------------------ STOP ------------------
    hl.stop()

if __name__ == "__main__":
  main()
