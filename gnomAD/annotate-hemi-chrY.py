# Purpose: Create nhemi fields from AC fields for chrY variants
# Mireia Marin i Ginestar

import hail as hl
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python rename_chry_nhomalt.py <vcf_path>")
        sys.exit(1)

    # ------------------ INIT ------------------
    hl.init(log="/tmp/hail_rename_chry_nhomalt.log")

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

    # ------------------ IDENTIFY CHRY (AFTER FILTERING) ------------------
    is_chry = mt.locus.contig == 'chrY'

    # ------------------ FIND AC FIELDS ------------------
    # Get all info fields that start with 'AC_joint'
    ac_fields = [f for f in mt.info if f.startswith('AC_joint')]
    
    print(f"Found {len(ac_fields)} AC_joint fields to process:")
    print(ac_fields)

    # ------------------ CREATE NHEMI FIELDS ------------------
    new_annotations = {}
    
    # Process each AC field found
    for field in ac_fields:
        # Generate corresponding nhemi field name
        # AC_joint -> nhemi_joint
        # AC_joint_afr -> nhemi_joint_afr
        # AC_joint_afr_XY -> nhemi_joint_afr_XY
        nhemi_field = field.replace('AC', 'nhemi')
        
        # Check if nhemi field already exists
        nhemi_exists = hasattr(mt.info, nhemi_field)
        
        # Get the field value to determine its type
        field_value = mt.info[field]
        field_type = field_value.dtype
        
        # Create a zero value matching the field type
        if isinstance(field_type, hl.tarray):
            # For array types, create array of zeros with same length
            zero_value = hl.map(lambda x: hl.int32(0), field_value)
        else:
            # For scalar types
            zero_value = hl.int32(0)
        
        # If on chrY: nhemi = AC_joint value
        # If not on chrY: nhemi = existing value or 0
        if nhemi_exists:
            new_annotations[nhemi_field] = hl.if_else(
                is_chry,
                field_value,  # Copy AC_joint value to nhemi for chrY
                mt.info[nhemi_field]  # Keep existing nhemi value
            )
        else:
            new_annotations[nhemi_field] = hl.if_else(
                is_chry,
                field_value,  # Copy AC_joint value to nhemi for chrY
                zero_value  # Set to 0 (matching type) for non-chrY
            )

    # ------------------ APPLY ANNOTATIONS ------------------
    if new_annotations:
        mt = mt.annotate_rows(info=mt.info.annotate(**new_annotations))
        print(f"Successfully processed {len(ac_fields)} AC_joint fields")
        print(f"- For chrY: nhemi_joint_* = AC_joint_* values")
        print(f"- For non-chrY: nhemi_joint_* preserved or set to 0")
        nhemi_fields = sorted(set([k for k in new_annotations.keys()]))
        print(f"Created/updated {len(nhemi_fields)} nhemi_joint fields:")
        print(nhemi_fields)
    else:
        print("Warning: No AC_joint fields found to process.")

    # ------------------ WRITE ------------------
    out_path = f"{vcf_path}.chry_nhemi.vcf.bgz"
    hl.export_vcf(mt, out_path)
    print(f"Wrote: {out_path}")

    # ------------------ STOP ------------------
    hl.stop()

if __name__ == "__main__":
    main()s
