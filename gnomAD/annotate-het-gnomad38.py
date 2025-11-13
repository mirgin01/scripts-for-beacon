import hail as hl
import sys

chr_num = sys.argv[1]

# Mireia Marin i Ginestar 

# Initialize Hail
hl.init()


# Read the VCF file
#mt = hl.import_vcf('../../gnomAD/gnomad38/gnomad.joint.v4.1.sites.chr1.vcf.bgz', reference_genome='GRCh38')
mt = hl.import_vcf(f'../../gnomAD/gnomad38/gnomad.joint.v4.1.sites.chr{chr_num}.vcf.bgz', reference_genome='GRCh38')

# ======================== FILTER VARIANTS ====================================================
# Filter out variants with specific FILTER values
filters_to_remove = ['BOTH_FILTERED', 'EXOMES_FILTERED', 'GENOMES_FILTERED']

# Count before filtering
n_variants_before = mt.count_rows()
print(f"Variants before filtering: {n_variants_before}")

# Filter: keep only variants where FILTER is NOT in the list
# filters is a set in Hail, so we check if any of the unwanted filters are present
mt = mt.filter_rows(
    ~hl.any(lambda f: hl.literal(filters_to_remove).contains(f), mt.filters)
)

# Count after filtering
n_variants_after = mt.count_rows()
print(f"Variants after filtering: {n_variants_after}")
print(f"Variants removed: {n_variants_before - n_variants_after}")

# ======================== CALCULATE STRATIFIED AFs ============================================

# Define ancestries
ancestries = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'remaining', 'sas']

# Define sex categories
sex_categories = ['', '_XX', '_XY']  # '' for combined, _XX for female, _XY for male

# Create annotation dictionary
annotations = {}

# Add nhet for joint (no ancestry suffix)
for sex in sex_categories:
    sex_suffix = sex if sex else ''
    ac_field = f'AC_joint{sex_suffix}'
    nhomalt_field = f'nhomalt_joint{sex_suffix}'
    nhet_field = f'nhet_joint{sex_suffix}'
    
    # Check if fields exist and calculate
    if hasattr(mt.info, ac_field) and hasattr(mt.info, nhomalt_field):
        annotations[nhet_field] = mt.info[ac_field] - 2 * mt.info[nhomalt_field]

# Add nhet for each ancestry
for ancestry in ancestries:
    for sex in sex_categories:
        sex_suffix = sex if sex else ''
        ac_field = f'AC_joint_{ancestry}{sex_suffix}'
        nhomalt_field = f'nhomalt_joint_{ancestry}{sex_suffix}'
        nhet_field = f'nhet_joint_{ancestry}{sex_suffix}'
        
        # Check if fields exist and calculate
        if hasattr(mt.info, ac_field) and hasattr(mt.info, nhomalt_field):
            annotations[nhet_field] = mt.info[ac_field] - 2 * mt.info[nhomalt_field]

# Add nhet for raw
if hasattr(mt.info, 'AC_joint_raw') and hasattr(mt.info, 'nhomalt_joint_raw'):
    annotations['nhet_joint_raw'] = mt.info.AC_joint_raw - 2 * mt.info.nhomalt_joint_raw

# Annotate the matrix table with all new fields
mt = mt.annotate_rows(info=mt.info.annotate(**annotations))

# Write the output VCF
hl.export_vcf(
        mt,
        f"gnomad.joint.v4.1.sites.chr{chr_num}.subset2_nhet.vcf.bgz"
        )


print(f"Successfully added {len(annotations)} heterozygote count fields")
print("Fields added:", list(annotations.keys()))

# Stop Hail
hl.stop()

