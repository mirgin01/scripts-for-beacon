# Scripts to comfortably work with beacon 

## load-dataset-into-af-beacon.README

Instructions to correctly load a dataset into a beacon instance

## fill-in-populations.py

It will fill in populations.json file based on the header information. The current version is curated for gnomad38 so only major populations, females and males are inserted. 

### Chromosomes 1..2

Leave empty the hemizygous field: "alleleCountHemizygous": "",

### Chromosomes X and Y

Leave empty the heterozygous field: "alleleCountHeterozygous": "",

## gnomAD

### annotate-hemi-chrX.py

Import a chrX VCF and annotate hemizygous counts into INFO using gnomAD-style per-group AC fields.

Script Logic:

1. Hemizygous count = AC_XY[0] (the first allele count from males)

2. Uses `AC_XY[0] = nhemi_XY[0]` because on chrX, males are hemizygous (one X chromosome)

3. Sets nhemi_XX = 0 (females with two X chromosomes cannot be hemizygous)

### annotate-hemi-chrY.py

Import a chrY VCF and annotate hemizygous counts into INFO using gnomAD-style per-group AC fields.

Script Logic:

1. Hemizygous count = entire AC_joint value 

2. Uses the full AC_joint because on chrY, essentially all calls are hemizygous

3. Doesn't distinguish XY vs XX (chrY is male-only in practice)

### annotate-het-gnomad38.py

Given a gnomAD per-chromosome VCF, this script:

1. Removes variants flagged with BOTH_FILTERED, EXOMES_FILTERED, or GENOMES_FILTERED.

2. For each ancestry and sex stratum (and joint/raw), computes `nhet = AC − 2 * nhomalt`.

3. Adds those nhet_* fields into the INFO field.

4. Writes out a new filtered + annotated VCF for that chromosome.

### Extract_hemi_not_0.py

Extracts variants where hemizygous count is not 0 to benchmark the results against gnomad. 