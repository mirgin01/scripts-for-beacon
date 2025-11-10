#!/usr/bin/env python3
import sys
import os
import re
import json
import gzip
from argparse import ArgumentParser

"""
Input: a VCF with the het counts. 
Process:
1. parses the VCF/BCF header to discover population prefixes,
2. builds the JSON structure,
3. filters to joint* populations (keeping only joint_XX and joint_XY among sex splits),
4. fills genotypeHemizygous for every remaining population as nhemi_<population>,
Output: the final JSON (to a file if -o is given, otherwise to STDOUT).
"""

# --- Configurable source metadata ---
SOURCE = "The Genome Aggregation Database (gnomAD)"
SOURCE_REF = "gnomad.broadinstitute.org/"

ID_RE = re.compile(r'^##INFO=<ID=([^,>]+)')

def open_any(path):
    if path.lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_header_ids(vcf_path):
    """Return a set of INFO IDs from the VCF/BCF header."""
    ids = set()
    with open_any(vcf_path) as f:
        for line in f:
            if not line.startswith("##INFO=<ID="):
                # Stop early once we reach the column header
                if line.startswith("#CHROM"):
                    break
                continue
            m = ID_RE.match(line)
            if m:
                ids.add(m.group(1))
    return ids

def build_populations(ids):
    """
    Discover population 'prefixes' from AF_/AC_/AN_/nhom*/nhet_/nhemi_ keys,
    and build per-population entries with available fields.
    """
    prefixes = set()

    def add_prefix_from(pattern):
        for s in ids:
            m = re.match(pattern, s)
            if m:
                prefixes.add(m.group(1))

    # Gather candidate prefixes from common fields
    add_prefix_from(r'^AF_(.+)$')
    add_prefix_from(r'^AC_(.+)$')
    add_prefix_from(r'^AN_(.+)$')
    add_prefix_from(r'^nhom_(.+)$')      # some datasets use nhom_<p>
    add_prefix_from(r'^nhomalt_(.+)$')   # gnomAD-style nhomalt_<p>
    add_prefix_from(r'^nhet_(.+)$')
    add_prefix_from(r'^nhemi_(.+)$')

    pops = []
    for p in sorted(prefixes):
        af   = f"AF_{p}"     if f"AF_{p}"     in ids else None
        ac   = f"AC_{p}"     if f"AC_{p}"     in ids else None
        an   = f"AN_{p}"     if f"AN_{p}"     in ids else None
        # Prefer nhomalt_<p>, otherwise accept nhom_<p> if present
        nhom = (f"nhomalt_{p}" if f"nhomalt_{p}" in ids
                else (f"nhom_{p}" if f"nhom_{p}" in ids else None))
        nhet = f"nhet_{p}"   if f"nhet_{p}"   in ids else None
        # We'll force-fill genotypeHemizygous later; here we only detect presence
        nhemi_present = f"nhemi_{p}" in ids

        # Require at least AF & AC to include the population
        if not (af and ac):
            continue

        pops.append({
            "population": p,
            "alleleFrequency": af,
            "alleleCount": ac,
            "genotypeHomozygous": nhom,
            "genotypeHeterozygous": nhet,
            "genotypeHemizygous": (f"nhemi_{p}" if nhemi_present else None),  # temporary; overwritten later
            "alleleNumber": an
        })

    return pops

def filter_joint_populations(populations):
    """
    Keep only populations whose name starts with 'joint'.
    Among sex-splits, keep only joint_XX and joint_XY.
    """
    keep = []
    for pop in populations:
        name = pop["population"]
        if not name.startswith("joint"):
            continue
        # Allow only joint_XX and joint_XY among *_XX/*_XY
        if re.search(r"_XX$", name) or re.search(r"_XY$", name):
            if name not in ("joint_XX", "joint_XY"):
                continue
        keep.append(pop)
    return keep

def finalize_hemi(populations):
    """
    Set genotypeHemizygous to 'nhemi_<population>' for every entry,
    regardless of whether the field was present in the header.
    """
    for pop in populations:
        pop_name = pop["population"]
        pop["genotypeHemizygous"] = f"nhemi_{pop_name}"

def assemble_json(populations):
    return {
        "numberOfPopulations": len(populations),
        "source": SOURCE,
        "sourceReference": SOURCE_REF,
        "populations": populations
    }

def main():
    ap = ArgumentParser(description="Build final gnomAD joint populations JSON from a VCF(.gz) header.")
    ap.add_argument("vcf", help="Input VCF/BCF file (optionally gz/bgz compressed)")
    ap.add_argument("-o", "--output", help="Write final JSON to this file; default is STDOUT")
    args = ap.parse_args()

    ids = parse_header_ids(args.vcf)
    populations = build_populations(ids)
    populations = filter_joint_populations(populations)
    finalize_hemi(populations)
    out_obj = assemble_json(populations)

    if args.output:
        with open(args.output, "w") as fh:
            json.dump(out_obj, fh, indent=4)
        # Quiet by default; print a brief note to stderr if you like:
        # print(f"Wrote {args.output} with {len(populations)} joint populations.", file=sys.stderr)
    else:
        json.dump(out_obj, sys.stdout, indent=4)
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
