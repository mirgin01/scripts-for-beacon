#!/usr/bin/env python3
import sys, os, re, json, gzip

# Configurable source metadata
SOURCE = "The Genome Aggregation Database (gnomAD)"
SOURCE_REF = "gnomad.broadinstitute.org/"

ID_RE = re.compile(r'^##INFO=<ID=([^,>]+)')

def open_any(path):
    if path.lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_header_ids(vcf_path):
    ids = set()
    with open_any(vcf_path) as f:
        for line in f:
            if not line.startswith("##INFO=<ID="):
                # Header can be long; stop once we reach column line
                if line.startswith("#CHROM"):
                    break
                continue
            m = ID_RE.match(line)
            if m:
                ids.add(m.group(1))
    return ids

def build_populations(ids):
    # Collect prefixes appearing after AF_/AC_/AN_/nhom_[alt]_ / nhet_ / nhemi_
    prefixes = set()

    def add_prefix_from(pattern):
        for s in ids:
            m = re.match(pattern, s)
            if m:
                prefixes.add(m.group(1))

    add_prefix_from(r'^AF_(.+)$')
    add_prefix_from(r'^AC_(.+)$')
    add_prefix_from(r'^AN_(.+)$')
    add_prefix_from(r'^nhom_(?:alt_)?(.+)$')
    add_prefix_from(r'^nhet_(.+)$')
    add_prefix_from(r'^nhemi_(.+)$')

    pops = []
    for p in sorted(prefixes):
        # Build field names for this prefix, only if present in header
        af   = f"AF_{p}"   if f"AF_{p}"   in ids else None
        ac   = f"AC_{p}"   if f"AC_{p}"   in ids else None
        an   = f"AN_{p}"   if f"AN_{p}"   in ids else None
        nhom = (f"nhom_alt_{p}" if f"nhom_alt_{p}" in ids
                else f"nhom_{p}" if f"nhom_{p}" in ids else None)
        nhet = f"nhet_{p}" if f"nhet_{p}" in ids else None
        # nhemi is optional in your schema; collect if present
        nhemi = f"nhemi_{p}" if f"nhemi_{p}" in ids else None

        # Require at least AF & AC to include the population
        if not (af and ac):
            continue

        entry = {
            "population": p,  # use the detected prefix as the population name
            "alleleFrequency": af,
            "alleleCount": ac,
            "genotypeHomozygous": nhom,
            "genotypeHeterozygous": nhet,
            "genotypeHemizygous": nhemi,
            "alleleNumber": an
        }
        # If you also want to include hemizygous counts when present, uncomment:
        # if nhemi: entry["alleleCountHemizygous"] = nhemi

        pops.append(entry)

    return pops

def main(vcf_path):
    ids = parse_header_ids(vcf_path)
    populations = build_populations(ids)

    out = {
        "numberOfPopulations": len(populations),
        "source": SOURCE,
        "sourceReference": SOURCE_REF,
        "populations": populations
    }

    out_path = os.path.join(os.path.dirname(vcf_path), "populations.json")
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=4)
    print(f"Wrote {out_path} with {len(populations)} populations.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <vcf(.gz)>")
    main(sys.argv[1])

