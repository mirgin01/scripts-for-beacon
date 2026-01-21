# script to convert all the files in a folder to BFF
# now configured for 1..22 --> Before transforming sex chromosomes update populations.json to read hemizygous counts
# between each VCF conversion the mongoDB is emptied. 1 VCF - 1 BFF relation

#!/bin/bash
SOURCE_DIR="/bio-scratch/mireia/beacon-datasets/gnomAD.v4/trying-loop"
DEST_DIR="/bio-scratch/mireia/beacon2_RI_Tools_v2/beacon-data-tools/files/vcf/files_to_read"
OUTPUT_DIR="/bio-scratch/mireia//beacon-datasets/gnomAD.v4/BFFs/"  # Ajusta esta ruta según necesites

# Create output directory if doesnt exist
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="$OUTPUT_DIR/conversion_log_$(date +%Y%m%d_%H%M%S).txt"

echo "Starting conversion from VCF to BFF - $(date)" | tee "$LOG_FILE"
echo "=============================================" | tee -a "$LOG_FILE"

# Loop for chr 1-22
for CHR in {1..22}; do
    echo "" | tee -a "$LOG_FILE"
    echo "Processing chromosome $CHR - $(date)" | tee -a "$LOG_FILE"
    echo "-------------------------------------------" | tee -a "$LOG_FILE"
   
    # Move VCF
    VCF_FILE="gnomad.joint.v4.1.sites.chr${CHR}._nhet_subset.vcf.gz"
    SOURCE_PATH="${SOURCE_DIR}/${VCF_FILE}"
    
    # Check the file exists
    if [ ! -f "$SOURCE_PATH" ]; then
        echo "ERROR: File not found: $SOURCE_PATH" | tee -a "$LOG_FILE"
        continue
    fi
    
    # 1. Move file
    echo "1. Moving file $VCF_FILE..." | tee -a "$LOG_FILE"
    mv "$SOURCE_PATH" "$DEST_DIR/" 2>&1 | tee -a "$LOG_FILE"
    if [ $? -ne 0 ]; then
        echo "ERROR: file couldnt be moved" | tee -a "$LOG_FILE"
        continue
    fi
    
    # 2. Convert to BFF
    echo "2. Executing genomicVariations_vcf.py..." | tee -a "$LOG_FILE"
    docker exec -it ri-tools python genomicVariations_vcf.py 2>&1 | tee -a "$LOG_FILE"
    if [ $? -ne 0 ]; then
        echo "ERROR: while executing container" | tee -a "$LOG_FILE"
        # Devolver el archivo a su ubicación original
        mv "${DEST_DIR}/${VCF_FILE}" "$SOURCE_PATH"
        continue
    fi
    
    # 3. Export JSON
    OUTPUT_JSON="$OUTPUT_DIR/genomicVariations_gnomadv4.chr${CHR}.json"
    echo "3. Export to JSON: $OUTPUT_JSON..." | tee -a "$LOG_FILE"
    docker exec ri-tools-mongo mongoexport --jsonArray \
        --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" \
        --collection genomicVariations > "$OUTPUT_JSON" 2>&1 | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: something went wrong while exporting to JSON" | tee -a "$LOG_FILE"
    else
        echo "Export complete: $(wc -l < "$OUTPUT_JSON"). N of lines:" | tee -a "$LOG_FILE"
    fi
    
    # 4. Clean MongoDB
    echo "4. Cleaning MongoDB..." | tee -a "$LOG_FILE"
    docker exec -i ri-tools-mongo mongosh <<EOF 2>&1 | tee -a "$LOG_FILE"
use admin
db.auth("root","example")
use beacon
db.genomicVariations.deleteMany({})
print("Delete variants. Total: " + db.genomicVariations.countDocuments({}))
EOF
    
    # Remove first and last lines (MongoDB metadata)
    echo "Removing MongoDB metadata from JSON..." | tee -a "$LOG_FILE"
    sed -i '1d;$d' "$OUTPUT_JSON"
    echo "Cleanup complete. Final file: $(wc -l < "$OUTPUT_JSON") lines" | tee -a "$LOG_FILE"

    # Devolver el archivo a su ubicación original (opcional)
    echo "5. Move VCF to original path..." | tee -a "$LOG_FILE"
    mv "${DEST_DIR}/${VCF_FILE}" "$SOURCE_PATH" 2>&1 | tee -a "$LOG_FILE"
    
    echo "Chromosome $CHR complete - $(date)" | tee -a "$LOG_FILE"
done

echo "" | tee -a "$LOG_FILE"
echo "=============================================" | tee -a "$LOG_FILE"
echo "Conversion finished - $(date)" | tee -a "$LOG_FILE"
echo "Log saved in: $LOG_FILE" | tee -a "$LOG_FILE"
