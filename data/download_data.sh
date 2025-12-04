#!/bin/bash
# Download raw count matrices from GEO for Melms et al. 2021
# GSE171524: A molecular single-cell lung atlas of lethal COVID-19

set -e

echo "========================================"
echo "Downloading COVID-19 snRNA-seq Data"
echo "GEO Accession: GSE171524"
echo "========================================"

# Base URL for GEO supplementary files
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"

# Array of sample IDs
CONTROL_SAMPLES=(
  "GSM5226574_C51ctr"
  "GSM5226575_C52ctr"
  "GSM5226576_C53ctr"
  "GSM5226577_C54ctr"
  "GSM5226578_C55ctr"
  "GSM5226579_C56ctr"
  "GSM5226580_C57ctr"
)

COVID_SAMPLES=(
  "GSM5226581_L01cov"
  "GSM5226582_L03cov"
  "GSM5226583_L04cov"
  "GSM5226584_L04covaddon"
  "GSM5226585_L05cov"
  "GSM5226586_L06cov"
  "GSM5226587_L07cov"
  "GSM5226588_L08cov"
  "GSM5226589_L09cov"
  "GSM5226590_L10cov"
  "GSM5226591_L11cov"
  "GSM5226592_L12cov"
  "GSM5226593_L13cov"
  "GSM5226594_L15cov"
  "GSM5226595_L16cov"
  "GSM5226596_L17cov"
  "GSM5226597_L18cov"
  "GSM5226598_L19cov"
  "GSM5226599_L21cov"
  "GSM5226600_L22cov"
)

# Function to download and extract
download_sample() {
  local sample_id=$1
  local gsm=$(echo $sample_id | cut -d'_' -f1)
  local gsm_dir=$(echo $gsm | sed 's/GSM\([0-9]\{3\}\).*/GSM\1nnn/')
  
  local url="${GEO_BASE}/${gsm_dir}/${gsm}/suppl/${sample_id}_raw_counts.csv.gz"
  
  echo "Downloading: $sample_id"
  
  if [ -f "${sample_id}_raw_counts.csv" ]; then
    echo "  → Already exists, skipping"
    return
  fi
  
  wget -q "$url" -O "${sample_id}_raw_counts.csv.gz"
  
  if [ $? -eq 0 ]; then
    echo "  → Extracting..."
    gunzip "${sample_id}_raw_counts.csv.gz"
    echo "  ✓ Complete"
  else
    echo "  ✗ Failed to download"
  fi
}

# Download all samples
echo ""
echo "Downloading Control samples (7)..."
for sample in "${CONTROL_SAMPLES[@]}"; do
  download_sample "$sample"
done

echo ""
echo "Downloading COVID-19 samples (20)..."
for sample in "${COVID_SAMPLES[@]}"; do
  download_sample "$sample"
done

echo ""
echo "========================================"
echo "✓ Download Complete"
echo "========================================"
echo "Total files: $(ls -1 *_raw_counts.csv 2>/dev/null | wc -l)"
echo ""
echo "Note: If downloads failed, you may need to:"
echo "  1. Visit GEO directly: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524"
echo "  2. Download files manually from 'Supplementary file' section"
echo "  3. Extract .gz files to this directory"
