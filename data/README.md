# Data Directory

This directory contains **placeholders only**. Raw data files are not included in this repository due to size constraints.

## Data Sources

### Original Publication
**Melms et al., Nature 2021**  
"A molecular single-cell lung atlas of lethal COVID-19"  
DOI: [10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1)

### Raw Count Matrices
**GEO Accession:** [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524)

**Samples used in this analysis:**
- 7 Control samples: `GSM5226574_C51ctr` through `GSM5226580_C57ctr`
- 20 COVID-19 samples: `GSM5226581_L01cov` through `GSM5226600_L22cov`

### Processed Data (Optional)
Pre-integrated Seurat object and processed data available on Zenodo:
- [Original Study Data Repository](https://zenodo.org/record/4698034)

## Download Instructions

### Option 1: Manual Download from GEO
1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524
2. Navigate to "Supplementary file" section
3. Download individual `*_raw_counts.csv.gz` files
4. Extract to this `data/` directory

### Option 2: Automated Download Script
```bash
# From repository root
bash data/download_data.sh
```

This will download all 27 raw count matrices (~10 GB total).

## Expected Directory Structure

After downloading, your `data/` directory should look like:

```
data/
├── README.md (this file)
├── download_data.sh
├── GSM5226574_C51ctr_raw_counts.csv
├── GSM5226575_C52ctr_raw_counts.csv
├── ...
└── GSM5226600_L22cov_raw_counts.csv
```

## Pre-trained SCVI Models

The pipeline expects pre-trained SCVI latent embeddings. If you have these from a previous run:
```
scvi_models/
├── GSM5226574_C51ctr/
│   ├── adata.h5ad
│   ├── latent.npy
│   └── model.pt
├── GSM5226575_C52ctr/
│   └── ...
...
```

If you don't have these, set `--train-scvi` flag (not recommended; requires ~48 hours on HPC).

## Data Not Included

⚠️ **Not available from GEO:**
- Raw FASTQ files
- Droplet-level metadata
- CellRanger outputs
- Cell type annotations

These limitations prevent full reproduction of QC steps (CellBender, Scrublet).

## Citation

If you use this data, please cite:
```
Melms, J.C., Biermann, J., Huang, H. et al. 
A molecular single-cell lung atlas of lethal COVID-19. 
Nature 595, 114–119 (2021). 
https://doi.org/10.1038/s41586-021-03569-1
```
