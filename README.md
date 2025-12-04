# COVID-19 snRNA-seq Reproduction Pipeline

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41586--021--03569--1-blue)](https://doi.org/10.1038/s41586-021-03569-1)

**Reproducibility project for:**  
*Melms et al., Nature 2021*  
"A molecular single-cell lung atlas of lethal COVID-19"

---

## 📋 Overview

This repository contains a computational pipeline to reproduce key analyses from the landmark COVID-19 single-nucleus RNA-seq study by Melms and colleagues. Using publicly available data from GEO (GSE171524), we implement integration, clustering, cell-type annotation, and differential expression analysis using **Python/Scanpy/scvi-tools**.

### Key Findings Reproduced:
✅ Integrated 81,000 high-quality nuclei across 27 samples (20 COVID, 7 control)  
✅ Identified broad cell populations: Fibroblast-ECM, Endothelial, Macrophage-like  
✅ Detected fibroblast subpopulations (ECM-high, Alveolar)  
✅ Differential expression reveals ECM remodeling and inflammatory signatures in COVID-19

**Limitations:** Due to limited gene coverage (161 HVGs vs. ~30K in original), fine-grained cell types (DATPs, CTHRC1+ fibroblasts, T cell subsets) could not be fully resolved.

---

## 🚀 Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/VedKulkarni01/covid19-snrnaseq-reproduction.git
cd covid19-snrnaseq-reproduction
```

### 2. Set Up Environment
```bash
# Option A: Conda (recommended)
conda env create -f environment.yml
conda activate covid19-snrnaseq

# Option B: pip
pip install -r requirements.txt
```

### 3. Download Data
```bash
cd data
bash download_data.sh
cd ..
```
> **Note:** Downloads ~10 GB of raw count matrices. See [`data/README.md`](data/README.md) for manual download instructions.

### 4. Run Pipeline

#### Option A: Jupyter Notebook (interactive)
```bash
jupyter lab notebooks/Project_sc_rna_seq.ipynb
```

#### Option B: Command-line script
```bash
# Dry run (check pipeline steps)
python code/run_pipeline.py --dry-run

# Full analysis
python code/run_pipeline.py \
  --data-dir data/ \
  --scvi-dir scvi_models/ \
  --output-dir output/
```

---

## 📊 Results

### UMAP Visualization
![UMAP colored by cell type](figures/final/umap_combined.png)

### Cell Type Distribution
| Cell Type | Control | COVID-19 | Total |
|-----------|---------|----------|-------|
| Fibroblast-ECM | 6,804 | 22,725 | 29,529 |
| Endothelial | 17,258 | 22,004 | 39,262 |
| Macrophage-like | 2,969 | 9,226 | 12,195 |

### Differential Expression
- **Fibroblasts:** COL1A1, COL3A1, FN1 upregulated (ECM remodeling)
- **Endothelial:** VWF pathway activation
- **Macrophages:** Inflammatory signatures consistent with IL-1β-driven pathology

Full results: [`output/DE_results/`](output/DE_results/)

---

## 🗂️ Repository Structure

```
covid19-snrnaseq-reproduction/
├── code/                    # Analysis scripts
│   └── run_pipeline.py      # Main pipeline
├── notebooks/               # Jupyter notebooks
│   └── Project_sc_rna_seq.ipynb
├── data/                    # Raw data (not tracked)
│   ├── README.md            # Download instructions
│   └── download_data.sh     # Automated download
├── scvi_models/             # Pre-trained models (optional)
├── output/                  # Generated results
│   ├── figures/
│   └── DE_results/
├── docs/                    # Documentation & reports
│   ├── individual_writeup.pdf
│   └── presentation.pdf
├── environment.yml          # Conda environment
├── requirements.txt         # pip requirements
└── README.md               # This file
```

---

## 🔬 Methods

### Integration & Batch Correction
- **scVI (single-cell Variational Inference)** for latent space embedding
- 10-dimensional latent representation across 27 samples
- UMAP for visualization, Leiden for clustering

### Cell Type Annotation
- Marker-based scoring (COL1A1, VWF, MRC1, etc.)
- Limited to broad cell types due to HVG constraints

### Differential Expression
- Wilcoxon rank-sum test on raw counts
- COVID vs. Control comparison per cell type
- Volcano plots & heatmaps for visualization

---

## 📦 Requirements

**Core Dependencies:**
- Python 3.11
- scanpy 1.11.5
- scvi-tools 1.4.0
- PyTorch 2.9.0 (CUDA 12.8)
- anndata 0.12.6
- pandas 2.3.3
- numpy 2.3.3

See [`environment.yml`](environment.yml) for complete list.

**Computational Requirements:**
- **RAM:** 32 GB minimum (64 GB recommended)
- **Storage:** ~15 GB for data + models
- **GPU:** Optional (speeds up SCVI training)
- **Estimated runtime:** 2-4 hours (without SCVI retraining)

---

## 🎓 Educational Use

This project was developed as part of **BINF 6430: Single-Cell Genomics** at Northeastern University. It demonstrates:
- Reproducible single-cell RNA-seq workflows
- Integration of public datasets
- Best practices in computational biology (version control, environment management, documentation)

**Instructor:** [Your Instructor's Name]  
**Semester:** Fall 2024

---

## 📖 Citation

If you use this code or build upon this work, please cite:

**Original Study:**
```bibtex
@article{melms2021molecular,
  title={A molecular single-cell lung atlas of lethal COVID-19},
  author={Melms, Johannes C and Biermann, Julia and Huang, Haoyang and others},
  journal={Nature},
  volume={595},
  number={7865},
  pages={114--119},
  year={2021},
  publisher={Nature Publishing Group},
  doi={10.1038/s41586-021-03569-1}
}
```

**This Reproduction:**
```bibtex
@misc{kulkarni2024covid19,
  author={Kulkarni, Vedant},
  title={COVID-19 snRNA-seq Reproduction Pipeline},
  year={2024},
  publisher={GitHub},
  url={https://github.com/VedKulkarni01/covid19-snrnaseq-reproduction}
}
```

---

## 🐛 Issues & Contributions

Found a bug? Have suggestions? Please [open an issue](https://github.com/VedKulkarni01/covid19-snrnaseq-reproduction/issues).

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for guidelines on contributing.

---

## 📄 License

This project is licensed under the **MIT License** - see [`LICENSE`](LICENSE) file for details.

**Data License:** Original data from Melms et al. is subject to GEO data use policies.

---

## 🙏 Acknowledgments

- **Melms et al.** for publishing the original study and making data publicly available
- **scvi-tools developers** for the integration framework
- **Scanpy team** for single-cell analysis tools
- **Northeastern University** for HPC resources (Discovery/Explorer clusters)

---

## 📞 Contact

**Vedant Kulkarni**  
GitHub: [@VedKulkarni01](https://github.com/VedKulkarni01)  
Email: [your.email@northeastern.edu]

---

**Last Updated:** December 2024
