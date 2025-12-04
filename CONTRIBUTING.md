# Contributing Guide

Thank you for your interest in reproducing or building upon this COVID-19 snRNA-seq analysis!

## 🎯 Goals of This Repository

This project aims to:
1. **Reproduce** key findings from Melms et al. (Nature 2021)
2. **Demonstrate** best practices for computational reproducibility
3. **Provide** educational materials for single-cell RNA-seq analysis
4. **Enable** further exploration of COVID-19 lung pathology

---

## 🚀 Getting Started

### Prerequisites
- Python 3.11 (or 3.10+)
- 32 GB RAM minimum (64 GB recommended)
- ~15 GB disk space
- Optional: NVIDIA GPU with CUDA 12.8+ (for SCVI training)

### Setup Environment

**Option 1: Conda (Recommended)**
```bash
# Clone repository
git clone https://github.com/VedKulkarni01/covid19-snrnaseq-reproduction.git
cd covid19-snrnaseq-reproduction

# Create environment
conda env create -f environment.yml
conda activate covid19-snrnaseq

# Verify installation
python -c "import scanpy as sc, scvi; print(f'scanpy: {sc.__version__}, scvi: {scvi.__version__}')"
```

**Option 2: pip + virtualenv**
```bash
# Create virtual environment
python3.11 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "import scanpy as sc, scvi; print('OK')"
```

### Download Data
```bash
cd data
bash download_data.sh
cd ..
```

> **Note:** This downloads ~10 GB. See `data/README.md` for manual download.

---

## 🔬 Running the Analysis

### Jupyter Notebook (Interactive)
```bash
jupyter lab notebooks/Project_sc_rna_seq.ipynb
```

**Before running:**
1. Update paths in first cell (if not using default structure)
2. Set `RUN_SCVI_FROM_SCRATCH = False` (unless retraining models)
3. Run cells sequentially

### Command-Line Script
```bash
# Dry run (no computation, just checks)
python code/run_pipeline.py --dry-run

# Full analysis (uses pre-computed SCVI latents)
python code/run_pipeline.py \
  --data-dir data/ \
  --scvi-dir scvi_models/ \
  --output-dir output/

# Custom parameters
python code/run_pipeline.py \
  --resolution 0.8 \
  --n-top-genes 3000 \
  --output-dir ./my_results
```

### On HPC (Explorer/Discovery)
```bash
# Load modules
module load python/3.11
module load cuda/12.8

# Activate environment
source /path/to/venv/bin/activate

# Run pipeline
python code/run_pipeline.py \
  --data-dir /scratch/$USER/data \
  --scvi-dir /scratch/$USER/scvi_models \
  --output-dir /scratch/$USER/results
```

---

## 🛠️ Making Changes

### Code Style
- Follow **PEP 8** guidelines
- Use **type hints** for function arguments
- Add **docstrings** for all functions
- Keep functions **< 50 lines** when possible

**Example:**
```python
def calculate_qc_metrics(adata: ad.AnnData, mt_threshold: float = 20.0) -> ad.AnnData:
    """
    Calculate and filter QC metrics for single-cell data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    mt_threshold : float, default 20.0
        Mitochondrial percentage threshold
        
    Returns
    -------
    AnnData
        Filtered AnnData object
    """
    # Implementation
    return adata
```

### Testing Changes
```bash
# Run dry-run to check pipeline
python code/run_pipeline.py --dry-run

# Lint code
ruff check code/
# or
flake8 code/

# Format code
black code/
```

### Committing Changes
```bash
# Use conventional commit messages
git commit -m "feat: add new clustering method"
git commit -m "fix: correct UMAP parameter typo"
git commit -m "docs: update README quick start"
git commit -m "chore: update dependencies"
```

**Commit message format:**
- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation only
- `style:` Code formatting (no logic change)
- `refactor:` Code restructuring
- `test:` Adding tests
- `chore:` Maintenance tasks

---

## 📝 Submitting Changes

### Pull Request Process

1. **Fork** the repository
2. **Create branch** from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Make changes** and commit
4. **Push** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
5. **Open Pull Request** on GitHub with description:
   - What changed
   - Why the change is needed
   - Testing performed
   - Any breaking changes

### PR Checklist
- [ ] Code follows style guidelines
- [ ] Docstrings added for new functions
- [ ] `--dry-run` passes without errors
- [ ] No hardcoded absolute paths
- [ ] No large files added (check `.gitignore`)
- [ ] README updated if needed

---

## 🐛 Reporting Issues

### Bug Reports
Use the [issue tracker](https://github.com/VedKulkarni01/covid19-snrnaseq-reproduction/issues) with:
- **Title:** Short, descriptive summary
- **Environment:** OS, Python version, conda/pip
- **Steps to reproduce:** Exact commands run
- **Expected behavior:** What should happen
- **Actual behavior:** What actually happened
- **Error logs:** Full traceback if applicable

**Template:**
```markdown
**Environment:**
- OS: Ubuntu 22.04
- Python: 3.11.5
- scanpy: 1.11.5
- scvi-tools: 1.4.0

**Steps:**
1. Run `python code/run_pipeline.py`
2. Error occurs at UMAP step

**Error:**
```
KeyError: 'X_scVI' not found in adata.obsm
```

**Expected:** UMAP should compute on latent embedding
```

### Feature Requests
Describe:
- **Use case:** What problem does this solve?
- **Proposed solution:** How would you implement it?
- **Alternatives:** Other approaches considered

---

## 📚 Additional Resources

### Learning Materials
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [Single-cell Best Practices](https://www.sc-best-practices.org/)

### Related Projects
- [Scanpy GitHub](https://github.com/scverse/scanpy)
- [scvi-tools GitHub](https://github.com/scverse/scvi-tools)
- [Original Study Data](https://zenodo.org/record/4698034)

---

## 💬 Communication

- **GitHub Issues:** Bug reports, feature requests
- **GitHub Discussions:** General questions, show & tell
- **Email:** [your.email@northeastern.edu] for private inquiries

---

## 📄 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

## 🙏 Acknowledgments

Thank you for helping improve reproducibility in computational biology!

**Special thanks to:**
- Original authors (Melms et al.)
- scvi-tools development team
- Scanpy community

---

**Questions?** Open an issue or reach out via email!
