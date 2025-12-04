# Documentation

This directory contains project documentation, reports, and presentations.

## Files

### 📄 `individual_writeup.pdf`
**Individual Project Report**  
Full written report including:
- Introduction & background
- Methods & pipeline description
- Results & discussion
- Limitations & future directions
- Retrospective analysis

### 📊 `presentation.pdf`
**Final Project Presentation Slides**  
Presentation covering:
- Project motivation
- Dataset & methods
- Key findings
- Visualizations (UMAPs, volcano plots, heatmaps)
- Lessons learned

### 📋 `project_summary.pdf` (optional)
One-page summary with:
- 3-4 key findings
- Main figures
- Citation information

## Viewing Recommendations

**Presentation:**
- Best viewed in PDF reader with presentation mode
- Contains 23 slides with embedded figures
- Estimated presentation time: 15-20 minutes

**Write-up:**
- Comprehensive report (3 pages)
- Includes technical details and code snippets
- References the Jupyter notebook for reproducibility

## Export Commands

If you need to regenerate PDFs from source:

```bash
# From PowerPoint/Google Slides
# File → Export → PDF

# From Word/Google Docs
# File → Download → PDF Document

# From Jupyter Notebook (if creating summary)
jupyter nbconvert --to pdf notebooks/Project_sc_rna_seq.ipynb
```

## Additional Resources

- **Original Paper:** [Melms et al., Nature 2021](https://doi.org/10.1038/s41586-021-03569-1)
- **GEO Dataset:** [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524)
- **scvi-tools Documentation:** [https://docs.scvi-tools.org](https://docs.scvi-tools.org)
