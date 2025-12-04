#!/usr/bin/env python3
"""
COVID-19 snRNA-seq Analysis Pipeline
Reproduction of Melms et al., Nature 2021

This script reproduces the main analysis steps using pre-computed SCVI latent embeddings.
"""

import os
import sys
import argparse
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run COVID-19 snRNA-seq reproduction pipeline'
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default='/scratch/roy.suc/final_project/data',
        help='Directory containing raw CSV files'
    )
    parser.add_argument(
        '--scvi-dir',
        type=str,
        default='/scratch/roy.suc/final_project/scvi_models',
        help='Directory containing pre-trained SCVI models'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./output',
        help='Output directory for results'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print pipeline steps without executing heavy computation'
    )
    parser.add_argument(
        '--skip-de',
        action='store_true',
        help='Skip differential expression analysis'
    )
    parser.add_argument(
        '--n-top-genes',
        type=int,
        default=2000,
        help='Number of highly variable genes (default: 2000)'
    )
    parser.add_argument(
        '--resolution',
        type=float,
        default=0.5,
        help='Leiden clustering resolution (default: 0.5)'
    )
    
    return parser.parse_args()


def load_integrated_data(scvi_dir, dry_run=False):
    """Load pre-integrated AnnData with SCVI latent embeddings."""
    print("=" * 60)
    print("STEP 1: Loading Pre-integrated Data")
    print("=" * 60)
    
    if dry_run:
        print("[DRY RUN] Would load SCVI outputs from:", scvi_dir)
        return None
    
    sample_dirs = sorted([
        d for d in os.listdir(scvi_dir)
        if os.path.isdir(os.path.join(scvi_dir, d))
    ])
    
    print(f"Found {len(sample_dirs)} SCVI sample folders")
    
    adatas = []
    for sample in sample_dirs:
        folder = os.path.join(scvi_dir, sample)
        adata_path = os.path.join(folder, "adata.h5ad")
        latent_path = os.path.join(folder, "latent.npy")
        
        if not os.path.exists(adata_path) or not os.path.exists(latent_path):
            print(f"⚠️  Skipping {sample} (missing files)")
            continue
        
        a = sc.read_h5ad(adata_path)
        z = np.load(latent_path)
        
        if a.n_obs != z.shape[0]:
            print(f"⚠️  Skipping {sample} (dimension mismatch)")
            continue
        
        a.obsm["X_scVI"] = z
        a.obs["sample"] = sample
        adatas.append(a)
    
    # Concatenate
    adata = ad.concat(
        adatas,
        join="inner",
        label="batch",
        keys=[a.obs["sample"][0] for a in adatas],
        index_unique=None
    )
    
    print(f"\n✓ Integrated: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"✓ Latent dimensions: {adata.obsm['X_scVI'].shape}")
    
    return adata


def quality_control(adata, dry_run=False):
    """Apply quality control filters."""
    print("\n" + "=" * 60)
    print("STEP 2: Quality Control")
    print("=" * 60)
    
    if dry_run:
        print("[DRY RUN] Would apply QC filters:")
        print("  - Mitochondrial RNA < 20%")
        print("  - Ribosomal RNA < 2%")
        print("  - Gene count 98th percentile filter")
        return adata
    
    # Calculate metrics
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT_"))
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPL", "RPS"))
    
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # Filter
    upper_ngenes = np.quantile(adata.obs["n_genes_by_counts"].values, 0.98)
    n_before = adata.n_obs
    
    adata = adata[adata.obs["n_genes_by_counts"] < upper_ngenes, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()
    adata = adata[adata.obs["pct_counts_ribo"] < 2, :].copy()
    
    print(f"✓ Cells before QC: {n_before}")
    print(f"✓ Cells after QC: {adata.n_obs}")
    print(f"✓ Removed: {n_before - adata.n_obs} cells")
    
    return adata


def clustering_umap(adata, resolution=0.5, dry_run=False):
    """Compute neighbors, UMAP, and Leiden clustering."""
    print("\n" + "=" * 60)
    print("STEP 3: Clustering & UMAP")
    print("=" * 60)
    
    if dry_run:
        print(f"[DRY RUN] Would compute:")
        print(f"  - Neighbors (k=15) on SCVI latent")
        print(f"  - UMAP embedding")
        print(f"  - Leiden clustering (resolution={resolution})")
        return adata
    
    # Neighbors on SCVI latent
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Leiden
    sc.tl.leiden(adata, resolution=resolution)
    
    print(f"✓ UMAP computed: {adata.obsm['X_umap'].shape}")
    print(f"✓ Leiden clusters: {len(adata.obs['leiden'].unique())}")
    
    return adata


def annotate_cell_types(adata, dry_run=False):
    """Annotate broad cell types using marker scoring."""
    print("\n" + "=" * 60)
    print("STEP 4: Cell Type Annotation")
    print("=" * 60)
    
    if dry_run:
        print("[DRY RUN] Would score markers for:")
        print("  - Fibroblast_ECM")
        print("  - Endothelial")
        print("  - Macrophage_like")
        return adata
    
    marker_dict = {
        "Fibroblast_ECM": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A3", "FN1", "DCN", "LUM"],
        "Endothelial": ["VWF", "PECAM1", "KDR"],
        "Macrophage_like": ["MRC1", "PPARG", "MARCO"]
    }
    
    gene_symbols = adata.var_names.str.upper()
    
    for celltype, genes in marker_dict.items():
        genes_upper = [g.upper() for g in genes]
        genes_present = list(set(genes_upper) & set(gene_symbols))
        
        if len(genes_present) == 0:
            print(f"⚠️  No markers found for {celltype}")
            continue
        
        print(f"✓ {celltype}: {len(genes_present)} markers")
        
        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=f"{celltype}_score",
            use_raw=False
        )
    
    # Assign best cell type
    score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
    if len(score_cols) > 0:
        score_mat = adata.obs[score_cols].to_numpy()
        best_idx = np.nanargmax(score_mat, axis=1)
        best_labels = np.array([c.replace("_score", "") for c in score_cols])[best_idx]
        adata.obs["broad_cell_type"] = best_labels
        
        print(f"\n✓ Cell type distribution:")
        print(adata.obs["broad_cell_type"].value_counts())
    
    return adata


def assign_condition(adata, dry_run=False):
    """Assign COVID vs Control condition."""
    print("\n" + "=" * 60)
    print("STEP 5: Condition Assignment")
    print("=" * 60)
    
    if dry_run:
        print("[DRY RUN] Would assign condition based on sample names")
        return adata
    
    def label_condition(sample):
        s = str(sample).lower()
        if "ctr" in s:
            return "Control"
        elif "cov" in s:
            return "COVID"
        else:
            return "Unknown"
    
    adata.obs["condition"] = adata.obs["sample"].apply(label_condition)
    
    print(f"✓ Condition distribution:")
    print(adata.obs["condition"].value_counts())
    
    return adata


def save_results(adata, output_dir, dry_run=False):
    """Save processed AnnData and key results."""
    print("\n" + "=" * 60)
    print("STEP 6: Saving Results")
    print("=" * 60)
    
    if dry_run:
        print(f"[DRY RUN] Would save to: {output_dir}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save AnnData
    adata_path = os.path.join(output_dir, "adata_integrated.h5ad")
    adata.write(adata_path)
    print(f"✓ Saved: {adata_path}")
    
    # Save metadata
    meta_path = os.path.join(output_dir, "cell_metadata.csv")
    adata.obs.to_csv(meta_path)
    print(f"✓ Saved: {meta_path}")
    
    # Save UMAP coordinates
    if "X_umap" in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"],
            index=adata.obs_names,
            columns=["UMAP_1", "UMAP_2"]
        )
        umap_path = os.path.join(output_dir, "umap_coordinates.csv")
        umap_df.to_csv(umap_path)
        print(f"✓ Saved: {umap_path}")
    
    print(f"\n✓ All results saved to: {output_dir}")


def main():
    """Main pipeline execution."""
    args = parse_args()
    
    print("\n" + "=" * 60)
    print("COVID-19 snRNA-seq Reproduction Pipeline")
    print("Melms et al., Nature 2021")
    print("=" * 60)
    print(f"Data directory: {args.data_dir}")
    print(f"SCVI directory: {args.scvi_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Dry run: {args.dry_run}")
    print("=" * 60)
    
    # Pipeline steps
    adata = load_integrated_data(args.scvi_dir, dry_run=args.dry_run)
    
    if not args.dry_run:
        adata = quality_control(adata, dry_run=args.dry_run)
        adata = clustering_umap(adata, resolution=args.resolution, dry_run=args.dry_run)
        adata = annotate_cell_types(adata, dry_run=args.dry_run)
        adata = assign_condition(adata, dry_run=args.dry_run)
        save_results(adata, args.output_dir, dry_run=args.dry_run)
    
    print("\n" + "=" * 60)
    print("✓ Pipeline Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
