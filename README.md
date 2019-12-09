# Compare celltypes for BMI brain

Scripts to identify similar celltypes by transferring labels using Seurat 3, and to compute correlations between common- and rare-variant enrichment scores across celltypes.

## Quick start
1. Clone the repository: `git clone --recurse-submodules https://github.com/perslab/19-BMI-brain-compare-celltypes.git`
2. Convert the Zeisel et al (2018) 'Mouse Nervous System' dataset from .loom to .csv format: edit paths in `mousebrain-preprocessing-for-compare-celltypes.R` and call `Rscript mousebrain-preprocessing-for-compare-celltypes.R`
3. Use Seurat label transfer workflow to compare celltypes:
    1. CELLECT prioritized hypothalamus celltypes with each other: edit paths in `run_seurat_compare_hypo_CELLECT_hypo_CELLECT.sh` and call `bash run_seurat_compare_hypo_CELLECT_hypo_CELLECT.sh`
    2. CELLECT prioritized Mouse Nervous System celltypes with CELLECT prioritized hypothalamus celltypes: edit paths in `run_seurat_compare_MNS_CELLECT_hypo_CELLECT.sh` and call `run_seurat_compare_MNS_CELLECT_hypo_CELLECT.sh`
    3. Rare/mendelian variant Mouse Nervous System celltypes with CELLECT prioritized hypothalamus celltypes: edit paths in `run_seurat_compare_MNS_rarevariant_hypo_CELLECT.sh` and call `run_seurat_compare_MNS_rarevariant_hypo_CELLECT.sh`
4. Analyse the correlations between rare/mendelian BMI variants and CELLECT scores for 10 GWAS phenotypes: Update the paths in `run_correlate_celltypes.sh` and call `bash run_correlate_celltypes.sh`


