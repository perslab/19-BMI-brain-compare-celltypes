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

## Session info

pander(sessionInfo())

**R version 3.5.3 (2019-03-11)**

**Platform:** x86_64-pc-linux-gnu (64-bit)

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:**
_parallel_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:**
_pander(v.0.6.3)_, _devtools(v.2.2.1)_, _usethis(v.1.5.0)_, _magrittr(v.1.5)_, _loomR(v.0.2.0)_, _itertools(v.0.1-3)_, _iterators(v.1.0.12)_, _R6(v.2.4.1)_, _Seurat(v.2.3.4)_, _cowplot(v.1.0.0)_, _ggplot2(v.3.2.1)_, _hdf5r(v.1.0.1)_, _Matrix(v.1.2-17)_, _here(v.0.1)_, _data.table(v.1.12.6)_, _dplyr(v.0.8.3)_ and _optparse(v.1.6.4)_

**loaded via a namespace (and not attached):**
_Rtsne(v.0.15)_, _colorspace(v.1.4-1)_, _ellipsis(v.0.3.0)_, _class(v.7.3-15)_, _modeltools(v.0.2-22)_, _ggridges(v.0.5.1)_, _mclust(v.5.4.5)_, _rprojroot(v.1.3-2)_, _htmlTable(v.1.13.2)_, _fs(v.1.3.1)_, _base64enc(v.0.1-3)_, _rstudioapi(v.0.10)_, _proxy(v.0.4-23)_, _npsurv(v.0.4-0)_, _remotes(v.2.1.0)_, _getopt(v.1.20.3)_, _flexmix(v.2.3-15)_, _bit64(v.0.9-7)_, _codetools(v.0.2-16)_, _splines(v.3.5.3)_, _R.methodsS3(v.1.7.1)_, _lsei(v.1.2-0)_, _robustbase(v.0.93-5)_, _knitr(v.1.26)_, _pkgload(v.1.0.2)_, _zeallot(v.0.1.0)_, _Formula(v.1.2-3)_, _jsonlite(v.1.6)_, _ica(v.1.0-2)_, _cluster(v.2.1.0)_, _kernlab(v.0.9-29)_, _png(v.0.1-7)_, _R.oo(v.1.23.0)_, _compiler(v.3.5.3)_, _httr(v.1.4.1)_, _backports(v.1.1.5)_, _assertthat(v.0.2.1)_, _lazyeval(v.0.2.2)_, _cli(v.1.1.0)_, _prettyunits(v.1.0.2)_, _lars(v.1.2)_, _acepack(v.1.4.1)_, _htmltools(v.0.4.0)_, _tools(v.3.5.3)_, _igraph(v.1.2.4.1)_, _gtable(v.0.3.0)_, _glue(v.1.3.1)_, _RANN(v.2.6.1)_, _reshape2(v.1.4.3)_, _Rcpp(v.1.0.3)_, _vctrs(v.0.2.0)_, _gdata(v.2.18.0)_, _ape(v.5.3)_, _nlme(v.3.1-142)_, _fpc(v.2.2-3)_, _gbRd(v.0.4-11)_, _lmtest(v.0.9-37)_, _xfun(v.0.11)_, _stringr(v.1.4.0)_, _ps(v.1.3.0)_, _testthat(v.2.3.0)_, _lifecycle(v.0.1.0)_, _irlba(v.2.3.3)_, _gtools(v.3.8.1)_, _DEoptimR(v.1.0-8)_, _MASS(v.7.3-51.4)_, _zoo(v.1.8-6)_, _scales(v.1.1.0)_, _doSNOW(v.1.0.18)_, _RColorBrewer(v.1.1-2)_, _memoise(v.1.1.0)_, _reticulate(v.1.13)_, _pbapply(v.1.4-2)_, _gridExtra(v.2.3)_, _rpart(v.4.1-15)_, _segmented(v.1.0-0)_, _latticeExtra(v.0.6-28)_, _stringi(v.1.4.3)_, _desc(v.1.2.0)_, _foreach(v.1.4.7)_, _checkmate(v.1.9.4)_, _caTools(v.1.17.1.2)_, _pkgbuild(v.1.0.6)_, _bibtex(v.0.4.2)_, _Rdpack(v.0.11-0)_, _SDMTools(v.1.1-221.1)_, _rlang(v.0.4.1)_, _pkgconfig(v.2.0.3)_, _dtw(v.1.21-3)_, _prabclus(v.2.3-1)_, _bitops(v.1.0-6)_, _lattice(v.0.20-38)_, _ROCR(v.1.0-7)_, _purrr(v.0.3.3)_, _htmlwidgets(v.1.5.1)_, _processx(v.3.4.1)_, _bit(v.1.1-14)_, _tidyselect(v.0.2.5)_, _plyr(v.1.8.4)_, _snow(v.0.4-3)_, _gplots(v.3.0.1.1)_, _Hmisc(v.4.2-0)_, _pillar(v.1.4.2)_, _foreign(v.0.8-72)_, _withr(v.2.1.2)_, _fitdistrplus(v.1.0-14)_, _mixtools(v.1.1.0)_, _survival(v.3.1-7)_, _nnet(v.7.3-12)_, _tsne(v.0.1-3)_, _tibble(v.2.1.3)_, _crayon(v.1.3.4)_, _KernSmooth(v.2.23-16)_, _grid(v.3.5.3)_, _callr(v.3.3.2)_, _metap(v.1.1)_, _digest(v.0.6.22)_, _diptest(v.0.75-7)_, _tidyr(v.1.0.0)_, _R.utils(v.2.9.0)_, _stats4(v.3.5.3)_, _munsell(v.0.5.0)_ and _sessioninfo(v.1.1.1)_>
