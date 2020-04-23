#!/usr/bin/bash

#### GWAS (Fig 3A) ####

# BMI (BMI_UKBB_Loh2018)
# Educational Attainment (EA3_Lee2018)
# Intelligence (INTELLIGENCE_Savage2018)
# Schizophrenia (SCZ_Pardinas2018)
# Insomnia (INSOMNIA_Jansen2018)
# Multiple Schlerosis (MS_Patsopoulos2011) #TODO check 
# Rheumatoid Arthritis (RA_Okada2014) #TODO check 
# Low Density Lipoprotein (LIPIDS_LDL_Teslovich2010) #TODO check 
# Waist-hip Ratio (adj. BMI) (WHRadjBMI_UKBB_Loh2018) 
# Height (HEIGHT_UKBB_Loh2018) 

  Rscript ./correlate_celltypes.R --path_CELLECT_results /projects/timshel/sc-genetics/timshel-bmicelltypes2019/results/cellect_ldsc/prioritization.csv \
                                 --dir_geneset_results /projects/jonatan/pub-perslab/19-BMI-brain-genesettests/output/ \
                               --id_datExpr mousebrain \
                               --geneset_results_regex 200422 \
                               --vec_geneset_name 'c("protein_monogenic_extreme","syndromic")' \
                               --vec_GWAS 'c("BMI_UKBB_Loh2018")' \
                               --method pearson \
                               --output_label cor_CELLECT_geneset_MNS \
                               --append_results F

# vec_GWAS 'c("BMI_UKBB_Loh2018", "EA3_Lee2018", "INTELLIGENCE_Savage2018", "SCZ_Pardinas2018",  "INSOMNIA_Jansen2018", "MS_Patsopoulos2011","RA_Okada2014", "LIPIDS_LDL_Teslovich2010","WHRadjBMI_UKBB_Loh2018","HEIGHT_UKBB_Loh2018")' \
echo "bash script done"
