#!/usr/bin/bash

Rscript ./normalize-expr.R --dir_data tmp-data/expression \
                           --vec_id_data "c('campbell2017','mikkelsen2019','moffitt2018','mousebrain')" \
                           --regex_raw_data umi \
                           --regex_metadata metadata \
                           --vec_metadata_cell_id_col "c('cell_id','cell_id','cell_id','cell_id')" \
                           --filename_out SCT_seuratObj \
                           --regressOut_percent_mt TRUE \
                           --dir_out tmp-data/expression

echo "bash script done"
