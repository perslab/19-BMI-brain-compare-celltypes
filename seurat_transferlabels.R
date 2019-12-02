############### SYNOPSIS ###################
# Run Seurat transfer labels analysis

### OUTPUT:
#' csv table, example:
#' id_datExpr_test,cell_type_test,id_datExpr_ref,predicted.id,n_cells,median.prediction.score
#' moffitt2018.umi,e8_Cck_Ebf3,moffitt2018.umi,i9_Gaba,12,0.565629814207739

### REMARKS:
# ....

### REFERENCE:
# * 

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--dir_data", type="character", default=NULL,
              help = "character, relative path from project top dir to dir containing raw test input expression data, gene*cell table with gene names in the first column, [default %default]"), 
  make_option("--id_datExpr_test", type="character",
              help = "character, used to match raw expression matrix filename in dir_data, gene*cell table with gene names in the first column, also used as a label, e.g. 'campbell2017.umi'"), 
  # make_option("--path_seuratObj_test", type="character", default=NULL,
  #             help = "character, relative path to seuratObj with SCT or normalized data, variable genes and scale.data already exists, use this instead of path_datExpr_test to save time. Assumes that say and Idents should be set correctly, and that percent.mt has already been calculated. If path is provided to a file that does not exist, the script will save an object to the path provided for subsequent reuse. Not that, if conducting an analysis with a different subset of the same dataset, the seurat object should be saved under a new filename  [default %default]"), 
  # make_option("--id_metadata_test", type="character", default = NULL,
  #             help = "character, relative path from project top dir to test metadata table."), 
  # make_option("--metadata_test_cell_id_col", type="character", default = "cell_id",
  #             help = "character, name of test metadata column containing individual cell id, [default %default]"),
  make_option("--metadata_test_annot_col", type="character", default = "cell_type",
              help = "character, name of test metadata column containing cell annotation [default %default]"),
  make_option("--vec_metadata_test_annots_subset", type="character", default= NULL,
              help = "quoted vector of celltype ids to map [default $default]"), 
  # make_option("--dir_data", type="character", default = NULL,
  #             help = "character, relative path from project top dir to raw reference input expression data, gene*cell table with gene names in the first column, also used as a label, [default %default]"),
  make_option("--id_datExpr_ref", type="character",
              help = "character, used to match raw expression matrix filename in dir_data [default %default]"),
  # make_option("--path_seuratObj_ref", type="character", default=NULL,
  #             help = "character, relative path to seuratObj with SCT or normalized data, variable genes and scale.data already exists, use this instead of path_datExpr_ref to save time. Assumes that DefaultAssay and Idents should be set correctly, and that percent.mt has already been calculated. If path is provided to a file that does not exist, the script will save an object to the path provided.[default %default]"), 
  # make_option("--id_metadata_ref", type="character", default = NULL,
  #             help = "character, relative path from project top dir to reference metadata table."),
  # make_option("--metadata_ref_cell_id_col", type="character", default = "cell_id",
  #             help = "character, name of reference metadata column containing individual cell id, [default %default]"), 
  make_option("--metadata_ref_annot_col", type="character", default = "cell_type",
              help = "character, name of reference metadata column containing cell annotation [default %default]"), 
  # make_option("--vec_metadata_ref_annots_subset", type="character", default= NULL,
  #             help = "quoted vector of celltype ids to map to [default $default]"), 
  make_option("--nComp", type="integer", default= 50L,
              help = "integer, number of components to use in the reduced dimension [default $default]"), 
  # make_option("--regressOut_percent_mt", type="logical", default= T,
  #             help = "logical, regress out percent.mt? Requires gene names to be symbol, [default $default]"), 
  make_option("--output_label", type="character",
              help = "character, label for output files"),
  make_option("--dir_out", type="character",
              help = "Project directory to which to write outputs, relative to main project directory"),
  make_option("--append_results", type="logical", default=TRUE,
              help = "If file already exists, append results? [default %default] ")
)


opt <- parse_args(OptionParser(option_list=option_list))

dir_data <- opt$dir_data
id_datExpr_test <- opt$id_datExpr_test
#path_seuratObj_test <- opt$path_seuratObj_test
#id_metadata_test <- opt$id_metadata_test
#metadata_test_cell_id_col <- opt$metadata_test_cell_id_col
metadata_test_annot_col <- opt$metadata_test_annot_col
vec_metadata_test_annots_subset <- opt$vec_metadata_test_annots_subset
if (!is.null(opt$vec_metadata_test_annots_subset)) vec_metadata_test_annots_subset <- eval(parse(text=vec_metadata_test_annots_subset))
#dir_data <- opt$dir_data
id_datExpr_ref <- opt$id_datExpr_ref
#path_seuratObj_ref <- opt$path_seuratObj_ref
#id_metadata_ref <- opt$id_metadata_ref
#metadata_ref_cell_id_col <- opt$metadata_ref_cell_id_col
metadata_ref_annot_col <- opt$metadata_ref_annot_col
# vec_metadata_ref_annots_subset <- opt$vec_metadata_ref_annots_subset
# if (!is.null(opt$vec_metadata_ref_annots_subset)) vec_metadata_ref_annots_subset <- eval(parse(text=vec_metadata_ref_annots_subset))
#regressOut_percent_mt <- opt$regressOut_percent_mt
nComp <- opt$nComp
output_label <- opt$output_label
dir_out <- opt$dir_out
append_results <- opt$append_results

######################################################################
######################### CHECK INPUT VALIDITY #######################
######################################################################

# stopifnot(any(sapply(c(dir_data, path_seuratObj_test), function(x)!is.null(x))))
# stopifnot(any(sapply(c(dir_data, path_seuratObj_ref), function(x)!is.null(x))))

######################################################################
########################### PACKAGES #################################
######################################################################

message("Loading packages")

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Seurat"))

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F, 
        use="pairwise.complete.obs", 
        warn=1, 
        verbose=F,
        mc.cores=40 # for parallel computation
) 

######################################################################
############################ CONSTANTS ###############################
######################################################################

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed = 12345
set.seed(randomSeed)

######################################################################
########################## LOAD AND REFORMAT FILES ###################
######################################################################

# if (!is.null(path_seuratObj_test)) {
#   if (file.exists(here(path_seuratObj_test))) {
#     warning(" Seurat object detected for test data, ignoring raw data")
#     dir_data <- NULL
#   }
# }

# if (!is.null(path_seuratObj_ref)) {
#   if (file.exists(here(path_seuratObj_ref))) {
#     warning(" Seurat object detected for ref data, ignoring raw data")
#     dir_data <- NULL
#   }
# }
message("Loading data")

# if (!is.null(dir_data)) {
#   
#   seuratObj_test <- NULL
#   
  # path_datExpr_test <- dir(path=here(dir_data), pattern=id_datExpr_test, full.names = T)
  # 
  # dt_datExpr_test <- fread(path_datExpr_test)
  # mat_datExpr_test <- as.matrix(dt_datExpr_test[,-1])
  # rownames(mat_datExpr_test) <- dt_datExpr_test[,1][[1]]
  # rm(dt_datExpr_test)
  # 
  # path_metadata_test <- dir(path=here(dir_data), pattern=id_metadata_test, full.names = T)
  # dt_metadata_test <- fread(path_metadata_test)
  # dt_metadata_test <- dt_metadata_test[,c(..metadata_test_cell_id_col,..metadata_test_annot_col),]
  # 
  # mat_datExpr_test <- mat_datExpr_test[,match(dt_metadata_test[,..metadata_test_cell_id_col,][[1]],colnames(mat_datExpr_test))]

# } else {
#  mat_datExpr_test <- NULL


path_seuratObj_test <- dir(path=dir_data, pattern=paste0(id_datExpr_test, ".*seuratObj.RDS.gz$"), full.names=T)
seuratObj_test <- readRDS(here(path_seuratObj_test))
Idents(seuratObj_test) <- seuratObj_test@meta.data[[metadata_test_annot_col]]
#}

# if (!is.null(dir_data)) {
#   
#   seuratObj_ref <- NULL
#   
#   path_datExpr_ref <- dir(path=here(dir_data), pattern=id_datExpr_ref, full.names = T)
#   dt_datExpr_ref <- fread(path_datExpr_ref)
#   mat_datExpr_ref <- as.matrix(dt_datExpr_ref[,-1])
#   rownames(mat_datExpr_ref) <- dt_datExpr_ref[,1][[1]]
#   rm(dt_datExpr_ref)
#   
#   path_metadata_ref <- dir(path=here(dir_data), pattern=id_metadata_ref, full.names = T)
#   dt_metadata_ref <- fread(path_metadata_ref)
#   dt_metadata_ref <- dt_metadata_ref[,c(..metadata_ref_cell_id_col,..metadata_ref_annot_col),]
#   
#   mat_datExpr_ref <- mat_datExpr_ref[,match(dt_metadata_ref[,..metadata_ref_cell_id_col,][[1]],colnames(mat_datExpr_ref))]
# 
#   
# } else {
#  mat_datExpr_ref <- NULL
path_seuratObj_ref <- dir(path=dir_data, pattern=paste0(id_datExpr_ref, ".*seuratObj.RDS.gz$"), full.names=T)
seuratObj_ref<- readRDS(here(path_seuratObj_ref))
Idents(seuratObj_ref) <- seuratObj_ref@meta.data[[metadata_ref_annot_col]]
#}

######################################################################
##################### QUANTIFY CONFOUNDERS ###########################
######################################################################
# if (!is.null(mat_datExpr_test)) {
#   # percent mitochrondria RNA
#   vec_logic_mt.genes_test <- grepl(pattern = "^MT-|^Mt-|^mt-", rownames(mat_datExpr_test))
#   vec_pct.mt_test <- (colSums(mat_datExpr_test[vec_logic_mt.genes_test,])*100) / colSums(mat_datExpr_test)
#   names(vec_pct.mt_test) <- colnames(mat_datExpr_test)
#   if (sum(vec_pct.mt_test)==0) warning(paste0("No mitochrondrial genes detected in ", id_datExpr_test))
# }
# 
# if (!is.null(mat_datExpr_ref)) {
#   # percent mitochrondria RNA
#   vec_logic_mt.genes_ref <- grepl(pattern = "^MT-|^Mt-|^mt-", rownames(mat_datExpr_ref))
#   vec_pct.mt_ref <- (colSums(mat_datExpr_ref[vec_logic_mt.genes_ref,])*100) / colSums(mat_datExpr_ref)
#   names(vec_pct.mt_ref) <- colnames(mat_datExpr_ref)
#   if (sum(vec_pct.mt_ref)==0) warning(paste0("No mitochrondrial genes detected in ", id_datExpr_ref))
# }
# 
# #####################################################################
# ########################## MAP GENES TO ENSEMBL? #####################
# ######################################################################
# 
# vec_genes_test <- if (!is.null(mat_datExpr_test)) rownames(mat_datExpr_test) else rownames(seuratObj_test)
# 
# vec_genes_ref <- if(!is.null(mat_datExpr_ref)) rownames(mat_datExpr_ref) else rownames(seuratObj_ref)
# 
# if (!any(grepl("ENSMUS",vec_genes_test)) & 
#     any(grepl("ENSMUS",vec_genes_ref))) {
#   
#   # For mapping symbol to ensembl
#   mapping_mm_filepath = here("src/compare-celltypes/gene_mapping", "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
#   # Synonyms
#   mapping_mm_synonyms_filepath = here("src/compare-celltypes/gene_mapping", "Mus_musculus.gene_info_symbol2ensembl.gz")
#   
#   # Step 1: direct mapping
#   mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
#   mapping = data.frame(symbol=vec_genes_test, 
#                        ensembl = mapping_direct$ensembl_gene_id[ match(toupper(gsub("-|_|\\.",".",vec_genes_test)), toupper(gsub("-|_|\\.",".",mapping_direct$gene_name_optimal))) ])
#   
#   # Step 2: map remaining using synonyms
#   mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
#   mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(gsub("-|_|\\.",".",mapping$symbol[which(is.na(mapping$ensembl)) ])) , toupper(gsub("-|_|\\.",".",mapping_synonyms$symbol)))]
#   
#   rm(mapping_direct, mapping_synonyms)
#   
#   # Make a log of unmapped genes 
#   log_not_mapped_filepath = here(dir_out, paste0(output_label,"_log_genes_not_mapped_testdata_", flag_date, ".tab"))
#   log_duplicate_filepath = here(dir_out, paste0(output_label,"_log_duplicate_ensembl_id_genenames_testdata_", flag_date, ".tab"))
#   
#   df_not_mapped = mapping[is.na(mapping$ensembl),]
#   write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
#   
#   # Make a log of duplicate genes
#   idx_duplicate_genes <- duplicated(mapping$ensembl)
#   df_duplicate <- mapping[idx_duplicate_genes,]
#   write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
#   
#   if (!is.null(mat_datExpr_test))  {
#    
#     # Filter out unmapped and duplicate genes from Seurat object
#     mat_datExpr_test <- mat_datExpr_test[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
#   
#     # rename rows where mapping was successful to ensembl ID
#     rownames(mat_datExpr_test) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
#   } else {
#     # Filter out unmapped and duplicate genes
#     SetAssayData(object=seuratObj_test, new.data = GetAssayData(seuratObj_test)[!is.na(mapping$ensembl) & !idx_duplicate_genes,])
#     # rename rows where mapping was successful to ensembl ID
#     rownames(seuratObj_test) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
#   }
#   
# } else if (any(grepl("ENSMUS",vec_genes_test)) & 
#            !any(grepl("ENSMUS",vec_genes_ref))) {
#   
#   # For mapping symbol to ensembl
#   mapping_mm_filepath = here("src/compare-celltypes/gene_mapping",  "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
#   
#   # Synonyms
#   mapping_mm_synonyms_filepath = here("src/compare-celltypes/gene_mapping",  "Mus_musculus.gene_info_symbol2ensembl.gz")
#   
#   # Step 1: direct mapping
#   mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
#   mapping = data.frame(symbol=vec_genes_ref, 
#                        ensembl = mapping_direct$ensembl_gene_id[ match(toupper(gsub("-|_|.","\\.",vec_genes_ref)), toupper(gsub("-|_|\\.",".",mapping_direct$gene_name_optimal))) ])
#   
#   # Step 2: map remaining using synonyms
#   mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
#   mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(gsub("-|_|\\.",".",mapping$symbol[which(is.na(mapping$ensembl)) ])) , toupper(gsub("-|_|\\.",".",mapping_synonyms$symbol)))]
#   rm(mapping_direct, mapping_synonyms)
#  
#   # Make a log of unmapped genes 
#   log_not_mapped_filepath = here(dir_out, paste0(output_label, "_log_genes_not_mapped_refdata_", flag_date, ".tab"))
#   log_duplicate_filepath = here(dir_out, paste0(output_label, "log_duplicate_ensembl_id_genenames_refdata_", flag_date, ".tab"))
#   
#   df_not_mapped = mapping[is.na(mapping$ensembl),]
#   write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
#   
#   # Make a log of duplicate genes
#   idx_duplicate_genes <- duplicated(mapping$ensembl)
#   df_duplicate <- mapping[idx_duplicate_genes,]
#   write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)
#   
#   if (!is.null(mat_datExpr_ref))  {
#     
#     # Filter out unmapped and duplicate genes
#     mat_datExpr_ref <- mat_datExpr_ref[!is.na(mapping$ensembl) & !idx_duplicate_genes,]
#     # rename rows where mapping was successful to ensembl ID
#     rownames(mat_datExpr_ref) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes] 
#   } else {
#     # Filter out unmapped and duplicate genes
#     SetAssayData(object=seuratObj_ref, new.data = GetAssayData(seuratObj_ref)[!is.na(mapping$ensembl) & !idx_duplicate_genes,])
#     
#     # rename rows where mapping was successful to ensembl ID
#     rownames(seuratObj_ref) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes] 
#   }
# }


######################################################################
############### SUBSET TEST AND / REF DATA BY ANNOTATIONS? ###########
######################################################################
# 
# if (!is.null(vec_metadata_test_annots_subset)) {
# 
#   if (!is.null(dir_data)) {
#     message(paste0("Subsetting test data ", id_datExpr_test))
#     vec_logical_condition  <- quote(dt_metadata_test[,..metadata_test_annot_col,][[1]] %in% vec_metadata_test_annots_subset)
#     dt_metadata_test <- dt_metadata_test[eval(vec_logical_condition),]
#     mat_datExpr_test <- mat_datExpr_test[,colnames(mat_datExpr_test) %in% dt_metadata_test[,..metadata_test_cell_id_col][[1]]]
#     } else  {
#     seuratObj_test <- subset(x=seuratObj_test, idents = vec_metadata_test_annots_subset)
#   }
# }
# 
# if (!is.null(vec_metadata_ref_annots_subset)) {
# 
#   if (!is.null(dir_data)) {
#     message(paste0("Subsetting reference data ", id_datExpr_ref))
#     vec_logical_condition  <- quote(dt_metadata_ref[,..metadata_ref_annot_col,][[1]] %in% vec_metadata_ref_annots_subset)
#     dt_metadata_ref <- dt_metadata_ref[eval(vec_logical_condition),,]
#     mat_datExpr_ref <- mat_datExpr_ref[,colnames(mat_datExpr_ref) %in% dt_metadata_ref[,..metadata_ref_cell_id_col][[1]]]
#   } else {
#     seuratObj_ref <- subset(x=seuratObj_ref, idents = vec_metadata_ref_annots_subset)
#   }
# }

######################################################################
#################### CREATE SEURAT OBJECTS ###########################
######################################################################

# if (!is.null(dir_data)) {
#   
#   message(paste0("Creating Seurat object for test data ",id_datExpr_test))
#   
#   df_metadata_test <- data.frame(
#     row.names = dt_metadata_test[[metadata_test_cell_id_col]],
#     "cell_type" = dt_metadata_test[[metadata_test_annot_col]],
#     "percent.mt" = vec_pct.mt_test[match(dt_metadata_test[[metadata_test_cell_id_col]], names(vec_pct.mt_test))])
#   colnames(df_metadata_test)[colnames(df_metadata_test)=="cell_type"] <- metadata_test_annot_col
#   seuratObj_test <- CreateSeuratObject(counts = mat_datExpr_test,assay = "RNA",meta.data = df_metadata_test)
#   Idents(seuratObj_test) <- seuratObj_test@meta.data[[metadata_test_annot_col]]
# }
# 
# if (!is.null(dir_data)) {
#   
#   message(paste0("Creating Seurat object for reference data ",id_datExpr_ref))
#   
#   df_metadata_ref <-data.frame(
#     row.names = dt_metadata_ref[[metadata_ref_cell_id_col]],
#     "cell_type" = dt_metadata_ref[[metadata_ref_annot_col]],
#     "percent.mt" = vec_pct.mt_ref[match(dt_metadata_ref[[metadata_ref_cell_id_col]], names(vec_pct.mt_ref))])
#   colnames(df_metadata_ref)[colnames(df_metadata_ref)=="cell_type"] <- metadata_ref_annot_col
#   seuratObj_ref <- CreateSeuratObject(counts = mat_datExpr_ref,assay = "RNA",meta.data = df_metadata_ref)
#   Idents(seuratObj_ref) <- seuratObj_ref@meta.data[[metadata_ref_annot_col]]
# }
######################################################################
#################### SC-TRANSFORM NORMALIZE THE DATA #################
######################################################################

# if (!is.null(dir_data)) {
#   # test data
#   
#   message(paste0("SCT normalising test data ",id_datExpr_test))
#   
#   seuratObj_test <- SCTransform(object=seuratObj_test,
#                                 do.correct.umi = T,
#                                 vars.to.regress = if (regressOut_percent_mt & sum(seuratObj_test$percent.mt) > 0) c("percent.mt") else NULL,
#                                 do.scale = F,
#                                 do.center = T,
#                                 seed.use = randomSeed,
#                                 verbose=F)
# }
# 
# if (!is.null(dir_data)) {
#   
#   message(paste0("SCT normalising reference data ",id_datExpr_ref))
#   
#   # reference data
#   seuratObj_ref <- SCTransform(object=seuratObj_ref,
#                               do.correct.umi = T,
#                               vars.to.regress = if (regressOut_percent_mt & sum(seuratObj_ref$percent.mt) > 0) c("percent.mt") else NULL,
#                               do.scale = F,
#                               do.center = T,
#                               seed.use = randomSeed,
#                               verbose=F)
#                             
# 
# }


######################################################################
###################### FIND TRANSFER ANCHORS #########################
######################################################################
# project test dataset onto reference PCA subspace

anchors <- FindTransferAnchors(reference = seuratObj_ref, 
                               reduction="pcaproject",
                               query = seuratObj_test, 
                               npcs = if (!is.null(seuratObj_ref@reductions[["pca"]])) NULL else nComp,
                               dims = 1:nComp, 
                               verbose=T)

######################################################################
############# SAVE SEURAT OBJECTS IF NOT ALREADY DONE ################
######################################################################

# if (!is.null(path_seuratObj_test)) {
#   if (!file.exists(here(path_seuratObj_test))) {
#     message(paste0("Saving normalized test data as Seurat Object to ", here(path_seuratObj_test)))
#     file.out.test <- here(path_seuratObj_test)
#     #if (!is.null(vec_metadata_test_annots_subset)) file.out.test <- gsub("\\.RDS\\.gz$", "_subset.RDS.gz", file.out.test)
#     saveRDS(seuratObj_test, file = file.out.test, compress="gzip")
#   }
# }
# 
# if (!is.null(path_seuratObj_ref)) {
#   if (!file.exists(here(path_seuratObj_ref))) {
#     message(paste0("Saving normalized ref data as Seurat Object to ", here(path_seuratObj_ref)))
#     file.out.ref <- here(path_seuratObj_ref)
#     #if (!is.null(vec_metadata_ref_annots_subset)) file.out.ref <- gsub("\\.RDS\\.gz$", "_subset.RDS.gz", file.out.ref)
#     saveRDS(seuratObj_ref, file = file.out.ref, compress="gzip")
#   }
# }

######################################################################
######################## MAKE CELLTYPE PREDICTIONS #######################
######################################################################

predictions <- TransferData(anchorset = anchors, 
                            refdata = seuratObj_ref@meta.data[[metadata_ref_annot_col]],  
                            dims = 1:nComp, 
                            verbose=T) 

colnames(predictions) <- paste0(metadata_ref_annot_col, "_", colnames(predictions))

######################################################################
########################## prepare results ###########################
######################################################################
message("preparing outputs")

predictions$cell_id <- colnames(seuratObj_test)
  
if (!is.null(vec_metadata_test_annots_subset))  {
  seuratObj_test <- subset(x=seuratObj_test, idents=vec_metadata_test_annots_subset)
  predictions <- predictions[predictions$cell_id %in% colnames(seuratObj_test)]
}

predictions <- data.table("cell_type_test"=seuratObj_test@meta.data[[metadata_test_annot_col]],
                          predictions)


predictions_long <- reshape2::melt(predictions, 
                                   id.vars= "cell_id", 
                                   measure.vars = colnames(predictions)[!colnames(predictions) %in% c("cell_id","cell_type_test", paste0(metadata_ref_annot_col, "_predicted.id"))])

predictions_long$cell_type_test = seuratObj_test@meta.data[[metadata_test_annot_col]][match(predictions_long$cell_id, rownames(seuratObj_test@meta.data))]
predictions_long$id_datExpr_test = id_datExpr_test
predictions_long$id_datExpr_ref = id_datExpr_ref
predictions_long$predicted.id <- predictions[[paste0(metadata_ref_annot_col, "_predicted.id")]][match(predictions_long$cell_id, predictions$cell_id)]
predictions_long$prediction.score.max <- predictions[[paste0(metadata_ref_annot_col, "_prediction.score.max")]][match(predictions_long$cell_id, predictions$cell_id)]

######################################################################
########################## prepare summary data.table ################
######################################################################

predictions_long <- setDT(predictions_long)

predictions_summary <- predictions_long[,
                                      .(n_cells = .N,
                                        median.prediction.score = median(prediction.score.max)),
                                      by = .(id_datExpr_test,cell_type_test,id_datExpr_ref,predicted.id),]           

predictions_summary <- predictions_summary[order(cell_type_test, -n_cells)]

######################################################################
#############################  WRAP UP  ##############################
######################################################################
#message("Writing out results")
#file.out.results.full <- here(dir_out, paste0(output_label, "_cf_celltypes_full.csv"))
#fwrite(x = predictions_long, file =file.out.results.full,  append = append_results, nThread=24, verbose=T)

#R.utils::gzip(file.out.results.full, overwrite=TRUE) # gzip

file.out.results.summary<- here(dir_out, paste0(output_label, "_cf_celltypes_sum.csv"))
if (!file.exists(file.out.results.summary)) append_results <- F
fwrite(x = predictions_summary, file = file.out.results.summary,  append = append_results, nThread=24, verbose=T)
message(paste0("Seurat ", id_datExpr_ref, " label transfer to ",  id_datExpr_test, " done!"))
# if (file.exists(file.out.results.summary)){
#   dt_existing <- fread(file.out.results.summary)
#   dt_combined <- rbind.data.table(l=list(dt_existing,predictions_summary))
# } else {
#   dt_combined <- dt_existing <- predictions_summary # this is just a hack to make the next condition return F
# }

# if (append_results & all.equal(which(duplicated(dt_combined)),(nrow(dt_existing)+1):nrow(dt_combined))) {
#   stop(paste0("comparisons already exist in ", file.out.results.summmary, ", will not overwrite"))
# } else {
#   fwrite(x = predictions_summary, file = file.out.results.summary,  append = append_results, nThread=24, verbose=T)
#   #R.utils::gzip(file.out.results.summary, overwrite=TRUE) # gzip
#   message(paste0("Seurat ", id_datExpr_ref, " label transfer to ",  id_datExpr_test, " done!"))
# }
  


