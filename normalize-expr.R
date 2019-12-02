############### SYNOPSIS ###################
# lognormalise datasets as preprarat

### INPUT: 
# raw count expression data matrices

### OUTPUT: 
# lognormalized expression data matrices

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ================================ OPTPARSE ============================== #
# ======================================================================= #

library(optparse)

option_list <- list(
  make_option("--dir_data", type="character", default = "tmp-data/expression",
              help = "Character, relative path to directory containing raw expression gene * sample matrices, relative to top working directory, default [$default]"),
  make_option("--vec_id_data", type="character",
              help = "**Quoted** character vector of dataset strings to match e.g. ''('mikkelsen2019','moffitt2018')''"),
  make_option("--regex_raw_data", type="character", default = "umi",
              help = "Character, string to match raw data files [default %default]"),
  make_option("--regex_metadata", type="character", default = "metadata",
              help = "Character, string to match metadata files [default %default]"),
  make_option("--vec_metadata_cell_id_col", type="character",
              help = "**Quoted** character, name of test metadata column containing individual cell id, [default %default]"),
  make_option("--regressOut_percent_mt", type="logical", default= T,
              help = "logical, regress out percent.mt? Requires gene names to be symbol, [default $default]"), 
  make_option("--filename_out", type="character", default = "SCT_seuratObj",
              help = "Character, generic name of output files, [default %default]"),
  make_option("--dir_out", type="character", default = "tmp-data/expression",
              help = "Character, relative path to directory to which to write files, relative to top working directory, default [$default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir_data <-opt$dir_data
vec_id_data <- eval(parse(text=opt$vec_id_data))
regex_raw_data <- opt$regex_raw_data
regex_metadata <- opt$regex_metadata
vec_metadata_cell_id_col <- eval(parse(text=opt$vec_metadata_cell_id_col))
regressOut_percent_mt <- opt$regressOut_percent_mt
filename_out<- opt$filename_out
dir_out <- opt$dir_out

# ======================================================================= #
# ========================== Packages ================================ #
# ======================================================================= #

library(here)
library(Seurat)
library(Matrix)
library(magrittr)
library(data.table)

# ======================================================================= #
# ========================== Options and constants ============================ #
# ======================================================================= #

options(stringsAsFactors = F, 
        use="pairwise.complete.obs", 
        warn=1, 
        verbose=F,
        mc.cores=40 # for parallel computation
) 
 
randomSeed = 12345
set.seed(randomSeed)


flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

# ======================================================================= #
# ========== lognormalize expression matrices and write out ==============#
# ======================================================================= #

fun = function(i) {  
  
  id_data = vec_id_data[i]
  
  file_matches = dir(path=here(dir_data), pattern=id_data, full.names=T)
  
  path_datExpr = grep(pattern = regex_raw_data, x = file_matches, value=T)
  
  if (!length(path_datExpr)) stop(paste0(id_data, ": no files match", regex_raw_data))
  if (length(path_datExpr)>1) stop(paste0(id_data, ": multiple files match ", regex_raw_data))
  
  message(paste0("Loading ", id_data, " raw data from disk"))
  
  # count data 
  dt_datExpr <- fread(input = path_datExpr)
  mat_datExpr <- as.matrix(dt_datExpr[,-1])
  rownames(mat_datExpr) <- dt_datExpr[,1,][[1]] 
  rm(dt_datExpr)
  
  # metadata 
  path_metadata <- grep(pattern = regex_metadata, x=file_matches, value=T)
  
  if (!length(path_metadata)) stop(paste0(id_data, ": no files match", regex_metadata))
  if (length(path_metadata)>1) stop(paste0(id_data, ": multiple files match ", regex_metadata))

  message(paste0("Loading ", id_data, " metadata from disk"))

  dt_metadata <- fread(path_metadata)
  metadata_cell_id_col <- vec_metadata_cell_id_col[i]
  df_metadata <- data.frame(dt_metadata, row.names =  dt_metadata[[metadata_cell_id_col]])
  rm(dt_metadata)
  
  mat_datExpr <- mat_datExpr[,match(rownames(df_metadata),colnames(mat_datExpr))]
  
  #####################################################################
  ############################## Percent.mt ###########################
  #####################################################################
  if (regressOut_percent_mt) {
    # percent mitochrondria RNA
    vec_logic_mt.genes <- grepl(pattern = "^MT-|^Mt-|^mt-", rownames(mat_datExpr))
    vec_pct.mt <- (colSums(mat_datExpr[vec_logic_mt.genes,])*100) / colSums(mat_datExpr)
    names(vec_pct.mt) <- colnames(mat_datExpr)
    if (sum(vec_pct.mt)>0) {
      df_metadata$percent.mt <- vec_pct.mt 
    } else {
      warning(paste0("No mitochrondrial genes detected in ", id_datExpr))
    }
  } 
  #####################################################################
  ########################## MAP GENES TO ENSEMBL #####################
  ######################################################################

  vec_genes <-  rownames(mat_datExpr) 

  if (!any(grepl("ENSMUS",vec_genes))) {
  
    # For mapping symbol to ensembl
    mapping_mm_filepath = here("src/compare-celltypes/gene_mapping", "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
    # Synonyms
    mapping_mm_synonyms_filepath = here("src/compare-celltypes/gene_mapping", "Mus_musculus.gene_info_symbol2ensembl.gz")

    # Step 1: direct mapping
    mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T, stringsAsFactors = F)
    mapping = data.frame(symbol=vec_genes,
                         ensembl = mapping_direct$ensembl_gene_id[ match(toupper(gsub("-|_|\\.",".",vec_genes)), toupper(gsub("-|_|\\.",".",mapping_direct$gene_name_optimal))) ])

    # Step 2: map remaining using synonyms
    mapping_synonyms = read.delim(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T, stringsAsFactors = F)
    mapping$ensembl[ which(is.na(mapping$ensembl)) ] = mapping_synonyms$ensembl[ match( toupper(gsub("-|_|\\.",".",mapping$symbol[which(is.na(mapping$ensembl)) ])) , toupper(gsub("-|_|\\.",".",mapping_synonyms$symbol)))]

    rm(mapping_direct, mapping_synonyms)

    # Make a log of unmapped genes
    log_not_mapped_filepath = here(dir_out, paste0(id_data,"_log_genes_not_mappeddata_", flag_date, ".tab"))
    log_duplicate_filepath = here(dir_out, paste0(id_data,"_log_duplicate_ensembl_id_genenamesdata_", flag_date, ".tab"))

    df_not_mapped = mapping[is.na(mapping$ensembl),]
    write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)

    # Make a log of duplicate genes
    idx_duplicate_genes <- duplicated(mapping$ensembl)
    df_duplicate <- mapping[idx_duplicate_genes,]
    write.table(df_duplicate,log_duplicate_filepath,quote=F,sep="\t",row.names=F)

    # Filter out unmapped and duplicate genes from Seurat object
    mat_datExpr <- mat_datExpr[!is.na(mapping$ensembl) & !idx_duplicate_genes,]

    # rename rows where mapping was successful to ensembl ID
    rownames(mat_datExpr) <- mapping$ensembl[!is.na(mapping$ensembl) & !idx_duplicate_genes]
  } 
  
  #####################################################################
  ######################### CREATE SEURAT OBJECT ######################
  #####################################################################
  
  seuratObj <- CreateSeuratObject(counts = mat_datExpr, 
                                  project= id_data,
                                  meta.data = df_metadata, 
                                  assay = "RNA", 
                                  min.cells = 0, 
                                  min.features = 0)
  
  rm(mat_datExpr, df_metadata)
  
  seuratObj <- SCTransform(object=seuratObj,
                           do.correct.umi = T,
                           vars.to.regress = if (!is.null(seuratObj$percent.mt)) "percent.mt" else NULL, # TODO if (regressOut_percent_mt & sum(seuratObj_test$percent.mt) > 0) c("percent.mt") else NULL,
                           do.scale = F,
                           do.center = T,
                           seed.use = randomSeed,
                           verbose=F)#NormalizeData(object = seuratObj)
  
  # write out to file
  message(paste0("saving ", id_data, " to disk"))
  
  saveRDS(object=seuratObj, file = here(dir_out, paste0(id_data,".", filename_out, ".RDS.gz")), compress="gzip")
  # GetAssayData(object=seuratObj, slot="data") %>% as.data.table -> dt_data 
  # dt_data <- data.table("gene"=rownames(seuratObj), dt_data)
  # rm(seuratObj)
  # 
  # print(dt_data[0:3,0:3])
  # 
  # gsub(".*/", "", path_datExpr) %>% 
  #   gsub(pattern = file_regex, replacement = output_label, x= .) %>% 
  #   gsub("\\.gz","",.)-> filename
  # 
  # file.out.data <- here(dir_out, filename)
  # data.table::fwrite(dt_data, file=file.out.data,  # fwrite cannot write gziped files
  #                    nThread=24, verbose=T)  
  # R.utils::gzip(file.out.data, overwrite=TRUE) # gzip
  
  message(paste0(id_data, " done!"))
}

safeParallel()