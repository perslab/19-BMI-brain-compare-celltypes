############### SYNOPSIS ###################
# Run correlation analysis n=on CELLECT and rare / mendelian variant results

### OUTPUT:
#' TODO
#' 
### REMARKS:
# ....

### REFERENCE:
# * 

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--path_CELLECT_results", type="character",
              help = "character, absolute path to prioritization.csv file with CELLECT outs, [default %default]"), 
  make_option("--dir_geneset_results", type="character", 
              help = "character, absolute path to dir containing ES geneset enrichment results, [default %default]"), 
  make_option("--id_datExpr", type="character",
              help = "e.g. 'mousebrain'"), 
  make_option("--geneset_results_regex", type="character",
              help = "string to use alongside id_datExpr to match correct file in dir_geneset_results"), 
  make_option("--col_geneset_results", type="character", default= "p.value",
              help = "which column of the geneset results should be correlated? "), 
  make_option("--col_CELLECT", type="character", default= "pvalue",
              help = "which column of the CELLECT results should be correlated? "), 
  make_option("--vec_geneset_name", type="character", 
              help = "quoted vector of character, select this in geneset_name column of the geneset results csv. If geneset_name column doesn't exist it uses all rows [default $default]"), 
  make_option("--vec_GWAS", type="character",
              help = "quoted vector of GWAS among CELLECT results to include [default $default]"), 
  make_option("--method", type="character", default= "pearson",
              help = "correlation coefficient, passed to stats::cor, one of pearson, kendall, spearman, [default $default]"), 
  make_option("--minLog10transform", type="logical", default= T,
              help = "transform values by -log10 before correlating? [default $default]"), 
  make_option("--output_label", type="character",
              help = "character, label for output files"),
  make_option("--append_results", type="logical", default=TRUE,
              help = "If file already exists, append results? [default %default] ")
)


opt <- parse_args(OptionParser(option_list=option_list))

path_CELLECT_results  <- opt$path_CELLECT_results
dir_geneset_results <- opt$dir_geneset_results
id_datExpr <- opt$id_datExpr
geneset_results_regex <- opt$geneset_results_regex
vec_geneset_name <- eval(parse(text=opt$vec_geneset_name))
col_geneset_results <- opt$col_geneset_results
col_CELLECT <- opt$col_CELLECT
vec_GWAS <- eval(parse(text=opt$vec_GWAS))
method <- opt$method
minLog10transform <- opt$minLog10transform
output_label <- opt$output_label
append_results <- opt$append_results


######################################################################
########################### PACKAGES #################################
######################################################################

message("Loading packages")

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("here"))

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
########################## LOAD FILES ################################
######################################################################

dt_CELLECT <- fread(path_CELLECT_results)

path_geneset_results <- dir(path=dir_geneset_results, 
                            pattern=geneset_results_regex, 
                            full.names=T)
path_geneset_results <- grep(pattern = id_datExpr, 
                             x= path_geneset_results,
                             value=T)
path_geneset_results <- grep(pattern = "\\.csv\\.gz", 
                             x= path_geneset_results,
                             value=T)
dt_geneset_results <- fread(path_geneset_results)

######################################################################
########################## FILTER AND REORDER  #######################
######################################################################

mat_cor <- sapply(vec_GWAS, function(GWAS) {
  sapply(vec_geneset_name, function(geneset_name){
    
    # cell_type,statistic,parameter,p.value,p.value_emp,alternative
    vec_geneset_value <- dt_geneset_results[[col_geneset_results]]
    names(vec_geneset_value) <- dt_geneset_results$cell_type
    
    if (!is.null(dt_geneset_results$geneset_name)) vec_geneset_value <- 
      vec_geneset_value[dt_geneset_results$geneset_name == geneset_name]
    
    # gwas,specificity_id,annotation,tau,se,pvalue
    vec_CELLECT_value <- dt_CELLECT[gwas == GWAS & specificity_id == id_datExpr,][[col_CELLECT]]
    names(vec_CELLECT_value) <- dt_CELLECT[gwas == GWAS & specificity_id == id_datExpr,][["annotation"]]
    
    # check the two vectors
    stopifnot(all(names(vec_CELLECT_value) %in%  names(vec_geneset_value)))
    stopifnot(length(vec_CELLECT_value) == length(vec_geneset_value))
    
    # put vectors in same order
    vec_CELLECT_value <- vec_CELLECT_value[match(names(vec_geneset_value),names(vec_CELLECT_value))]
    
    # logtransform p values
    if (minLog10transform) {
      vec_geneset_value <- -log10(vec_geneset_value)
      vec_CELLECT_value <- -log10(vec_CELLECT_value)
    }
    
    # compute correlation
    cor(x=vec_geneset_value,
        y=vec_CELLECT_value,
        method=method)
  })
})

mat_cor <- matrix(data=mat_cor, nrow=length(vec_geneset_name), ncol = length(vec_GWAS))
colnames(mat_cor) <- vec_GWAS

dt_cor <- data.table(
  geneset_name = vec_geneset_name,
  mat_cor 
)

dt_cor_long <- melt.data.table(dt_cor, 
                    id.vars= "geneset_name",
                    variable.name = "gwas",
                    value.name = paste0(method, "_cor"))

dt_cor_long$minLog10_vals <- minLog10transform
######################################################################
#############################  WRAP UP  ##############################
######################################################################

file.out.results <- here(paste0(output_label, "_cor_CELLECT_rarevariant.csv"))
fwrite(x = dt_cor_long, file =file.out.results,  append = if (file.exists(file.out.results)) append_results else F, nThread=24, verbose=T)

message("script done!")

