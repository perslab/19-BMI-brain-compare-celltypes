library("Matrix")
library("here")
library("data.table")
library("dplyr")
library("hdf5r", lib.loc = "/home/cbmr/lhv464/R/x86_64-pc-linux-gnu-library/3.5")
library("Seurat", lib.loc = "/projects/jonatan/R-lib-custom/")
library("loomR")

path_loomfile = "/scratch/data-for_fast_access/pub-others/zeisel-biorxiv-2018/l5_all.loom"
path_fileout_rawData = here("tmp-data","expression", "mousebrain.umi.csv")
path_fileout_metaData = here("tmp-data","expression", "mousebrain.metadata.csv")

ds = loomR::connect(path_loomfile, mode= "r")
#
dt_metaData = data.table(
    cell_id = ds[["col_attrs/CellID"]][],
  ClusterName = ds[["col_attrs/ClusterName"]][]
)
  
# NB: matrix is transposed in loomR object

mat_raw <- ds[["matrix"]][,]

colnames(mat_raw) <- ds[["row_attrs/Gene"]][]
rownames(mat_raw) <- ds[["col_attrs/CellID"]][]

# write to disk

mat_raw %>% t %>% as.data.table -> dt_raw
dt_raw <- data.table("gene"=colnames(mat_raw), dt_raw)
data.table::fwrite(dt_raw, file=path_fileout_rawData,  # fwrite cannot write gziped files
                   nThread=48, verbose=T) # write file ---> write to scratch
R.utils::gzip(path_fileout_rawData, overwrite=TRUE) # gzip

data.table::fwrite(dt_metaData, file=path_fileout_metaData,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch

# close the connection to the loom file
ds$close_all()

print("Mousebrain preprocessing done!")
