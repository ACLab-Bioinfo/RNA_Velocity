set_lib_paths("/home/patrick/miniconda3/envs/RNAvelocity/lib/R/library")
.libPaths()

library(Seurat)
library(SeuratDisk)
library(stringr)
library(velocyto.R)
library(SeuratWrappers)

##### Input loom. file
PW035_E002 <- ReadVelocity(file = "./patrick/data_20221116/PW035_E002/velocyto/PW035_E002.loom")

##### Rename the barcodes with sample id (cd ./public/processed/sample.meta.txt)
colnames(PW035_E002$spliced) <- paste0(sapply(str_split(colnames(PW035_E002$spliced), ":|x"), "[", 2), "-", "96")
colnames(PW035_E002$unspliced) <- paste0(sapply(str_split(colnames(PW035_E002$unspliced), ":|x"), "[", 2), "-", "96")
colnames(PW035_E002$ambiguous) <- paste0(sapply(str_split(colnames(PW035_E002$ambiguous), ":|x"), "[", 2), "-", "96")

##### Convert loom. to seurat object
seu.PW035_E002 <- as.Seurat(x = PW035_E002)

colnames(seu.PW035_E002)

##### Merge seuobj
seu.velocyto.231122 <- merge(seu.velocyto, y = c(seu.PW030_C000, seu.PW030_E000, seu.PW030_C002, seu.PW030_E002, seu.PW031_C000, seu.PW031_E000, seu.PW031_C002, seu.PW031_E002,
                                                 seu.PW032_C000, seu.PW032_E000, seu.PW034_C000, seu.PW034_E000, seu.PW034_C002, seu.PW034_E002, seu.PW035_C000, seu.PW035_E000,
                                                 seu.PW035_C002, seu.PW035_E002))
                                           

saveRDS(seu.velocyto.231122, "/public/processed/seu.velocyto.231122.rds")


colnames(seu.velocyto)


