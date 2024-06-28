
  load("data/cwas.test.data.rda")

  inputdata$PathwayScore1 <- inputdata$S.Score
  inputdata$PathwayScore2 <- inputdata$G2M.Score

  #inputdata <- readRDS("I:\\BaiduSyncdisk\\DT20230523已备份\\DT20230801\\20231104CWAS\\20231212axisFrame\\重新注释sc20231217.rds")
  cwas.test.data <-  cwas_build_model_data(inputdata = inputdata,
                                           Sample = "orig.ident", Phenotype = "Pathology",
                                           Celltype = "celltype",Group = "Tissue",Subgroup = "GSE",
                                           Control_label="Health",Disease_label ="Fibrosis",
                                           genelist = c("TGFB1","TGFB2","TGFB3"),
                                           geneset_score =  c("PathwayScore1","PathwayScore2")
                                           )
cwas.test.data[is.na(cwas.test.data)] <- "No Label"
cwas.test.data.single.model <- cwas_allmodel_cal(cwas.test.data)
cwas.test.data.single.model.best <- cwas_autoselected_model(cwas.test.data.single.model)


usethis::use_data(cwas.test,cwas.test.data,cwas.test.data.single.model,cwas.test.data.single.model.best,overwrite = T)

Target_names <- readRDS("I:\\BaiduSyncdisk\\DT20230523已备份\\DT20230801\\20231104CWAS\\20240519DrugTarget\\Target_names0522.rds")


usethis::use_data(Target_names,overwrite = T)


covid.data <- readRDS("COVID2024-01-17105618.rds")
covid.data$dataset
covid.data$sampleID
covid.data$cellratio <- covid.data$num/covid.data$allnum
covid.data$tissue
covid.data$severity
covid.data <- covid.data[,c("sampleID","celltype","cellratio","dataset","tissue","severity")]
usethis::use_data(covid.data,overwrite = T)


MNP.data <- readRDS("I:\\BaiduSyncdisk\\DT20230523已备份\\DT20230801\\20231104CWAS\\20240605DrugALLcal\\MNP2021_MNP_Verse.RDS")

metadata <- MNP.data@meta.data
unique(metadata$Tissue)
usethis::use_data(metadata,overwrite = T)


gene_names <- paste0("Gene", 1:1000)
cell_names <- paste0("Cell", 1:100)
expression_matrix <- matrix(rpois(1000 * 100, lambda = 5), nrow = 1000, ncol = 100)
rownames(expression_matrix) <- gene_names
colnames(expression_matrix) <- cell_names

# 创建 Seurat 对象
seurat_object <- CreateSeuratObject(counts = expression_matrix, project = "ExampleProject")
seurat_object@meta.data <- metadata

usethis::use_data(seurat_object,overwrite = T)











######################################################################依赖包######################

# 安装和加载 usethis 包（如果尚未安装）
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}



# 加载 usethis 包
library(usethis)

# 定义要添加为建议依赖的包列表
suggest_packages <- c(
  "factoextra",
  "ggpubr",
  "ggrepel",
  "grid",
  "gtools",
  "gtsummary",
  "leaps",
  "MuMIn",
  "simr",
  "stats",
  "reshape2",
  "bruceR"

)

import_packages <- c(
  "dplyr",
  "pheatmap",
  "ggplot2",
  "lme4",
  "forestploter",
  "tidyr",
  "doParallel",
  "plyr",
  "foreach",
  "iterators",
  "igraph",
  "bruceR"

)

# 使用 use_package() 函数将每个包添加为建议依赖
for (pkg in suggest_packages) {
  use_package(pkg, type = "Suggest")
}
# 使用 use_package() 函数将每个包添加为建议依赖
for (pkg in import_packages) {
  use_package(pkg, type = "Import")
}


##########安装###########################

detach("package:SPARKLE", unload = TRUE)

remove.packages("SPARKLE")
# ָ?????س??????ļ???·??

install.packages("rio")
install.packages("dplyr")
#install.packages("magrittr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("ggrepel")
install.packages("leaps")
install.packages("lme4")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("simr")
install.packages("forestploter")
install.packages("gtsummary")
install.packages("MuMIn")
install.packages("tidyr")
install.packages("Hmisc")
install.packages("gtools")
package_path <- "I:\\BaiduSyncdisk\\refdatabase20220619\\Rpackage\\SPARKLE_0.1.8.tar.gz"  # ?滻Ϊ?????ļ?·??

# ??װ???س?????
install.packages(package_path, repos = NULL, type = "source")


devtools::install_github("chenxi199506/SPARKLE")


