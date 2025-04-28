library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dendextend)
library(grid)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Results/Limma_Results")

# 手动设置颜色
drug_colors <- c(
  "#FFD33A", "#6F9AE8", "#8C9196", "#D4B17F", "#FC7444",
  "#D7E49D", "#1A7DC0", "#EE5A3E", "#467330", "#70D1F4",
  "#370332", "#0A5049", "#C50912", "#92321E", "#1A9893",
  "#FE8DBF"  # 16种颜色
)

organ_colors <- c(
  "#92321E","#D4B17F", "#8C9196", "#1A9893", "#D7E49D", "#C50912",
  "#370332", "#FFD33A", "#1A7DC0", "#70D1F4", "#FC7444",
  "#467330", "#FE8DBF", "#6F9AE8"
)

# 创建器官和药物名称的对应关系
organ_names <- c(
  "gan" = "Liver",
  "Lung" = "Lung",
  "xin" = "Heart",
  "nao" = "Brain",
  "gao" = "Testis",
  "pi" = "Skin",
  "Blood" = "Blood",
  "wei" = "Stomach",
  "gu" = "Bone",
  "ji" = "Skeletal Muscle",
  "shen" = "Kidney",
  "qian" = "Prostate",
  "guan" = "Blood Vessel",
  "jie" = "Colon"
)

drug_names <- c(
  "Fluorouracil" = "5-Fluorouracil",
  "Doxorubicine" = "Doxorubicine",
  "Paclitaxel" = "Paclitaxel",
  "Cisplatin" = "Cisplatin",
  "Cyclophosphamide" = "Cyclophosphamide",
  "Carboplatin" = "Carboplatin",
  "Capecitabine" = "Capecitabine",
  "Cytarabine" = "Cytarabine",
  "Docetaxel" = "Docetaxel",
  "Etoposide" = "Etoposide",
  "Gemcitabine" = "Gemcitabine",
  "Irinotecan" = "Irinotecan",
  "Methotrexate" = "Methotrexate",
  "Oxaliplatin" = "Oxaliplatin",
  "Temozolomide" = "Temozolomide",
  "Vincristine" = "Vincristine"
)

# 获取所有csv文件
files <- list.files(pattern = ".*_limma_results\\.csv$")

# 创建存储数据的列表
data_list <- list()

# 读取所有数据并提取logFC
for(file in files) {
  # 从文件名提取器官和药物信息
  info <- strsplit(file, "_")[[1]]
  organ <- info[1]
  drug <- info[2]
  
  # 读取数据
  data <- read.csv(file)
  colnames(data)[7] <- "cell_type"
  
  # 创建组合名称
  combination <- paste(organ, drug, sep="_")
  
  # 存储logFC值
  data_list[[combination]] <- data[, c("cell_type", "logFC")]
}

# 创建宽格式数据框
all_cells <- unique(unlist(lapply(data_list, function(x) x$cell_type)))

# 创建矩阵
mat <- matrix(NA, nrow = length(data_list), ncol = length(all_cells))
rownames(mat) <- names(data_list)
colnames(mat) <- all_cells

# 填充矩阵
for(comb in names(data_list)) {
  data <- data_list[[comb]]
  mat[comb, data$cell_type] <- data$logFC
}

# 处理NA值
mat[is.na(mat)] <- 0

# 对每列进行z分数标准化
mat_scaled <- scale(mat, center = TRUE, scale = TRUE)

# 将Z分数限制在[-3, 3]之间
mat_scaled <- pmin(pmax(mat_scaled, -3), 3)

# 清理列名
colnames(mat_scaled) <- gsub("_CIBERSORT", "", colnames(mat_scaled))
colnames(mat_scaled) <- gsub("_", " ", colnames(mat_scaled))

# 创建注释数据框
annotation_df <- data.frame(
  Drug = sub(".*_", "", rownames(mat_scaled)),
  Organ = sub("_.*", "", rownames(mat_scaled)),
  row.names = rownames(mat_scaled)
)

# 替换器官和药物名称
annotation_df$Organ <- organ_names[annotation_df$Organ]
annotation_df$Drug <- drug_names[annotation_df$Drug]

# 获取实际存在的器官和药物名称并按字母顺序排序
unique_organs <- sort(unique(annotation_df$Organ))
unique_drugs <- sort(unique(annotation_df$Drug))

# 为实际存在的类别分配颜色（按字母顺序）
organ_color_subset <- organ_colors[1:length(unique_organs)]
names(organ_color_subset) <- unique_organs

drug_color_subset <- drug_colors[1:length(unique_drugs)]
names(drug_color_subset) <- unique_drugs

# 修改注释数据框中的因子水平，以控制显示顺序
annotation_df$Organ <- factor(annotation_df$Organ, levels = unique_organs)
annotation_df$Drug <- factor(annotation_df$Drug, levels = unique_drugs)

# 创建注释颜色
annotation_colors <- list(
  Drug = drug_color_subset,
  Organ = organ_color_subset
)

# 绘制热图
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
max_abs <- max(abs(mat_scaled), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)

p2 <- pheatmap(mat_scaled,
               color = colors,
               breaks = breaks,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               fontsize = 10,
               fontsize_row = 8,
               main = "Drug-Organ Combination Effects on Immune Cell Populations (Z-score)",
               angle_col = 45,
               border_color = NA,
               annotation_row = annotation_df,
               annotation_colors = annotation_colors,
               annotation_names_row = TRUE,
               show_rownames = FALSE,
               treeheight_row = 30,  # 修改这里，设置为合适的值，比如30
               treeheight_col = 15,  # 列的聚类树高度也设置为相同值
               annotation_legend = TRUE,
               annotation_legend_side = "right",
               annotation_names_side = "right",
               gaps_row = NULL,
               annotation_position = "right")

# 保存热图
pdf("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/Results/Organ+Drug/drug_organ_immune_heatmap_scaled_tree.pdf", 
    width = 8,
    height = 15)
print(p2)
dev.off()