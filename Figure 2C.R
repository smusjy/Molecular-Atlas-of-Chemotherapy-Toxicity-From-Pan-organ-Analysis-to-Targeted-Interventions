#PART 1:数据清洗################################################################

#载入数据
library(readr)
exp_chemo <- read_csv("exp_chemo.csv")

#制作小鼠的ENSEMBL转SYMBOL表

# 读取GTF文件
library(rtracklayer)
library(dplyr)
gtf <- import('gencode.vM36.annotation.gtf')  # 替换为GTF文件路径

# 转换为数据框并提取需要的信息
gtf_df <- as.data.frame(gtf)

# 提取基因信息
gtf_df <- gtf_df[gtf_df$type == "gene",]

id_conversion <- gtf_df[,c(10,11,12)]
id_conversion$ENSEMBL <- id_conversion$gene_id
id_conversion$SYMBOL <- id_conversion$gene_name
id_conversion <- id_conversion[,c(2,4,5)]
table(id_conversion$gene_type)
id_conversion <- id_conversion[id_conversion$gene_type == "protein_coding",]
id_conversion <- id_conversion[,c(2,3)]

# 去除版本号（如果需要）
id_conversion$ENSEMBL <- sub("\\.\\d+$", "", id_conversion$ENSEMBL)

# 更换基因名
exp_chemo <- left_join(exp_chemo, id_conversion, by = "ENSEMBL")
exp_chemo <- exp_chemo[!is.na(exp_chemo$SYMBOL),]

organ_specific_top5_genes_simple <- read.csv("organ_specific_top5_genes_simple.csv")

exp_chemo <- exp_chemo[exp_chemo$SYMBOL %in% organ_specific_top5_genes_simple$gene,]

exp_chemo <- exp_chemo[,c(2:1360)]

save(exp_chemo, file = "exp_chemo.RData")

#PART 2: 导入数据，批量进行limma################################################
# 加载必要的包
library(tidyverse)
library(uwot)
library(ggplot2)
library(ggrepel)

# 创建空数据框存储logFC结果
logfc_df <- matrix(nrow = nrow(exp_chemo), ncol = 0)

# 获取所有器官名称
organs <- c("blood", "xin", "guan", "fei", "gan", "gao", "nao", 
            "gu", "jiechang", "ji", "wei", "pi", "qian", "shen")

# 主循环：对每个药物-器官组合进行分析
for(drug in unique(drug_organ$Drug)) {
  # 跳过PBS组
  if(drug == "PBS") next
  
  drug_info <- drug_organ[drug_organ$Drug == drug, ]
  
  for(organ in organs) {
    # 获取药物处理组的样本名
    drug_samples <- unlist(strsplit(drug_info[[organ]], ", "))
    
    # 特殊处理blood和jiechang的样本名
    if(organ == "blood") {
      drug_samples_suffix <- paste0(drug_samples, ".blood")
      pbs_samples <- c("EAP.PBS1", "EAP.PBS2", "EAP.PBS3", "EAP.PBS4")
    } else if(organ == "jiechang") {
      drug_samples_suffix <- paste0(drug_samples, ".jie")
      pbs_samples <- grep("\\.jie$", 
                          colnames(exp_chemo)[grep("^P", colnames(exp_chemo))],
                          value = TRUE)
    } else {
      drug_samples_suffix <- paste0(drug_samples, ".", organ)
      pbs_samples <- grep(paste0("\\.", organ, "$"), 
                          colnames(exp_chemo)[grep("^P", colnames(exp_chemo))],
                          value = TRUE)
    }
    
    # 创建样本名到实际列名的映射
    names(drug_samples_suffix) <- drug_samples
    
    # 检查样本是否存在于数据中
    drug_samples_suffix <- drug_samples_suffix[drug_samples_suffix %in% colnames(exp_chemo)]
    pbs_samples <- pbs_samples[pbs_samples %in% colnames(exp_chemo)]
    
    if(length(drug_samples_suffix) > 0 && length(pbs_samples) > 0) {
      # 计算对照组的平均表达值
      control_mean <- rowMeans(exp_chemo[, pbs_samples, drop=FALSE])
      
      # 对每个药物处理组样本单独计算logFC
      for(i in seq_along(drug_samples_suffix)) {
        # 获取单个样本的表达值
        drug_exp <- exp_chemo[, drug_samples_suffix[i]]
        
        # 计算log2FC
        logFC <- log2(drug_exp) - log2(control_mean)
        
        # 使用原始样本编号（不带后缀）作为标识符
        sample_id <- names(drug_samples_suffix)[i]
        
        # 创建新的列名：编号_药物_器官
        column_name <- paste(sample_id, drug, organ, sep="_")
        
        # 添加到数据框
        logfc_df <- cbind(logfc_df, logFC)
        colnames(logfc_df)[ncol(logfc_df)] <- column_name
      }
    }
  }
}

drug_colors <- c(
  "5-Fluorouracil" = "#FFD33A",
  "Capecitabine" = "#6F9AE8",
  "Carboplatin" = "#8C9196",
  "Cisplatin" = "#D4B17F",
  "Cyclophosphamide" = "#FC7444",
  "Cytarabine" = "#D7E49D",
  "Docetaxel" = "#1A7DC0",
  "Doxorubicine" = "#EE5A3E",
  "Etoposide" = "#467330",
  "Gemcitabine" = "#70D1F4",
  "Irinotecan" = "#370332",
  "Methotrexate" = "#0A5049",
  "Oxaliplatin" = "#C50912",
  "Paclitaxel" = "#92321E",
  "Temozolomide" = "#1A9893",
  "Vincristine" = "#FE8DBF"
)

# 从列名中创建metadata
metadata <- data.frame(
  sample = colnames(logfc_df),
  stringsAsFactors = FALSE) %>%
  separate(sample, c("sample_id", "drug", "organ"), sep="_")

# 修改药物名称
metadata$drug[metadata$drug == "Fluorouracil"] <- "5-Fluorouracil"

# 按照drug_colors的顺序排列数据
metadata$drug <- factor(metadata$drug, levels = names(drug_colors))

# 创建器官名称映射字典
organ_mapping <- c(
  "shen" = "Kidney",
  "nao" = "Brain",
  "ji" = "SkeletalMuscle",
  "gao" = "Testis",
  "fei" = "Lung",
  "wei" = "Stomach",
  "gan" = "Liver",
  "gu" = "Bone",
  "xin" = "Heart",
  "guan" = "BloodVessel",
  "jiechang" = "Colon",
  "qian" = "Prostate",
  "pi" = "Skin",
  "blood" = "Blood"
)

# 转换器官名称为英文
metadata$organ_en <- organ_mapping[metadata$organ]

# 转置矩阵使样本为行
data_for_umap <- t(logfc_df)

# 运行UMAP
set.seed(42)
umap_result <- umap(data_for_umap,
                    n_neighbors = 8,
                    min_dist = 1.0,
                    spread = 1.5,
                    repulsion_strength = 2,
                    negative_sample_rate = 10)

# 创建用于绘图的数据框
plot_data <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  Organ = metadata$organ_en,
  Drug = metadata$drug,
  Sample = metadata$sample_id
)

# 设置颜色
organ_colors <- setNames(
  c("#92321E","#D4B17F", "#8C9196", "#1A9893", "#D7E49D", "#C50912",
    "#370332", "#FFD33A", "#1A7DC0", "#70D1F4", "#FC7444",
    "#467330", "#FE8DBF", "#6F9AE8"),
  c("Blood", "BloodVessel", "Bone", "Brain", "Colon", "Heart", 
    "Kidney", "Liver", "Lung", "Prostate",
    "SkeletalMuscle","Skin", "Stomach", "Testis")
)

# 创建基础主题
base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.margin = unit(c(1,1,1,1), "cm")
  )

# 按器官绘图（无标签，有图例）
p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Organ)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = organ_colors) +
  labs(title = "Samples Grouped by Organ", 
       color = "Organ") +
  base_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

# 按药物绘图（无标签，有图例）
p2 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Drug)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = drug_colors) +
  labs(title = "Samples Grouped by Drug", 
       color = "Drug") +
  base_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

# 按器官绘图（有标签）
p3 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Organ)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = organ_colors) +
  geom_text_repel(aes(label = Sample), 
                  size = 2,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = Inf) +
  labs(title = "Samples Grouped by Organ", 
       color = "Organ") +
  base_theme +
  theme(legend.position = "right")

# 按药物绘图（有标签）
p4 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Drug)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = drug_colors) +
  geom_text_repel(aes(label = Sample), 
                  size = 2,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = Inf) +
  labs(title = "Samples Grouped by Drug", 
       color = "Drug") +
  base_theme +
  theme(legend.position = "right")

# 创建无图例版本的图
p1_no_legend <- p1 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")

# 保存图形
ggsave("UMAP_by_organ_with_legend.pdf", p1, width = 8, height = 6)
ggsave("UMAP_by_drug_with_legend.pdf", p2, width = 8, height = 6)
ggsave("UMAP_by_organ_no_legend.pdf", p1_no_legend, width = 6, height = 6)
ggsave("UMAP_by_drug_no_legend.pdf", p2_no_legend, width = 6, height = 6)


# 保存结果对象
save(logfc_df, metadata, plot_data, umap_result, file = "logFC_UMAP_results.RData")