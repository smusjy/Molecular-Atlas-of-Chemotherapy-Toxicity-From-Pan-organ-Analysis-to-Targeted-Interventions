#PART 1:数据清洗################################################################
#载入数据
load("E:/9. Chemo_AEs_altas/Analysis/Data/细胞系数据/Chemo_ICI_cell_count_data.Rdata")

shen <- exp_chemo_ici[1]
shen <- as.data.frame(shen)
colnames(shen) <- gsub("X091.1_Count.Rdata.", "", colnames(shen))

nao <- exp_chemo_ici[2]
nao <- as.data.frame(nao)
colnames(nao) <- gsub("X091.10_Count.Rdata.", "", colnames(nao))

ji <- exp_chemo_ici[3]
ji <- as.data.frame(ji)
colnames(ji) <- gsub("X091.11_Count.Rdata.", "", colnames(ji))

gao <- exp_chemo_ici[4]
gao <- as.data.frame(gao)
colnames(gao) <- gsub("X091.12_Count.Rdata.", "", colnames(gao))

fei <- exp_chemo_ici[5]
fei <- as.data.frame(fei)
colnames(fei) <- gsub("X091.13_Count.Rdata.", "", colnames(fei))

wei <- exp_chemo_ici[6]
wei <- as.data.frame(wei)
colnames(wei) <- gsub("X091.2_Count.Rdata.", "", colnames(wei))

gan <- exp_chemo_ici[7]
gan <- as.data.frame(gan)
colnames(gan) <- gsub("X091.3_Count.Rdata.", "", colnames(gan))

gu <- exp_chemo_ici[8]
gu <- as.data.frame(gu)
colnames(gu) <- gsub("X091.4_Count.Rdata.", "", colnames(gu))

xin <- exp_chemo_ici[9]
xin <- as.data.frame(xin)
colnames(xin) <- gsub("X091.5_Count.Rdata.", "", colnames(xin))

guan <- exp_chemo_ici[10]
guan <- as.data.frame(guan)
colnames(guan) <- gsub("X091.6_Count.Rdata.", "", colnames(guan))

jie <- exp_chemo_ici[11]
jie <- as.data.frame(jie)
colnames(jie) <- gsub("X091.7_Count.Rdata.", "", colnames(jie))

qian <- exp_chemo_ici[12]
qian <- as.data.frame(qian)
colnames(qian) <- gsub("X091.8_Count.Rdata.", "", colnames(qian))

pi <- exp_chemo_ici[13]
pi <- as.data.frame(pi)
colnames(pi) <- gsub("X091.9_Count.Rdata.", "", colnames(pi))

library(readxl)
drug_organ <- read_excel("drug_organ.xlsx")

head(drug_organ)

#PART 2: 数据分析和可视化################################################
# 加载必要的包
# 加载所需的包
# 加载必要的包
library(Matrix)      # 用于矩阵运算
library(irlba)       # 用于SVD分解
library(uwot)        # UMAP实现
library(tidyverse)   # 数据处理
library(ggplot2)     # 可视化

# 创建数据列表
organ_data_list <- list(
  "shen" = shen,
  "jiaozhi" = nao,
  "ji" = ji,
  "gao" = gao,
  "fei" = fei,
  "wei" = wei,
  "gan" = gan,
  "gu" = gu,
  "xin" = xin,
  "neipi" = guan,
  "jie" = jie,
  "qian" = qian,
  "biao" = pi
)

# 初始化一个空列表来存储每个器官的结果
all_logfc_results <- list()

# 对每个器官进行计算
for(organ_name in names(organ_data_list)) {
  # 获取当前器官的数据
  current_data <- organ_data_list[[organ_name]]
  
  # 找出PBS样本 (E8, E9, E10)
  pbs_samples <- c("E8", "E9", "E10")
  pbs_cols <- grep(paste(pbs_samples, collapse="|"), colnames(current_data), value = TRUE)
  
  # 计算PBS样本的平均值
  pbs_mean <- rowMeans(current_data[, pbs_cols, drop = FALSE])
  
  # 初始化当前器官的结果矩阵
  organ_results <- matrix(nrow = nrow(current_data), ncol = 0)
  rownames(organ_results) <- rownames(current_data)
  
  # 对每个非PBS样本计算logFC
  non_pbs_cols <- setdiff(colnames(current_data), pbs_cols)
  
  for(col in non_pbs_cols) {
    # 获取样本数据
    sample_data <- current_data[, col]
    
    # 计算logFC
    logfc <- log2((sample_data + 1)/(pbs_mean + 1))
    
    # 从drug_organ中找出对应的药物
    sample_id <- sub("_.*$", "", col)
    drug_match <- NULL
    
    for(drug in unique(drug_organ$Drug)) {
      if(drug != "PBS" && !is.na(drug_organ[[organ_name]][drug_organ$Drug == drug])) {
        if(grepl(sample_id, drug_organ[[organ_name]][drug_organ$Drug == drug])) {
          drug_match <- drug
          break
        }
      }
    }
    
    if(!is.null(drug_match)) {
      col_name <- paste(sample_id, drug_match, organ_name, sep="_")
      organ_results <- cbind(organ_results, logfc)
      colnames(organ_results)[ncol(organ_results)] <- col_name
    }
  }
  
  # 存储结果
  if(ncol(organ_results) > 0) {
    all_logfc_results[[organ_name]] <- organ_results
  }
}

# 检查每个器官的结果维度
for(organ_name in names(all_logfc_results)) {
  cat("Organ:", organ_name, "- Dimensions:", dim(all_logfc_results[[organ_name]]), "\n")
}

# 找出共同基因
common_genes <- Reduce(intersect, lapply(all_logfc_results, rownames))
cat("Number of common genes:", length(common_genes), "\n")

# 创建最终的合并矩阵
final_logfc_matrix <- do.call(cbind, lapply(all_logfc_results, function(x) x[common_genes,]))

# 运行UMAP
set.seed(42)
umap_result <- umap(t(final_logfc_matrix),
                    n_neighbors = 8,
                    min_dist = 1.0,
                    spread = 1.5,
                    repulsion_strength = 2,
                    negative_sample_rate = 10)

# 从列名中提取信息创建metadata
metadata <- data.frame(
  sample = colnames(final_logfc_matrix),
  stringsAsFactors = FALSE) %>%
  separate(sample, c("sample_id", "drug", "organ"), sep="_")

# 创建绘图数据框 - 直接使用umap_result
plot_data <- data.frame(
  UMAP1 = umap_result[,1],  # 直接使用矩阵
  UMAP2 = umap_result[,2],  # 直接使用矩阵
  Drug = metadata$drug,
  Organ = metadata$organ,
  Sample = metadata$sample_id
)

# 创建器官名称映射字典
organ_mapping <- c(
  "shen" = "Kidney",
  "jiaozhi" = "Brain",
  "ji" = "SkeletalMuscle",
  "gao" = "Testis",
  "fei" = "Lung",
  "wei" = "Stomach",
  "gan" = "Liver",
  "gu" = "Bone",
  "xin" = "Heart",
  "neipi" = "BloodVessel",
  "jie" = "Colon",
  "qian" = "Prostate",
  "biao" = "Skin"
)

# 转换器官名称为英文
metadata$organ_en <- organ_mapping[metadata$organ]

# 更新plot_data中的器官名称
plot_data$Organ <- organ_mapping[plot_data$Organ]

# 修改药物名称
plot_data$Drug <- gsub("Doxorubicin hydrochloride", "Doxorubicine", plot_data$Drug)

# 设置药物颜色
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

# 设置器官颜色
organ_colors <- setNames(
  c("#D4B17F", "#8C9196", "#1A9893", "#D7E49D", "#C50912",
    "#370332", "#FFD33A", "#1A7DC0", "#70D1F4", "#FC7444",
    "#467330", "#FE8DBF", "#6F9AE8"),
  c("BloodVessel", "Bone", "Brain", "Colon", "Heart", 
    "Kidney", "Liver", "Lung", "Prostate",
    "SkeletalMuscle", "Skin", "Stomach", "Testis")
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

# 绘图部分保持不变，只是使用新的配色方案
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

# 创建无图例版本
p1_no_legend <- p1 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")

# 保存图形
ggsave("UMAP_by_organ_with_legend.pdf", p1, width = 8, height = 6)
ggsave("UMAP_by_drug_with_legend.pdf", p2, width = 8, height = 6)
ggsave("UMAP_by_organ_no_legend.pdf", p1_no_legend, width = 6, height = 6)
ggsave("UMAP_by_drug_no_legend.pdf", p2_no_legend, width = 6, height = 6)


# 1. 检查数据格式和内容
print("数据维度:")
print(dim(final_logfc_matrix))

print("数据类型:")
print(class(final_logfc_matrix))

print("是否有NA值:")
print(sum(is.na(final_logfc_matrix)))

print("数据概要:")
print(summary(as.vector(final_logfc_matrix)))
