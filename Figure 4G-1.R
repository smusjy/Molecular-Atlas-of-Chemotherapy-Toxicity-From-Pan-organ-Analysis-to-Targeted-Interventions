# 加载必要的包
library(limma)
library(dplyr)
library(readxl)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 4/ssGSEA"

# -----------------------------
# 加载各组织ssGSEA数据并添加器官后缀
# -----------------------------
# 1. 肺
load(file.path(base_dir, "肺/ssGSEA.Rdata"))
lung_ssgsea <- ssgsea
lung_ssgsea$pathway <- row.names(lung_ssgsea)
colnames(lung_ssgsea)[1:(ncol(lung_ssgsea)-1)] <- paste0(colnames(lung_ssgsea)[1:(ncol(lung_ssgsea)-1)], "_fei")

# 2. 肾
load(file.path(base_dir, "肾/ssGSEA.Rdata"))
kidney_ssgsea <- ssgsea
kidney_ssgsea$pathway <- row.names(kidney_ssgsea)
colnames(kidney_ssgsea)[1:(ncol(kidney_ssgsea)-1)] <- paste0(colnames(kidney_ssgsea)[1:(ncol(kidney_ssgsea)-1)], "_shen")

# 3. 脑
load(file.path(base_dir, "脑/ssGSEA.Rdata"))
brain_ssgsea <- ssgsea
brain_ssgsea$pathway <- row.names(brain_ssgsea)
colnames(brain_ssgsea)[1:(ncol(brain_ssgsea)-1)] <- paste0(colnames(brain_ssgsea)[1:(ncol(brain_ssgsea)-1)], "_nao")

# 4. 肌肉
load(file.path(base_dir, "肌肉/ssGSEA.Rdata"))
muscle_ssgsea <- ssgsea
muscle_ssgsea$pathway <- row.names(muscle_ssgsea)
colnames(muscle_ssgsea)[1:(ncol(muscle_ssgsea)-1)] <- paste0(colnames(muscle_ssgsea)[1:(ncol(muscle_ssgsea)-1)], "_ji")

# 5. 睾丸
load(file.path(base_dir, "睾丸/ssGSEA.Rdata"))
testis_ssgsea <- ssgsea
testis_ssgsea$pathway <- row.names(testis_ssgsea)
colnames(testis_ssgsea)[1:(ncol(testis_ssgsea)-1)] <- paste0(colnames(testis_ssgsea)[1:(ncol(testis_ssgsea)-1)], "_gao")

# 6. 胃
load(file.path(base_dir, "胃/ssGSEA.Rdata"))
stomach_ssgsea <- ssgsea
stomach_ssgsea$pathway <- row.names(stomach_ssgsea)
colnames(stomach_ssgsea)[1:(ncol(stomach_ssgsea)-1)] <- paste0(colnames(stomach_ssgsea)[1:(ncol(stomach_ssgsea)-1)], "_wei")

# 7. 肝
load(file.path(base_dir, "肝/ssGSEA.Rdata"))
liver_ssgsea <- ssgsea
liver_ssgsea$pathway <- row.names(liver_ssgsea)
colnames(liver_ssgsea)[1:(ncol(liver_ssgsea)-1)] <- paste0(colnames(liver_ssgsea)[1:(ncol(liver_ssgsea)-1)], "_gan")

# 8. 骨
load(file.path(base_dir, "骨/ssGSEA.Rdata"))
bone_ssgsea <- ssgsea
bone_ssgsea$pathway <- row.names(bone_ssgsea)
colnames(bone_ssgsea)[1:(ncol(bone_ssgsea)-1)] <- paste0(colnames(bone_ssgsea)[1:(ncol(bone_ssgsea)-1)], "_gu")

# 9. 心
load(file.path(base_dir, "心/ssGSEA.Rdata"))
heart_ssgsea <- ssgsea
heart_ssgsea$pathway <- row.names(heart_ssgsea)
colnames(heart_ssgsea)[1:(ncol(heart_ssgsea)-1)] <- paste0(colnames(heart_ssgsea)[1:(ncol(heart_ssgsea)-1)], "_xin")

# 10. 血管
load(file.path(base_dir, "血管/ssGSEA.Rdata"))
vessel_ssgsea <- ssgsea
vessel_ssgsea$pathway <- row.names(vessel_ssgsea)
colnames(vessel_ssgsea)[1:(ncol(vessel_ssgsea)-1)] <- paste0(colnames(vessel_ssgsea)[1:(ncol(vessel_ssgsea)-1)], "_guan")

# 11. 结肠
load(file.path(base_dir, "结肠/ssGSEA.Rdata"))
colon_ssgsea <- ssgsea
colon_ssgsea$pathway <- row.names(colon_ssgsea)
colnames(colon_ssgsea)[1:(ncol(colon_ssgsea)-1)] <- paste0(colnames(colon_ssgsea)[1:(ncol(colon_ssgsea)-1)], "_jie")

# 12. 前列腺
load(file.path(base_dir, "前列腺/ssGSEA.Rdata"))
prostate_ssgsea <- ssgsea
prostate_ssgsea$pathway <- row.names(prostate_ssgsea)
colnames(prostate_ssgsea)[1:(ncol(prostate_ssgsea)-1)] <- paste0(colnames(prostate_ssgsea)[1:(ncol(prostate_ssgsea)-1)], "_qian")

# 13. 皮肤
load(file.path(base_dir, "皮肤/ssGSEA.Rdata"))
skin_ssgsea <- ssgsea
skin_ssgsea$pathway <- row.names(skin_ssgsea)
colnames(skin_ssgsea)[1:(ncol(skin_ssgsea)-1)] <- paste0(colnames(skin_ssgsea)[1:(ncol(skin_ssgsea)-1)], "_pi")

# -----------------------------
# 合并所有数据
# -----------------------------
all_ssgsea <- lung_ssgsea %>%
  inner_join(kidney_ssgsea, by = "pathway") %>%
  inner_join(brain_ssgsea, by = "pathway") %>%
  inner_join(muscle_ssgsea, by = "pathway") %>%
  inner_join(testis_ssgsea, by = "pathway") %>%
  inner_join(stomach_ssgsea, by = "pathway") %>%
  inner_join(liver_ssgsea, by = "pathway") %>%
  inner_join(bone_ssgsea, by = "pathway") %>%
  inner_join(heart_ssgsea, by = "pathway") %>%
  inner_join(vessel_ssgsea, by = "pathway") %>%
  inner_join(colon_ssgsea, by = "pathway") %>%
  inner_join(prostate_ssgsea, by = "pathway") %>%
  inner_join(skin_ssgsea, by = "pathway")

# 检查合并结果
dim(all_ssgsea)
head(all_ssgsea)

row.names(all_ssgsea) <- all_ssgsea$pathway
# 这里假设数据前96列为相关ssGSEA数值，最后一列为pathway
all_ssgsea <- all_ssgsea[, c(1:96, 98:ncol(all_ssgsea))]

# -----------------------------
# 导入药物-器官分组表
# 表中应包含 "Drug" 字段与其每个器官的有效样本编号（多个样本以 ", " 分隔）
# -----------------------------
drug_organ <- read_xlsx("drug_organ.xlsx")

# -----------------------------
# 获取所有有效样本（只提取样本编号部分）
# -----------------------------
all_valid_samples <- drug_organ %>%
  select(-Drug) %>%
  unlist() %>%
  strsplit(., ", ") %>%
  unlist() %>%
  unique()

# -----------------------------
# 获取所有器官后缀（样本名后缀以 "_" 分隔）
# -----------------------------
organs <- unique(sub(".*_(.*)", "\\1", colnames(all_ssgsea)))

# -----------------------------
# 创建完整的有效样本名列表（样本编号_器官）
# -----------------------------
valid_full_names <- unlist(lapply(all_valid_samples, function(x) paste0(x, "_", organs)))

# -----------------------------
# 筛选 all_ssgsea，只保留有效样本的列
# -----------------------------
all_ssgsea_filtered <- all_ssgsea[, colnames(all_ssgsea) %in% valid_full_names]

# -----------------------------
# 获取药物名称列表，过滤掉PBS（只针对化疗药物分析）
# -----------------------------
drugs <- unique(drug_organ$Drug)
drugs <- drugs[drugs != "PBS"]

# 设置平移时的安全边际 delta
delta <- 0.1

# -----------------------------
# 固定提取PBS组样本
# -----------------------------
pbs_samples <- drug_organ %>%
  filter(Drug == "PBS") %>%
  select(-Drug) %>%
  unlist() %>%
  strsplit(., ", ") %>%
  unlist() %>%
  unique()

pbs_full_names <- unlist(lapply(pbs_samples, function(x) paste0(x, "_", organs)))
pbs_matrix <- all_ssgsea_filtered[, colnames(all_ssgsea_filtered) %in% pbs_full_names, drop = FALSE]
pbs_mean_orig <- rowMeans(pbs_matrix)

# -----------------------------
# 循环分析各化疗药物，每种药物均与固定的PBS控制组比较
# -----------------------------
all_results <- list()

for(drug in drugs) {
  # 提取当前药物对应的样本
  drug_samples <- drug_organ %>%
    filter(Drug == drug) %>%
    select(-Drug) %>%
    unlist() %>%
    strsplit(., ", ") %>%
    unlist() %>%
    unique()
  
  drug_full_names <- unlist(lapply(drug_samples, function(x) paste0(x, "_", organs)))
  
  # 从 all_ssgsea_filtered 中提取药物组数据
  drug_matrix <- all_ssgsea_filtered[, colnames(all_ssgsea_filtered) %in% drug_full_names, drop = FALSE]
  
  # ----------- 根据泛器官方法计算 log₂FC -----------
  # 合并当前药物与固定PBS数据计算全局最小值
  all_values <- c(as.numeric(as.matrix(drug_matrix)),
                  as.numeric(as.matrix(pbs_matrix)))
  min_val <- min(all_values, na.rm = TRUE)
  
  if (min_val <= 0) {
    k <- -min_val + delta
  } else {
    k <- 0
  }
  
  # 对两组数据加上平移常数 k
  adjusted_drug <- drug_matrix + k
  adjusted_pbs  <- pbs_matrix + k
  
  # 计算加平移后两组数据的均值，并基于 log₂ 转换计算 log₂FC
  drug_mean_adj <- rowMeans(adjusted_drug)
  pbs_mean_adj  <- rowMeans(adjusted_pbs)
  
  log2FC <- log2(drug_mean_adj) - log2(pbs_mean_adj)
  # ----------- 结束 log₂FC 计算 -----------
  
  # 同时计算原始（未加 k）数据的组均值，用于后续参考
  drug_mean_orig <- rowMeans(drug_matrix)
  
  # 构建最终结果数据框：包括通路、log₂FC、处理组原始均值及PBS组原始均值
  results <- data.frame(
    pathway = rownames(all_ssgsea_filtered),
    log2FC = log2FC,
    Mean_Treatment = drug_mean_orig,
    Mean_PBS = pbs_mean_orig,
    stringsAsFactors = FALSE
  )
  
  # 保存结果到列表中，并写入 CSV 文件
  all_results[[drug]] <- results
  output_file <- paste0(drug, "_results.csv")
  write.csv(results, file = output_file, row.names = FALSE)
  
  # 打印进度
  cat("Completed analysis for:", drug, "\n")
}

cat("All drug analyses completed.\n")

# -------------------------------
# 新增：验算各个通路在每个器官中 E8、E9、E10 三列的均值
# -------------------------------
# 思路：对于 all_ssgsea（合并后原始数据），在每个器官中提取列名中包含 E8、E9 和 E10 的样本，
#        然后对每个通路计算这三列的均值。结果保存为一个数据框，每一列对应一个器官。

# 创建结果数据框，第一列为-pathway-
organ_verification <- data.frame(pathway = rownames(all_ssgsea), stringsAsFactors = FALSE)

for(org in organs) {
  # 构造正则表达式：列名必须以 "E8", "E9" 或 "E10" 开头，接着 "_" 和当前器官后缀
  pattern <- paste0("^(E8|E9|E10)_", org, "$")
  sel_cols <- grep(pattern, colnames(all_ssgsea), value = TRUE)
  
  if(length(sel_cols) > 0) {
    # 对选中列计算行均值（相当于将 E8, E9, E10 三列加起来后求均值）
    organ_verification[[org]] <- rowMeans(all_ssgsea[, sel_cols, drop = FALSE])
  } else {
    organ_verification[[org]] <- NA
  }
}

# 写入验证结果 CSV 文件
write.csv(organ_verification, "organ_verification.csv", row.names = FALSE)
cat("Verification results saved to organ_verification.csv\n")