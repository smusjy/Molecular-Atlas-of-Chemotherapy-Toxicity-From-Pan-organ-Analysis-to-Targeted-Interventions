#1. 计算TCR#####################################################################
library(immunarch)
library(dplyr)
library(tools)

# 设置工作目录
setwd("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 3/TCR/Mm_TCR_data/Chemo_Data_tsv")

# 获取所有TSV文件
tsv_files <- list.files(pattern = "\\.tsv$")

# 创建一个空的数据框来存储结果
results <- data.frame(
  File_Name = character(),
  Shannon_Index = numeric(),
  Simpson_Index = numeric(),
  Chao1 = numeric(),
  stringsAsFactors = FALSE
)

# 处理每个文件
for(file in tsv_files) {
  # 读取TSV文件
  data <- try(read.delim(file, header=TRUE, sep="\t", stringsAsFactors = FALSE))
  
  if(!inherits(data, "try-error")) {
    # TCR数据格式化
    tcr_formatted <- data %>% 
      filter(grepl("^TR[A-Z][A-Z]", V) | grepl("^TR[A-Z]C", C)) %>%
      select(
        V.name = V,
        J.name = J,
        D.name = D,
        Clones = X.count,
        Proportion = frequency,
        CDR3.nt = CDR3nt,
        CDR3.aa = CDR3aa
      ) %>%
      mutate(Sample = "TCR_Sample")
    
    tcr_formatted <- unique(tcr_formatted)
    
    # 检查数据框是否为空或只包含NA
    if(nrow(tcr_formatted) > 0 && !all(is.na(tcr_formatted$Clones))) {
      # 重新计算并标准化proportion
      tcr_formatted_norm <- tcr_formatted %>%
        mutate(
          Proportion = Clones / sum(Clones)
        )
      
      # 计算Shannon指数
      shannon_index <- tryCatch({
        tcr_formatted_norm %>%
          select(Proportion) %>%
          pull(Proportion) %>%
          entropy(., .base = exp(1), .norm = FALSE, .do.norm = NA, .laplace = 1e-12)
      }, error = function(e) NA)
      
      # 计算Simpson指数
      simpson_index <- tryCatch({
        repDiversity(tcr_formatted_norm, "gini.simp")
      }, error = function(e) NA)
      
      # 计算Chao1指数
      chao1_value <- tryCatch({
        div_result <- repDiversity(tcr_formatted_norm, .method = "chao1")
        as.numeric(div_result[1])
      }, error = function(e) NA)
      
    } else {
      shannon_index <- NA
      simpson_index <- NA
      chao1_value <- NA
    }
    
    # 添加结果到数据框
    results <- rbind(results, data.frame(
      File_Name = tools::file_path_sans_ext(file),
      Shannon_Index = shannon_index,
      Simpson_Index = simpson_index,
      Chao1 = chao1_value,
      stringsAsFactors = FALSE
    ))
  }
}

# 打印结果
print(results)

# 保存结果到CSV文件
write.csv(results, "diversity_indices_results.csv", row.names = FALSE)

# 查看前几行结果
head(results)

#2. 差异分析Pandrug#############################################################
library(dplyr)
library(tidyr)

# 读取drug_organ.csv
drug_organ <- read.csv("drug_organ.csv", stringsAsFactors = FALSE)

# 转换为长格式
drug_organ_long <- drug_organ %>%
  pivot_longer(
    cols = -Drug,
    names_to = "Organ_CN",
    values_to = "Identifiers"
  ) %>%
  separate_rows(Identifiers, sep = ",\\s*") %>%
  filter(Identifiers != "")

# 英文器官名称对照
organ_mapping <- c(
  blood = "Blood",
  xin = "Heart",
  guan = "BloodVessel",
  fei = "Lung",
  gan = "Liver",
  gao = "Testis",
  nao = "Brain",
  gu = "Bone",
  jiechang = "Colon",
  ji = "SkeletalMuscle",
  wei = "Stomach",
  pi = "Skin",
  qian = "Prostate",
  shen = "Kidney"
)

drug_organ_long <- drug_organ_long %>%
  mutate(Organ = organ_mapping[Organ_CN])

# 读取多样性指数结果
diversity_results <- read.csv("diversity_indices_results.csv", stringsAsFactors = FALSE)

# 解析文件名
parsed_samples <- diversity_results %>%
  mutate(
    base_name = sub("_report$", "", File_Name),
    parts = strsplit(base_name, "-"),
    Identifier = sapply(parts, function(x) x[1]),
    Organ_CN = sapply(parts, function(x) ifelse(length(x) > 1, x[2], "whole"))
  ) %>%
  select(-base_name, -parts)

parsed_samples$Organ_CN[parsed_samples$Organ_CN == "whole"] <- "blood"

# 合并药物信息
merged_data <- parsed_samples %>%
  left_join(
    drug_organ_long,
    by = c("Identifier" = "Identifiers", "Organ_CN" = "Organ_CN")
  ) %>%
  filter(!is.na(Drug))

# 计算每个器官和分组的平均值和log2FC
grouped_means <- merged_data %>%
  mutate(Group = ifelse(Drug == "PBS", "PBS", "Chemo")) %>%
  group_by(Organ, Group) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index, na.rm = TRUE),
    Mean_Simpson = mean(Simpson_Index, na.rm = TRUE),
    Mean_Chao1 = mean(Chao1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Mean_Shannon, Mean_Simpson, Mean_Chao1)
  ) %>%
  mutate(
    Shannon_log2FC = log2(Mean_Shannon_Chemo / Mean_Shannon_PBS),
    Simpson_log2FC = log2(Mean_Simpson_Chemo / Mean_Simpson_PBS),
    Chao1_log2FC = log2(Mean_Chao1_Chemo / Mean_Chao1_PBS)
  ) %>%
  arrange(Organ)

# 保存结果
write.csv(grouped_means, "Organ_Diversity_log2FC.csv", row.names = FALSE)

#3. 差异分析Panorgan############################################################
library(dplyr)
library(tidyr)

# 使用之前已经处理好的merged_data

# 计算每种药物的整体多样性指数
drug_diversity <- merged_data %>%
  group_by(Drug) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index, na.rm = TRUE),
    Mean_Simpson = mean(Simpson_Index, na.rm = TRUE),
    Mean_Chao1 = mean(Chao1, na.rm = TRUE)
  ) %>%
  ungroup()

# 获取PBS的平均值
pbs_values <- drug_diversity %>%
  filter(Drug == "PBS") %>%
  select(Mean_Shannon, Mean_Simpson, Mean_Chao1)

# 计算log2FC
drug_log2fc <- drug_diversity %>%
  mutate(
    Shannon_log2FC = log2(Mean_Shannon / pbs_values$Mean_Shannon),
    Simpson_log2FC = log2(Mean_Simpson / pbs_values$Mean_Simpson),
    Chao1_log2FC = log2(Mean_Chao1 / pbs_values$Mean_Chao1)
  ) %>%
  filter(Drug != "PBS") %>%  # 移除PBS组
  select(Drug, Shannon_log2FC, Simpson_log2FC, Chao1_log2FC) %>%
  arrange(Drug)

# 保存结果
write.csv(drug_log2fc, "Drug_Overall_Diversity_log2FC.csv", row.names = FALSE)