#1. 合并1#######################################################################
# 定义器官列表
organs <- c("肺", "肾", "脑", "肌肉", "睾丸", "胃", "肝", "骨", "心", "血管", "结肠", "前列腺", "皮肤","血")

# 基础路径
base_path <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results"

# 对每个器官进行处理
for(organ in organs) {
  # 设置当前器官的工作目录
  organ_path <- file.path(base_path, organ)
  
  # 检查目录是否存在
  if(!dir.exists(organ_path)) {
    cat("警告：目录不存在：", organ_path, "\n")
    next
  }
  
  # 设置工作目录
  setwd(organ_path)
  
  # 获取所有文件
  files <- list.files(pattern = "gsvascore.csv$")
  
  if(length(files) == 0) {
    cat("警告：在", organ, "目录中没有找到CSV文件\n")
    next
  }
  
  # 获取唯一的药物名称
  drugs <- unique(sapply(strsplit(files, "_"), `[`, 1))
  
  cat("开始处理", organ, "器官的数据...\n")
  
  # 对每种药物进行处理
  for(drug in drugs) {
    # 找到该药物的所有文件
    drug_files <- files[grep(paste0("^", drug, "_"), files)]
    
    # 创建空列表存储数据框
    dfs <- list()
    
    # 读取每个文件
    for(file in drug_files) {
      df <- read.csv(file)
      dfs[[file]] <- df
    }
    
    # 合并所有数据框
    merged_df <- do.call(rbind, dfs)
    
    # 保存合并后的文件
    output_filename <- paste0(drug, "_merged_gsvascore.csv")
    write.csv(merged_df, output_filename, row.names = FALSE)
    
    cat("  已完成", drug, "在", organ, "中的文件合并，保存为", output_filename, "\n")
  }
  
  cat("完成", organ, "器官的所有处理\n\n")
}

cat("所有器官的文件处理完成！\n")

#2. 合并2#######################################################################

library(readr)
# 肺
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺/Cyclophosphamide_merged_gsvascore.csv")

lung_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                      list(file1, file2, file3, file4, file5, file6, file7,
                           file8, file9, file10, file11, file12, file13, file14,
                           file15, file16, file17))
colnames(lung_ssgsea)[1] <- "pathway"
colnames(lung_ssgsea)[2:(ncol(lung_ssgsea))] <- paste0(colnames(lung_ssgsea)[2:(ncol(lung_ssgsea))], "_fei")

# 肾
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾/Cyclophosphamide_merged_gsvascore.csv")

kidney_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                        list(file1, file2, file3, file4, file5, file6, file7,
                             file8, file9, file10, file11, file12, file13, file14,
                             file15, file16, file17))
colnames(kidney_ssgsea)[1] <- "pathway"
colnames(kidney_ssgsea)[2:(ncol(kidney_ssgsea))] <- paste0(colnames(kidney_ssgsea)[2:(ncol(kidney_ssgsea))], "_shen")

# 脑
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑/Cyclophosphamide_merged_gsvascore.csv")

brain_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                       list(file1, file2, file3, file4, file5, file6, file7,
                            file8, file9, file10, file11, file12, file13, file14,
                            file15, file16, file17))
colnames(brain_ssgsea)[1] <- "pathway"
colnames(brain_ssgsea)[2:(ncol(brain_ssgsea))] <- paste0(colnames(brain_ssgsea)[2:(ncol(brain_ssgsea))], "_nao")

# 肌肉
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉/Cyclophosphamide_merged_gsvascore.csv")

muscle_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                        list(file1, file2, file3, file4, file5, file6, file7,
                             file8, file9, file10, file11, file12, file13, file14,
                             file15, file16, file17))
colnames(muscle_ssgsea)[1] <- "pathway"
colnames(muscle_ssgsea)[2:(ncol(muscle_ssgsea))] <- paste0(colnames(muscle_ssgsea)[2:(ncol(muscle_ssgsea))], "_ji")

# 心
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心/Cyclophosphamide_merged_gsvascore.csv")

heart_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                       list(file1, file2, file3, file4, file5, file6, file7,
                            file8, file9, file10, file11, file12, file13, file14,
                            file15, file16, file17))
colnames(heart_ssgsea)[1] <- "pathway"
colnames(heart_ssgsea)[2:(ncol(heart_ssgsea))] <- paste0(colnames(heart_ssgsea)[2:(ncol(heart_ssgsea))], "_xin")

# 肝
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝/Cyclophosphamide_merged_gsvascore.csv")

liver_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                       list(file1, file2, file3, file4, file5, file6, file7,
                            file8, file9, file10, file11, file12, file13, file14,
                            file15, file16, file17))
colnames(liver_ssgsea)[1] <- "pathway"
colnames(liver_ssgsea)[2:(ncol(liver_ssgsea))] <- paste0(colnames(liver_ssgsea)[2:(ncol(liver_ssgsea))], "_gan")

# 睾丸
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸/Cyclophosphamide_merged_gsvascore.csv")

testis_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                       list(file1, file2, file3, file4, file5, file6, file7,
                            file8, file9, file10, file11, file12, file13, file14,
                            file15, file16, file17))
colnames(testis_ssgsea)[1] <- "pathway"
colnames(testis_ssgsea)[2:(ncol(testis_ssgsea))] <- paste0(colnames(testis_ssgsea)[2:(ncol(testis_ssgsea))], "_gao")

# 胃
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃/Cyclophosphamide_merged_gsvascore.csv")

stomach_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                        list(file1, file2, file3, file4, file5, file6, file7,
                             file8, file9, file10, file11, file12, file13, file14,
                             file15, file16, file17))
colnames(stomach_ssgsea)[1] <- "pathway"
colnames(stomach_ssgsea)[2:(ncol(stomach_ssgsea))] <- paste0(colnames(stomach_ssgsea)[2:(ncol(stomach_ssgsea))], "_wei")

# 骨
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨/Cyclophosphamide_merged_gsvascore.csv")

bone_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                         list(file1, file2, file3, file4, file5, file6, file7,
                              file8, file9, file10, file11, file12, file13, file14,
                              file15, file16, file17))
colnames(bone_ssgsea)[1] <- "pathway"
colnames(bone_ssgsea)[2:(ncol(bone_ssgsea))] <- paste0(colnames(bone_ssgsea)[2:(ncol(bone_ssgsea))], "_gu")


# 血管
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管/Cyclophosphamide_merged_gsvascore.csv")

vessel_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                      list(file1, file2, file3, file4, file5, file6, file7,
                           file8, file9, file10, file11, file12, file13, file14,
                           file15, file16, file17))
colnames(vessel_ssgsea)[1] <- "pathway"
colnames(vessel_ssgsea)[2:(ncol(vessel_ssgsea))] <- paste0(colnames(vessel_ssgsea)[2:(ncol(vessel_ssgsea))], "_guan")

# 结肠
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠/Cyclophosphamide_merged_gsvascore.csv")

colon_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                        list(file1, file2, file3, file4, file5, file6, file7,
                             file8, file9, file10, file11, file12, file13, file14,
                             file15, file16, file17))
colnames(colon_ssgsea)[1] <- "pathway"
colnames(colon_ssgsea)[2:(ncol(colon_ssgsea))] <- paste0(colnames(colon_ssgsea)[2:(ncol(colon_ssgsea))], "_jie")

# 前列腺
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺/Cyclophosphamide_merged_gsvascore.csv")

prostate_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                       list(file1, file2, file3, file4, file5, file6, file7,
                            file8, file9, file10, file11, file12, file13, file14,
                            file15, file16, file17))
colnames(prostate_ssgsea)[1] <- "pathway"
colnames(prostate_ssgsea)[2:(ncol(prostate_ssgsea))] <- paste0(colnames(prostate_ssgsea)[2:(ncol(prostate_ssgsea))], "_qian")

# 皮肤
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤/Cyclophosphamide_merged_gsvascore.csv")

skin_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                          list(file1, file2, file3, file4, file5, file6, file7,
                               file8, file9, file10, file11, file12, file13, file14,
                               file15, file16, file17))
colnames(skin_ssgsea)[1] <- "pathway"
colnames(skin_ssgsea)[2:(ncol(skin_ssgsea))] <- paste0(colnames(skin_ssgsea)[2:(ncol(skin_ssgsea))], "_pi")

# 血
file1 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Temozolomide_merged_gsvascore.csv")
file2 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Vincristine_merged_gsvascore.csv")
file3 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Irinotecan_merged_gsvascore.csv")
file4 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Methotrexate_merged_gsvascore.csv")
file5 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Oxaliplatin_merged_gsvascore.csv")
file6 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Paclitaxel_merged_gsvascore.csv")
file7 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/PBS_merged_gsvascore.csv")
file8 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Cytarabine_merged_gsvascore.csv")
file9 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Docetaxel_merged_gsvascore.csv")
file10 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Doxorubicine_merged_gsvascore.csv")
file11 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Etoposide_merged_gsvascore.csv")
file12 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Fluorouracil_merged_gsvascore.csv")
file13 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Gemcitabine_merged_gsvascore.csv")
file14 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Capecitabine_merged_gsvascore.csv")
file15 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Carboplatin_merged_gsvascore.csv")
file16 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Cisplatin_merged_gsvascore.csv")
file17 <- read_csv("E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血/Cyclophosphamide_merged_gsvascore.csv")

blood_ssgsea <- Reduce(function(x, y) merge(x, y, by = "X"), 
                      list(file1, file2, file3, file4, file5, file6, file7,
                           file8, file9, file10, file11, file12, file13, file14,
                           file15, file16, file17))
colnames(blood_ssgsea)[1] <- "pathway"
colnames(blood_ssgsea)[2:(ncol(blood_ssgsea))] <- paste0(colnames(blood_ssgsea)[2:(ncol(blood_ssgsea))], "_blood")

library(dplyr)
# 合并所有数据
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
  inner_join(skin_ssgsea, by = "pathway") %>%
  inner_join(blood_ssgsea, by = "pathway")

# 检查结果
dim(all_ssgsea)
head(all_ssgsea)

row.names(all_ssgsea) <- all_ssgsea$pathway
all_ssgsea <- all_ssgsea[,c(2:1019)]

library(tidyr)  # 加载tidyr包
library(limma)
library(dplyr)
library(readxl)

     