##1. 肺#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肺"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,5)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_lung_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_lung_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Lung_results.csv")

##2. 血#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,2)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_blood_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_blood_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Blood_results.csv")
##3. 心#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/心"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,3)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_xin_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_xin_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Heart_results.csv")
##4. 血管#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/血管"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,4)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_guan_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_guan_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_BloodVessel_results.csv")
##5. 肝#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肝"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,6)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_gan_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_gan_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Liver_results.csv")
##6. 睾丸#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/睾丸"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,7)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_gao_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_gao_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Testis_results.csv")
##7. 脑#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/脑"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,8)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_nao_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_nao_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Brain_results.csv")
##8. 骨#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/骨"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,9)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_gu_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_gu_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Bone_results.csv")
##9. 结肠#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/结肠"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,10)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_jiechang_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_jiechang_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Colon_results.csv")
##10. 肌肉#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肌肉"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,11)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_ji_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_ji_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_SkeletalMuscle_results.csv")
##11. 胃#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/胃"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,12)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_wei_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_wei_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Stomach_results.csv")


##12. 皮肤#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/皮肤"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,13)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_pi_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_pi_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Skin_results.csv")


##13. 前列腺#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/前列腺"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,14)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_qian_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_qian_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Prostate_results.csv")


##14. 肾#########################################################################
library(limma)
library(readr)
library(tidyverse)

# 设置工作目录
base_dir <- "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/ssGSEA/Results/肾"

# 导入药物-器官分组表
drug_organ <- read_csv("drug_organ.csv")
current_drug_organ <- drug_organ[,c(1,15)]

treat <- current_drug_organ$Drug[c(1:13,15:17)]
treat_organ <- current_drug_organ[c(1:13,15:17),]
control_organ <- current_drug_organ[14,]

# 读取并合并所有数据
all_data <- NULL
file_types <- c("m2.cp.reactome", "m5.go", "mh.all")

# 先处理对照组数据
control_data <- NULL
for(type in file_types) {
  pattern <- paste0("PBS_", type, ".*_shen_gsvascore.csv$")
  control_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  
  if(length(control_file) > 0) {
    temp_data <- read.csv(control_file, check.names = FALSE)
    rownames(temp_data) <- temp_data[,1]
    temp_data <- temp_data[,-1]
    
    if(is.null(control_data)) {
      control_data <- temp_data
    } else {
      control_data <- rbind(control_data, temp_data)
    }
  }
}

if(is.null(control_data)) {
  stop("Control data not found!")
}

# 处理治疗组数据
for(drug in treat) {
  drug_data <- NULL
  
  for(type in file_types) {
    pattern <- paste0(drug, "_", type, ".*_shen_gsvascore.csv$")
    drug_file <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if(length(drug_file) > 0) {
      temp_data <- read.csv(drug_file, check.names = FALSE)
      rownames(temp_data) <- temp_data[,1]
      temp_data <- temp_data[,-1]
      
      if(is.null(drug_data)) {
        drug_data <- temp_data
      } else {
        drug_data <- rbind(drug_data, temp_data)
      }
    }
  }
  
  if(!is.null(drug_data)) {
    if(is.null(all_data)) {
      all_data <- drug_data
    } else {
      all_data <- cbind(all_data, drug_data)
    }
  }
}

# 合并处理组和对照组数据
final_data <- cbind(all_data, control_data)

# 创建分组信息
group <- factor(c(rep("treatment", ncol(all_data)), 
                  rep("control", ncol(control_data))))

# 运行limma
design <- model.matrix(~group)
fit <- lmFit(final_data, design)
fit <- eBayes(fit)

# 获取完整的limma结果
results <- topTable(fit, coef=2, number=Inf)

# 保存结果
write.csv(results, "E:/9. Chemo_AEs_altas/Analysis/Results/Figure 2/Figure 2F/all_drugs_vs_PBS_limma_Kidney_results.csv")

