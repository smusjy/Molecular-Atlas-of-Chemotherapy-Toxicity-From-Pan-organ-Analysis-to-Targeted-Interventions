library(readr)
SOC_ROR <- read_csv("SOC_ROR.csv")

library(ggplot2)

# 创建分段的标签和颜色
ror_breaks <- c("ROR > 16", "8 < ROR ≤ 16", "4 < ROR ≤ 8", "1 < ROR ≤ 4", "Not Significant")
ror_colors <- c("#A80C3A", "#D24A44", "#FF8B8B", "#FFC1C1", "#C1C1C1")

# 创建ROR分段的颜色函数
get_ror_color <- function(ror, significance) {
  if(significance == "No") return("#C1C1C1")
  else if(ror > 16) return("#A80C3A")
  else if(ror > 8) return("#D24A44")
  else if(ror > 4) return("#FF8B8B")
  else return("#FFC1C1")
}

# 添加颜色列和分类列
SOC_ROR$point_color <- mapply(get_ror_color, 
                              SOC_ROR$ROR, 
                              SOC_ROR$Significance)

SOC_ROR$ror_category <- factor(
  ifelse(SOC_ROR$Significance == "No", "Not Significant",
         ifelse(SOC_ROR$ROR > 16, "ROR > 16",
                ifelse(SOC_ROR$ROR > 8, "8 < ROR ≤ 16",
                       ifelse(SOC_ROR$ROR > 4, "4 < ROR ≤ 8",
                              "1 < ROR ≤ 4")))),
  levels = ror_breaks
)

# 创建排序用的交互项
drugs <- sort(unique(SOC_ROR$Drug))
databases <- unique(SOC_ROR$DataBase)
ordered_levels <- vector()
for(drug in drugs) {
  for(db in databases) {
    ordered_levels <- c(ordered_levels, paste(drug, db))
  }
}

# 创建新的排序因子
SOC_ROR$x_axis <- factor(paste(SOC_ROR$Drug, SOC_ROR$DataBase),
                         levels = ordered_levels)

ggplot(SOC_ROR, aes(x = x_axis, 
                    y = reorder(SOC, -rank(SOC)), # 修改这一行来反转顺序
                    size = Cases)) +
  geom_hline(yintercept = 1:length(unique(SOC_ROR$SOC)), 
             color = "grey90", 
             size = 0.5) +
  geom_point(aes(fill = ror_category), 
             alpha = 1,  
             shape = 21,        
             stroke = 0.25,      
             color = "black") + 
  scale_fill_manual(values = setNames(ror_colors, ror_breaks),
                    name = "ROR Range") +
  scale_size_continuous(range = c(1, 6),
                        breaks = c(500, 5000, 10000, 30000, 50000)) +
  labs(x = "Drug / Database", 
       y = "System Organ Class (SOC)",
       size = "Number of Cases") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))