# 安装并加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# 创建数据框
data <- read_excel("Figure 1B.xlsx")

# 设置自定义颜色
group_colors <- c(
  "Anti-metabolites" = "#F9CA2B",
  "Anthracyclines" = "#CA931D",
  "Alkylating agents" = "#00A9B4",
  "Taxanes" = "#B4D1A9",
  "Topoisomerase inhibitors" = "#E16965",
  "Vinca alkaloids" = "#8089BA"
)

# 分别创建FAERS和VigiBase的数据框
data_faers <- data %>%
  filter(Database == "FAERS") %>%
  arrange(Cases)

data_vigibase <- data %>%
  filter(Database == "VigiBase") %>%
  mutate(Drugname = factor(Drugname, levels = data_faers$Drugname))

# 创建FAERS横向条形图
p1 <- ggplot(data_faers, aes(x = reorder(Drugname, Cases), y = Cases, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +  # 减小柱子宽度
  geom_text(aes(label = format(Cases, big.mark = ",")), 
            hjust = -0.1,    # 调整文字位置
            size = 3) +      # 调整文字大小
  coord_flip() +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "",
    y = "Number of cases",
    title = "FAERS database"
  ) +
  scale_y_continuous(labels = scales::comma) +
  expand_limits(y = max(data_faers$Cases) * 1.15)  # 扩展y轴范围以适应标签

# 创建VigiBase横向条形图
p2 <- ggplot(data_vigibase, aes(x = Drugname, y = Cases, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +  # 减小柱子宽度
  geom_text(aes(label = format(Cases, big.mark = ",")), 
            hjust = -0.1,    # 调整文字位置
            size = 3) +      # 调整文字大小
  coord_flip() +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "",
    y = "Number of cases",
    title = "VigiBase database"
  ) +
  scale_y_continuous(labels = scales::comma) +
  expand_limits(y = max(data_vigibase$Cases) * 1.15)  # 扩展y轴范围以适应标签

# 分别保存图片
ggsave("Figure 1B_FAERS.pdf", p1, width = 5, height = 7, dpi = 300)
ggsave("Figure 1B_VigiBase.pdf", p2, width = 5, height = 7, dpi = 300)