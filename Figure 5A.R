library(readr)
library(dplyr)
library(ggplot2)

# 读取数据
wilcox <- read_csv("wilcox.csv")

# 处理特殊字符
wilcox$Indicator <- gsub("\xa6\xc3-GT", "γ-GT", wilcox$Indicator)

# 计算每个指标的显著药物数量
sig_counts <- wilcox %>%
  group_by(Indicator) %>%
  summarize(sig_count = sum(P_value < 0.05, na.rm = TRUE))

# 找到γ-GT的位置，并分别获取左右两侧的指标
all_indicators <- unique(wilcox$Indicator)
gamma_gt_pos <- which(all_indicators == "γ-GT")
left_indicators <- all_indicators[1:gamma_gt_pos]
right_indicators <- all_indicators[(gamma_gt_pos):length(all_indicators)]

# 对左右两侧分别按显著数量排序
left_ordered <- sig_counts %>%
  filter(Indicator %in% left_indicators) %>%
  arrange(desc(sig_count)) %>%
  pull(Indicator)

right_ordered <- sig_counts %>%
  filter(Indicator %in% right_indicators[-1]) %>%  # 排除γ-GT
  arrange(desc(sig_count)) %>%
  pull(Indicator)

# 合并排序后的结果
xorder <- c(left_ordered, right_ordered)

# 将Indicator转换为factor并指定顺序
wilcox$Indicator <- factor(wilcox$Indicator, levels = xorder)

# 创建log2FC的分类变量
wilcox <- wilcox %>%
  mutate(log2FC_cat = case_when(
    log2FC < -3 ~ "< -3",
    log2FC >= -3 & log2FC < -2 ~ "-3 to -2",
    log2FC >= -2 & log2FC < -1 ~ "-2 to -1",
    log2FC >= -1 & log2FC < 0 ~ "-1 to 0",
    log2FC >= 0 & log2FC < 1 ~ "0 to 1",
    log2FC >= 1 & log2FC < 2 ~ "1 to 2",
    log2FC >= 2 & log2FC < 3 ~ "2 to 3",
    log2FC >= 3 ~ "> 3"
  ))

# 获取Drug的唯一值，把5-Fluorouracil放在最前面
drug_levels <- unique(wilcox$Drug)
drug_levels <- c("5-Fluorouracil", drug_levels[drug_levels != "5-Fluorouracil"])

# 将Drug转换为factor并反转顺序
wilcox$Drug <- factor(wilcox$Drug, levels = rev(drug_levels))

# 绘制图形
ggplot(wilcox, aes(x=Indicator, y=Drug)) +
  geom_tile(fill="white", color="black",, size=0.4) +
  geom_point(data=subset(wilcox, P_value < 0.05 & !is.na(log2FC)),
             aes(size=-log10(P_value),
                 fill=log2FC_cat),
             shape=22,
             color="black")  +
  scale_fill_manual(values=c(
    "< -3" = "#08306B",
    "-3 to -2" = "#2171B5",
    "-2 to -1" = "#6BAED6",
    "-1 to 0" = "#C6DBEF",
    "0 to 1" = "#FEE0D2",
    "1 to 2" = "#FC9272",
    "2 to 3" = "#DE2D26",
    "> 3" = "#A50F15"
  )) +
  scale_size_continuous(range=c(2,6)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid=element_blank(),
        axis.text=element_text(size=8)) +
  labs(fill="log2FC",
       size="-log10(P-value)")

# 创建简化的上调/下调分类数据
plot_data <- wilcox %>%
  filter(P_value < 0.05) %>%
  mutate(regulation = case_when(
    log2FC > 0 ~ "Up-regulated",
    log2FC < 0 ~ "Down-regulated"
  )) %>%
  group_by(Indicator, regulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  ungroup()

# 计算每个Indicator的总数，并保留regulation信息
total_counts <- plot_data %>%
  group_by(Indicator) %>%
  summarise(total = sum(count)) %>%
  mutate(regulation = "total")  # 添加一个虚拟的regulation分类

# 计算每个柱子中标签的位置
plot_data <- plot_data %>%
  group_by(Indicator) %>%
  arrange(Indicator, regulation) %>%
  mutate(
    pos = cumsum(count),
    pos_middle = pos - (count/2)
  )

# 使用与热图相同的Indicator顺序
plot_data$Indicator <- factor(plot_data$Indicator, levels = xorder)  # 使用xorder而不是xorder_new
total_counts$Indicator <- factor(total_counts$Indicator, levels = xorder)  # 使用xorder而不是xorder_new

ggplot(plot_data, aes(x = Indicator, y = count, fill = regulation)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) + # 将linewidth改为0.2
  geom_text(data = total_counts,
            aes(y = total, label = total),
            vjust = -0.5,
            size = 3) +
  geom_text(aes(y = pos_middle, label = count),
            size = 3,
            color = "black") +
  scale_fill_manual(values = c(
    "Up-regulated" = "#FC9272",
    "Down-regulated" = "#6BAED6"
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Indicator",
    y = "Number of Significant Drugs",
    fill = "Regulation"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))