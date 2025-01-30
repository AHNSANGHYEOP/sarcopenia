### beta diversity

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(cowplot)

## 1. beta diversity 분석법

# beta diversity=ratio species diversity 샘플 "간" 다양성 

# Jaccard distance : 종의 존재유무
# Bray-Curtis distance : 종의 총량, Abundance 고려
# unweighted UniFrac distance : 미생물의 계통 거리에 영향을 받으면서,  
# 집단 간 미생물의 존재/부존재에 따라 측정
# weighted UniFrac distance : 미생물의 계통 거리에 영향을 받으면서, 
# 동시에 분포에 따라 측정

## 2. beta diversity index functions

#vegdist()
#method="bray" or method="jaccard"
#cmdscale()

## 3. PCoA plot try

# beta diversity의 결과는 PCoA(Principal Coordinates Analysis) plot으로 
# 표시할 수 있다. PCoA plot이란 차원 축소를 통해 샘플 사이의 거리를 좌표에
# 표시한 것이다. 좌표상에서 거리가 가까울수록 샘플 간의 다양성이 비슷하다는 
# 것을 알 수 있다. 

help(vegdist)
# Dissimilarity index, partial match to "manhattan", "euclidean", 
# "canberra", "clark", "bray", "kulczynski", "jaccard", "gower",
# "altGower", "morisita", "horn", "mountford", "raup", "binomial", 
# "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", 
# "aitchison", or "robust.aitchison".

# 위 metrics 중 자주 사용되는건 "bray", "jaccard"

# binary metrics는 rare taxa에서의 variations에 민감하고 
# (Jaccard, unweighted UniFrac), quantitative metrics는 abundant taxa에서의 
# variation에 민감하다 (Bray-Curtis, weighted UniFrac). 
# 그래서 샘플에 rare taxa가 많다면 전자, 아니라면 후자를 사용한다.
# jaccard는 이진 데이터의 유사성에, bray는 두 샘플 간 양적 차이에 기반한다.

# UniFrac의 경우 phylogenetic tree 고려하여 계산,
# phylogenetic tree = evolutionary relationship
# vegan package 대신에 phyloseq 또는 GUniFrac package 사용해야 한다.

# bray_6r

bray_dist_6r<-vegdist(level_6r[,c(11:1415)], method="bray")

bray_result_6r <- cmdscale(bray_dist_6r, eig = TRUE, k = 2)  

bray_df_6r <- data.frame(
  Sample = rownames(level_6r),
  PC1 = bray_result_6r$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_6r$points[, 2],  # 두 번째 주성분
  Group = level_6r$group 
)

eig_values_bray_6r <- bray_result_6r$eig / sum(bray_result_6r$eig)

ggplot(bray_df_6r, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance",
    x = paste0("PC1 (", round(eig_values_bray_6r[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_6r[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard_6r

jaccard_dist_6r<-vegdist(level_6r[,c(11:1415)], method="jaccard")

jaccard_result_6r <- cmdscale(jaccard_dist_6r, eig = TRUE, k = 2)  

jaccard_df_6r <- data.frame(
  Sample = rownames(level_6r),
  PC1 = jaccard_result_6r$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_6r$points[, 2],  # 두 번째 주성분
  Group = level_6r$group 
)

eig_values_ja_6r <- jaccard_result_6r$eig / sum(jaccard_result_6r$eig)

ggplot(jaccard_df_6r, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance",
    x = paste0("PC1 (", round(eig_values_ja_6r[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_ja_6r[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# bray_7r
View(level_7r)
bray_dist_7r<-vegdist(level_7r[,c(11:3832)], method="bray")

bray_result_7r <- cmdscale(bray_dist_7r, eig = TRUE, k = 2)  

bray_df_7r <- data.frame(
  Sample = rownames(level_7r),
  PC1 = bray_result_7r$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_7r$points[, 2],  # 두 번째 주성분
  Group = level_7r$group 
)

eig_values_bray_7r <- bray_result_7r$eig / sum(bray_result_7r$eig)

ggplot(bray_df_7r, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance",
    x = paste0("PC1 (", round(eig_values_bray_7r[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_7r[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard_7r

jaccard_dist_7r<-vegdist(level_7r[,c(11:3832)], method="jaccard")

jaccard_result_7r <- cmdscale(jaccard_dist_7r, eig = TRUE, k = 2)  

jaccard_df_7r <- data.frame(
  Sample = rownames(level_7r),
  PC1 = jaccard_result_7r$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_7r$points[, 2],  # 두 번째 주성분
  Group = level_7r$group 
)

eig_values_ja_7r <- jaccard_result_7r$eig / sum(jaccard_result_7r$eig)

ggplot(jaccard_df_7r, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance",
    x = paste0("PC1 (", round(eig_values_ja_7r[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_ja_7r[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# bray 6r 0week

level_6r_0week<-level_6r|>
  filter(visit=="0week")

bray_dist_6r_0week<-vegdist(level_6r_0week[,c(11:1415)], method="bray")

bray_result_6r_0week <- cmdscale(bray_dist_6r_0week, eig = TRUE, k = 2)  

bray_df_6r_0week <- data.frame(
  Sample = level_6r_0week$ID,
  PC1 = bray_result_6r_0week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_6r_0week$points[, 2],  # 두 번째 주성분
  Group = level_6r_0week$group 
)

eig_values_bray_6r_0week <- bray_result_6r_0week$eig / sum(bray_result_6r_0week$eig)

PCoA_6r_bray_0week<-ggplot(bray_df_6r_0week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 6r 0week",
    x = paste0("PC1 (", round(eig_values_bray_6r_0week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_6r_0week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")


# bray 6r 6week

level_6r_6week<-level_6r|>
  filter(visit=="6week")

bray_dist_6r_6week<-vegdist(level_6r_6week[,c(11:1415)], method="bray")

bray_result_6r_6week <- cmdscale(bray_dist_6r_6week, eig = TRUE, k = 2)  

bray_df_6r_6week <- data.frame(
  Sample = level_6r_6week$ID,
  PC1 = bray_result_6r_6week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_6r_6week$points[, 2],  # 두 번째 주성분
  Group = level_6r_6week$group 
)

eig_values_bray_6r_6week <- bray_result_6r_6week$eig / sum(bray_result_6r_6week$eig)

PCoA_6r_bray_6week<-ggplot(bray_df_6r_6week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 6r 6week",
    x = paste0("PC1 (", round(eig_values_bray_6r_6week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_6r_6week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")


# bray 6r 12week

level_6r_12week<-level_6r|>
  filter(visit=="12week")

bray_dist_6r_12week<-vegdist(level_6r_12week[,c(11:1415)], method="bray")

bray_result_6r_12week <- cmdscale(bray_dist_6r_12week, eig = TRUE, k = 2)  

bray_df_6r_12week <- data.frame(
  Sample = level_6r_12week$ID,
  PC1 = bray_result_6r_12week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_6r_12week$points[, 2],  # 두 번째 주성분
  Group = level_6r_12week$group 
)

eig_values_bray_6r_12week <- bray_result_6r_12week$eig / sum(bray_result_6r_12week$eig)

PCoA_6r_bray_12week<-ggplot(bray_df_6r_12week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 6r 12week",
    x = paste0("PC1 (", round(eig_values_bray_6r_12week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_6r_12week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 6r 0week

jaccard_dist_6r_0week<-vegdist(level_6r_0week[,c(11:1415)], method="jaccard")

jaccard_result_6r_0week <- cmdscale(jaccard_dist_6r_0week, eig = TRUE, k = 2)  

jaccard_df_6r_0week <- data.frame(
  Sample = level_6r_0week$ID,
  PC1 = jaccard_result_6r_0week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_6r_0week$points[, 2],  # 두 번째 주성분
  Group = level_6r_0week$group 
)

eig_values_jaccard_6r_0week <- jaccard_result_6r_0week$eig / sum(jaccard_result_6r_0week$eig)

PCoA_6r_jaccard_0week<-ggplot(jaccard_df_6r_0week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance in 6r 0week",
    x = paste0("PC1 (", round(eig_values_jaccard_6r_0week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_6r_0week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 6r 6week

jaccard_dist_6r_6week<-vegdist(level_6r_6week[,c(11:1415)], method="jaccard")

jaccard_result_6r_6week <- cmdscale(jaccard_dist_6r_6week, eig = TRUE, k = 2)  

jaccard_df_6r_6week <- data.frame(
  Sample = level_6r_6week$ID,
  PC1 = jaccard_result_6r_6week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_6r_6week$points[, 2],  # 두 번째 주성분
  Group = level_6r_6week$group 
)

eig_values_jaccard_6r_6week <- jaccard_result_6r_6week$eig / sum(jaccard_result_6r_6week$eig)

PCoA_6r_jaccard_6week<-ggplot(jaccard_df_6r_6week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance in 6r 6week",
    x = paste0("PC1 (", round(eig_values_jaccard_6r_6week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_6r_6week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 6r 12week

jaccard_dist_6r_12week<-vegdist(level_6r_12week[,c(11:1415)], method="jaccard")

jaccard_result_6r_12week <- cmdscale(jaccard_dist_6r_12week, eig = TRUE, k = 2)  

jaccard_df_6r_12week <- data.frame(
  Sample = level_6r_12week$ID,
  PC1 = jaccard_result_6r_12week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_6r_12week$points[, 2],  # 두 번째 주성분
  Group = level_6r_12week$group 
)

eig_values_jaccard_6r_12week <- jaccard_result_6r_12week$eig / sum(jaccard_result_6r_12week$eig)

PCoA_6r_jaccard_12week<-ggplot(jaccard_df_6r_12week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Jaccard Distance in 6r 12week",
    x = paste0("PC1 (", round(eig_values_jaccard_6r_12week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_6r_12week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")


# bray 7r 0week

View(level_7r)
level_7r_0week<-level_7r|>
  filter(visit=="0week")

bray_dist_7r_0week<-vegdist(level_7r_0week[,c(11:3832)], method="bray")

bray_result_7r_0week <- cmdscale(bray_dist_7r_0week, eig = TRUE, k = 2)  

bray_df_7r_0week <- data.frame(
  Sample = level_7r_0week$ID,
  PC1 = bray_result_7r_0week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_7r_0week$points[, 2],  # 두 번째 주성분
  Group = level_7r_0week$group 
)

eig_values_bray_7r_0week <- bray_result_7r_0week$eig / sum(bray_result_7r_0week$eig)

PCoA_7r_bray_0week<-ggplot(bray_df_7r_0week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 7r 0week",
    x = paste0("PC1 (", round(eig_values_bray_7r_0week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_7r_0week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# bray 7r 6week

level_7r_6week<-level_7r|>
  filter(visit=="6week")

bray_dist_7r_6week<-vegdist(level_7r_6week[,c(11:3832)], method="bray")

bray_result_7r_6week <- cmdscale(bray_dist_7r_6week, eig = TRUE, k = 2)  

bray_df_7r_6week <- data.frame(
  Sample = level_7r_6week$ID,
  PC1 = bray_result_7r_6week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_7r_6week$points[, 2],  # 두 번째 주성분
  Group = level_7r_6week$group 
)

eig_values_bray_7r_6week <- bray_result_7r_6week$eig / sum(bray_result_7r_6week$eig)

PCoA_7r_bray_6week<-ggplot(bray_df_7r_6week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 7r 6week",
    x = paste0("PC1 (", round(eig_values_bray_7r_6week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_7r_6week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# bray 7r 12week

level_7r_12week<-level_7r|>
  filter(visit=="12week")

bray_dist_7r_12week<-vegdist(level_7r_12week[,c(11:3832)], method="bray")

bray_result_7r_12week <- cmdscale(bray_dist_7r_12week, eig = TRUE, k = 2)  

bray_df_7r_12week <- data.frame(
  Sample = level_7r_12week$ID,
  PC1 = bray_result_7r_12week$points[, 1],  # 첫 번째 주성분
  PC2 = bray_result_7r_12week$points[, 2],  # 두 번째 주성분
  Group = level_7r_12week$group 
)

eig_values_bray_7r_12week <- bray_result_7r_12week$eig / sum(bray_result_7r_12week$eig)

PCoA_7r_bray_12week<-ggplot(bray_df_7r_12week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Distance in 7r 12week",
    x = paste0("PC1 (", round(eig_values_bray_7r_12week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_bray_7r_12week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 7r 0week

jaccard_dist_7r_0week<-vegdist(level_7r_0week[,c(11:3832)], method="jaccard")

jaccard_result_7r_0week <- cmdscale(jaccard_dist_7r_0week, eig = TRUE, k = 2)  

jaccard_df_7r_0week <- data.frame(
  Sample = level_7r_0week$ID,
  PC1 = jaccard_result_7r_0week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_7r_0week$points[, 2],  # 두 번째 주성분
  Group = level_7r_0week$group 
)

eig_values_jaccard_7r_0week <- jaccard_result_7r_0week$eig / sum(jaccard_result_7r_0week$eig)

PCoA_7r_jaccard_0week<-ggplot(jaccard_df_7r_0week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Jaccard in 7r 0week",
    x = paste0("PC1 (", round(eig_values_jaccard_7r_0week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_7r_0week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 7r 6week

jaccard_dist_7r_6week<-vegdist(level_7r_6week[,c(11:3832)], method="jaccard")

jaccard_result_7r_6week <- cmdscale(jaccard_dist_7r_6week, eig = TRUE, k = 2)  

jaccard_df_7r_6week <- data.frame(
  Sample = level_7r_6week$ID,
  PC1 = jaccard_result_7r_6week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_7r_6week$points[, 2],  # 두 번째 주성분
  Group = level_7r_6week$group 
)

eig_values_jaccard_7r_6week <- jaccard_result_7r_6week$eig / sum(jaccard_result_7r_6week$eig)

PCoA_7r_jaccard_6week<-ggplot(jaccard_df_7r_6week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Jaccard in 7r 6week",
    x = paste0("PC1 (", round(eig_values_jaccard_7r_6week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_7r_6week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# jaccard 7r 12week

jaccard_dist_7r_12week<-vegdist(level_7r_12week[,c(11:3832)], method="jaccard")

jaccard_result_7r_12week <- cmdscale(jaccard_dist_7r_12week, eig = TRUE, k = 2)  

jaccard_df_7r_12week <- data.frame(
  Sample = level_7r_12week$ID,
  PC1 = jaccard_result_7r_12week$points[, 1],  # 첫 번째 주성분
  PC2 = jaccard_result_7r_12week$points[, 2],  # 두 번째 주성분
  Group = level_7r_12week$group 
)

eig_values_jaccard_7r_12week <- jaccard_result_7r_12week$eig / sum(jaccard_result_7r_12week$eig)

PCoA_7r_jaccard_12week<-ggplot(jaccard_df_7r_12week, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCoA Plot based on Bray Jaccard in 7r 12week",
    x = paste0("PC1 (", round(eig_values_jaccard_7r_12week[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(eig_values_jaccard_7r_12week[2] * 100, 1), "%)")
  ) +
  stat_ellipse(linetype="dashed")+
  theme(legend.position = "right")

# all 6r bray

plot_grid(
  PCoA_6r_bray_0week,PCoA_6r_bray_6week,PCoA_6r_bray_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# all 6r jaccard

plot_grid(
  PCoA_6r_jaccard_0week,PCoA_6r_jaccard_6week,PCoA_6r_jaccard_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# all 7r bray

plot_grid(
  PCoA_7r_bray_0week,PCoA_7r_bray_6week,PCoA_7r_bray_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# all 7r jaccard

plot_grid(
  PCoA_7r_jaccard_0week,PCoA_7r_jaccard_6week,PCoA_7r_jaccard_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 두 축의 eigenvalues합이 6r bray의 각 주차에서 약 40%, 6r jaccard는 약 30%
# 7r bray는 약 30%, 7r jaccard는 약 20%이다.
# 좋은 PCoA plot은 2~3개의 축이 50%이상을 설명해줘야한다.
# bray보다 eigvalues의 합이 약 10%정도 낮으며, 절대적으로도 값이 높다고
# 볼 수 없는 jaccard를 결과 해석에 포함 시켜도 될지 문제.

# vegdist() 함수 자체는 alpha diversity분석에서 shannon, simpson처럼
# beta diversity의 요약값을 설명해주지 않는다.
# betadisper() 사용해서 추가 분석을 진행해서 분산 중심적 척도를 구하면
# 이를 beta diversity의 요약값을 설명해주는 값으로 사용해도 괜찮은지.