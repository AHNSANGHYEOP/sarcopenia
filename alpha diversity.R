### new way alpha diversity

View(level_6r)
View(level_7r)

rm(alpha_6r)
rm(alpha_7r)
rm(alpha_diversity_6r)
rm(alpha_diversity_7r)

## 1. new data frame

# level_6r alpha diversity

level_6r$Shannon_Diversity <- NA
level_6r$Simpson_Diversity <- NA
level_6r$Species_Richness <- NA

for (i in 1:nrow(level_6r)) {
  level_6r$Shannon_Diversity[i] <- diversity(level_6r[i, 11:1415], index = "shannon")
  level_6r$Simpson_Diversity[i] <- diversity(level_6r[i, 11:1415], index = "simpson")
  level_6r$Species_Richness[i] <- specnumber(level_6r[i, 11:1415])
}

View(level_6r)

alpha_6r<-level_6r[,-c(11:1415)]
View(alpha_6r)

# level_7r alpha diversity
ncol(level_7r)

level_7r$Shannon_Diversity <- NA
level_7r$Simpson_Diversity <- NA
level_7r$Species_Richness <- NA

for (i in 1:nrow(level_7r)) {
  level_7r$Shannon_Diversity[i] <- diversity(level_7r[i, 11:3832], index = "shannon")
  level_7r$Simpson_Diversity[i] <- diversity(level_7r[i, 11:3832], index = "simpson")
  level_7r$Species_Richness[i] <- specnumber(level_7r[i, 11:3832])
}

View(level_7r)

alpha_7r<-level_7r[,-c(11:3832)]
View(alpha_7r)

## 2. use filter() for ggplot visualization

# 0week alpha_6r

alpha_6r_0week <- alpha_6r |>
  filter(visit == "0week")

# 6week alpha_6r

alpha_6r_6week <- alpha_6r |>
  filter(visit == "6week")

# 12week alpha_6r

alpha_6r_12week <- alpha_6r |>
  filter(visit == "12week")

# 0week alpha_7r

alpha_7r_0week <- alpha_7r |>
  filter(visit == "0week")

# 6week alpha_7r

alpha_7r_6week <- alpha_7r |>
  filter(visit == "6week")

# 12week alpha_7r

alpha_7r_12week <- alpha_7r |>
  filter(visit == "12week")

## 3. use boxplot and point(position_jitter) for visualization

# 0week alpha_6r

plot_shannon_alpha_6r_0week<-ggplot(alpha_6r_0week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 6r 0week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_6r_0week<-ggplot(alpha_6r_0week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 6r 0week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_6r_0week<-ggplot(alpha_6r_0week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 6r 0week by Group",
       x = "Group",
       y = "Species_Richness")

# 6week alpha_6r

plot_shannon_alpha_6r_6week<-ggplot(alpha_6r_6week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 6r 6week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_6r_6week<-ggplot(alpha_6r_6week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 6r 6week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_6r_6week<-ggplot(alpha_6r_6week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 6r 6week by Group",
       x = "Group",
       y = "Species_Richness")

# 12week alpha_6r

plot_shannon_alpha_6r_12week<-ggplot(alpha_6r_12week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 6r 12week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_6r_12week<-ggplot(alpha_6r_12week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 6r 12week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_6r_12week<-ggplot(alpha_6r_12week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 6r 12week by Group",
       x = "Group",
       y = "Species_Richness")

# 0week alpha_7r

plot_shannon_alpha_7r_0week<-ggplot(alpha_7r_0week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 7r 0week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_7r_0week<-ggplot(alpha_7r_0week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 7r 0week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_7r_0week<-ggplot(alpha_7r_0week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 7r 0week by Group",
       x = "Group",
       y = "Species_Richness")

# 6week alpha_7r

plot_shannon_alpha_7r_6week<-ggplot(alpha_7r_6week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 7r 6week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_7r_6week<-ggplot(alpha_7r_6week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 7r 6week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_7r_6week<-ggplot(alpha_7r_6week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 7r 6week by Group",
       x = "Group",
       y = "Species_Richness")

# 12week alpha_7r

plot_shannon_alpha_7r_12week<-ggplot(alpha_7r_12week, aes(x = group, y = Shannon_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Shannon Index 7r 12week by Group",
       x = "Group",
       y = "Shannon Diversity Index")

plot_simpson_alpha_7r_12week<-ggplot(alpha_7r_12week, aes(x = group, y = Simpson_Diversity)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Simpson Index 7r 12week by Group",
       x = "Group",
       y = "Simpson Diversity Index")

plot_richness_alpha_7r_12week<-ggplot(alpha_7r_12week, aes(x = group, y = Species_Richness)) +
  geom_boxplot(fill = "skyblue") +
  geom_point(position = position_jitter(width = 0.2), color = "darkblue", alpha = 0.6)+
  theme_minimal() +
  labs(title = "Species_Richness 7r 12week by Group",
       x = "Group",
       y = "Species_Richness")

# 6r shannon

plot_grid(
  plot_shannon_alpha_6r_0week,plot_shannon_alpha_6r_6week,plot_shannon_alpha_6r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 6r simpson

plot_grid(
  plot_simpson_alpha_6r_0week,plot_simpson_alpha_6r_6week,plot_simpson_alpha_6r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 6r richness

plot_grid(
  plot_richness_alpha_6r_0week,plot_richness_alpha_6r_6week,plot_richness_alpha_6r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 7r shannon

plot_grid(
  plot_shannon_alpha_7r_0week,plot_shannon_alpha_7r_6week,plot_shannon_alpha_7r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 7r simpson

plot_grid(
  plot_simpson_alpha_7r_0week,plot_simpson_alpha_7r_6week,plot_simpson_alpha_7r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 7r richness

plot_grid(
  plot_richness_alpha_7r_0week,plot_richness_alpha_7r_6week,plot_richness_alpha_7r_12week,
  labesl="AUTO",
  ncol=2,
  align="hv"
)

# 정량적 계산 해보기
# 정량적 계산은 max와 min, mean 값 정도만 구해서 비교하면 되는지
# 따로 함수를 써서 계산해야 하는지.











