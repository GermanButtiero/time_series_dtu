source("assignment4/functions/functions_exercise2.R")
library(ggplot2)
library(gridExtra)

data <- read.csv("assignment4/transformer_data.csv")
data

#1.1 
# Panel 1: Transformer Temperature
png("assignment4/plots/ex2_1.png", width = 800, height = 600)
p1 <- ggplot(data, aes(x = time, y = Y)) +
  geom_line(color = "red") +
  labs(
    title = "Transformer Temperature (Yt)",
    x = "Time (hours)",
    y = "Temperature (°C)"
  ) +
  theme_bw()

# Panel 2: Air Temperature
p2 <- ggplot(data, aes(x = time, y = Ta)) +
  geom_line(color = "blue") +
  labs(
    title = "Air Temperature (Ta,t)",
    x = "Time (hours)",
    y = "Temperature (°C)"
  ) +
  theme_bw()

# Panel 3: Solar Radiation
p3 <- ggplot(data, aes(x = time, y = S)) +
  geom_line(color = "orange") +
  labs(
    title = "Solar Radiation (Φs,t)",
    x = "Time (hours)",
    y = "Radiation (W/m²)"
  ) +
  theme_bw()

# Panel 4: Load
p4 <- ggplot(data, aes(x = time, y = I)) +
  geom_line(color = "purple") +
  labs(
    title = "Load (ΦI,t)",
    x = "Time (hours)",
    y = "Load (kA)"
  ) +
  theme_bw()

grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
