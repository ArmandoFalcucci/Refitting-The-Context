# Packages used ####

library(tibble)
library(dplyr)
library(forcats)
library(ggplot2)
library(reshape)
library(igraph)


# Build the correlation table ####

correlation.connections <- read.csv("data/RF_Connections_Correlation.csv")

correlation.connections.table <- correlation.connections %>%
  janitor::tabyl(US.a, US.b)

correlation.connections.table <- correlation.connections.table %>%
  add_column(D1c = "0") %>%
  add_column(`D3a+b` = "0")

correlation.connections.table$`D3a+b` <- as.numeric(correlation.connections.table$`D3a+b`)
correlation.connections.table$D1c <- as.numeric(correlation.connections.table$D1c)

correlation.connections.table <- correlation.connections.table %>%
  add_row(US.a = "D3a+b sabbie", A2 = 0, A2R = 0, A1 = 0, D6 = 0, `D3d base` = 0, D3d = 0, `D3b alpha` = 0, D3b = 0, `D3+D6` = 0, `D3a+b sabbie` = 0, `D3a+b` = 0, D1c = 0)

correlation.connections.table$US.a <- as.factor(correlation.connections.table$US.a)

correlation.connections.table$`D3a+b sabbie` <- as.numeric(correlation.connections.table$`D3a+b sabbie`)

correlation.connections.table.reshaped <- correlation.connections.table %>%
  relocate(US.a, A2, A2R, A1, D6, `D3d base`, D3d, `D3b alpha`, D3b, `D3+D6`, D1c, `D3a+b`, `D3a+b sabbie`) %>%
  arrange(factor(US.a, levels = c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D1c", "D3a+b", "D3a+b sabbie"))) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))

correlation.connections.table.reshaped[-1] <- get.adjacency(
  graph_from_data_frame(
    get.data.frame(
      graph_from_adjacency_matrix(
        as.matrix(correlation.connections.table.reshaped[-1]), "directed"
      )
    ), FALSE
  ),
  type = "upper",
  sparse = FALSE
)

correlation.connections.table.reshaped.2 <- correlation.connections.table.reshaped[,-1]

correlation.connections.table.reshaped.2[lower.tri(correlation.connections.table.reshaped.2)] <- NA

correlation.connections.table.reshaped.2 <- as.data.frame(correlation.connections.table.reshaped.2)

US <- c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D1c", "D3a+b", "D3a+b sabbie")

correlation.connections.table.reshaped.2$US <- US

correlation.connections.table.reshaped.2 <- na.omit(melt(correlation.connections.table.reshaped.2, 'US'))

correlation.connections.table.reshaped.2$US <- factor(correlation.connections.table.reshaped.2$US, levels = c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D1c", "D3a+b", "D3a+b sabbie"))

correlation.connections.table.reshaped.2$variable <- factor(correlation.connections.table.reshaped.2$variable, levels = c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D1c", "D3a+b", "D3a+b sabbie"))


# Build the heatmap figure ####

heatmap.connections <- correlation.connections.table.reshaped.2 %>%
  filter(value > 0.1) %>%
  mutate(US = fct_relevel(US, c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D3a+b sabbie", "D3a+b", "D1c"))) %>%
  mutate(variable = fct_relevel(variable, c("A2", "A2R", "A1", "D6", "D3d base", "D3d", "D3b alpha", "D3b", "D3+D6", "D3a+b sabbie", "D3a+b", "D1c"))) %>%
  ggplot(aes(US, variable)) +
  ggtitle('') +
  ylab('SU') +
  xlab('SU') +
  labs(fill = "Log. value") +
  geom_tile(aes(fill = log.value), color='white', alpha = 0.8) +
  scale_fill_viridis_c(direction = 1) +
  scale_x_discrete(drop=F) + #used to add the values with 0 connections
  geom_text(aes(US, variable, label = value), color = "white", size = 5) +
  #theme(panel.background = element_blank()) +
  # theme(panel.background = element_blank()) +
  theme_bw() +
  #theme_minimal()
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), 
        legend.text=element_text(size=11), legend.title=element_text(size=11), legend.key.size=unit(0.8,"cm")) +
  theme(text = element_text(size = 16, face = "bold"))

heatmap.connections


# Save figure ####

ggsave(heatmap.connections, 
       filename = file.path("output/figures/raw_figures", "Heatmap.raw.tiff"),
       width = 22, 
       height = 14,
       dpi = 300, 
       units = "cm", 
       device = 'tiff')

# **Important note**: Since I was not able to add the squares with 0 values to complete the diagonal, I further modified this figure to add them manually, leaving the PDF file with the modified layers well visible.