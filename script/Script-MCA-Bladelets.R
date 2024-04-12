#### MCA of complete bladelets ####

# Packages used ####
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(WRS2)
library(DescTools)

# Selecting the lithic subsets and creating the dataset ####

# A2int ####
dataset_A2int <- read_excel("../data/RF_Dataset_Blanks.xlsx") %>%
  filter(unit == "A2") %>%
  filter(square %in% c("95", "100", "105", "110","115", "120", "125")) %>%
  rename(SU = "unit")

dataset_A2int$SU <- c("A2int")


# A2ext ####
dataset_A2ext <- read_excel("../data/RF_Dataset_Blanks.xlsx") %>%
  filter(unit == "A2") %>%
  filter(square %in% c("55","56", "65", "66")) %>%
  rename(SU = "unit")

dataset_A2ext$SU <- c("A2ext")


# A1ext ####
dataset_A1ext <- read_excel("../data/RF_Dataset_Blanks.xlsx") %>%
  filter(unit == "A1") %>%
  filter(square %in% c("55","56", "65", "66")) %>%
  rename(SU = "unit")

dataset_A1ext$SU <- c("A1ext")

# A1 east ####
dataset_A1east <- read_excel("../data/RF_Dataset_Blanks.xlsx") %>%
  filter(unit == "A1") %>%
  filter(square %in% c("62","63", "72", "73")) %>%
  rename(SU = "unit")

dataset_A1east$SU <- c("A1east")


# D3b alpha ####
dataset_D3balpha <- read_excel("../data/RF_Dataset_Blanks.xlsx") %>%
  rename(SU = "unit") %>%
  filter(SU == "D3b.alpha")

# Merging ####
dataset.blank <- rbind(dataset_A2ext, dataset_A2int, dataset_A1ext,dataset_A1east, dataset_D3balpha)

dataset.blank <- dataset.blank %>%
  mutate(SU = fct_relevel(SU, "A2int", "A2ext", "A1ext", "A1east", "D3b.alpha"))

# Data preparation ####
comparison.MCA <- c("A2int-A2ext", "A2int-A1ext", "A2int-A1east", "A2int-D3b alpha", "A2ext-A1ext", "A2ext-A1east", "A2ext-D3b alpha", "A1ext-A1east", "A1ext-D3b alpha", "A1east-D3b alpha")

MCA_bladelet.morpho <- dataset.blank %>%
  filter(class == "Blank", preservation == "Complete", blank == "Bladelet") %>%
  drop_na(length, width, thickness) %>%
  mutate(curvature = fct_lump(curvature, prop = 0.05)) %>%
  mutate(cross.section = fct_lump(cross.section, prop = 0.05)) %>%
  mutate(distal.end = fct_lump(distal.end, prop = 0.05)) %>%
  mutate(scar.pattern = fct_lump(scar.pattern, prop = 0.05)) %>%
  mutate(blank.shape = fct_lump(blank.shape, prop = 0.05)) %>%
  mutate(SU = recode(SU, D3b.alpha = "D3b alpha"))

MCA_bladelet.morpho$elongation <- RcmdrMisc::binVariable(MCA_bladelet.morpho$elongation, bins = 3, method = "natural", labels = c("low", "medium", "high"))
MCA_bladelet.morpho$robustness <- RcmdrMisc::binVariable(MCA_bladelet.morpho$robustness, bins = 3, method = "natural", labels = c("low", "medium", "high"))

MCA_bladelet.morpho <- MCA_bladelet.morpho %>%
  mutate(torsion.simplified = torsion) %>%
  mutate(torsion.simplified = recode(torsion.simplified, `1L` = "yes", `2L` = "yes", `1R` = "yes", `2R` = "yes", `0` = "no"))

MCA_bladelet.morpho <- MCA_bladelet.morpho %>%
  dplyr::select(SU, elongation, robustness, curvature, torsion.simplified, blank.shape, distal.end, scar.pattern) %>%
  na.omit() %>%
  droplevels() %>%
  mutate(SU = recode(SU, D3b.alpha = "D3b alpha"))

# MCA ####

MCA_bladelet.morpho.analysis <- MCA(MCA_bladelet.morpho, graph = FALSE, quali.sup = 1, method = "Burt")

MCA_bladelet.morpho.variable.category <- fviz_mca_var(MCA_bladelet.morpho.analysis, labelsize = 4, col.var = "contrib",
                                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                      repel = TRUE, invisible = "quali.sup",
                                                      ggtheme = theme_gray(), title = "") +
  theme_minimal() +
  labs(y= "Dimension 2 (12.7%)", x = "Dimension 1 (25.7%)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), text = element_text(size = 13))

MCA_bladelet.morpho.variable.category

# Scree plot #
MCA_bladelet.morpho.screeplot <- fviz_eig(MCA_bladelet.morpho.analysis, addlabels=T, labelsize=2,
                                          ggtheme=theme_grey(base_size = 11)) +ylim(0,26) +
  ggtitle("")

ggsave("./output/figures/Figure_S29.png", plot = MCA_bladelet.morpho.screeplot, width = 5, height = 5, units = "in", dpi = 300)

# Contribution of variables #
MCA_bladelet.morpho.contribution <- fviz_contrib(MCA_bladelet.morpho.analysis, choice = "var", axes = 1:2, top = 8,
                                                 ggtheme=theme_grey(base_size = 11)) +
  ggtitle("")

ggsave("./output/figures/Figure_S30.png", plot = MCA_bladelet.morpho.contribution, width = 5, height = 5, units = "in", dpi = 300)

# ANOVA, data preparation after Cascalheira (2019) ####

MCA_bladelet.morpho.analysis.coord <- as.data.frame(MCA_bladelet.morpho.analysis$ind)
context <- as.data.frame(MCA_bladelet.morpho$SU)

MCA_bladelet.morpho.analysis.coord <- MCA_bladelet.morpho.analysis.coord %>%
  dplyr::select(1:2)

MCA_bladelet.morpho.analysis.coord <- bind_cols(context, MCA_bladelet.morpho.analysis.coord)

MCA_bladelet.morpho.analysis.coord <- MCA_bladelet.morpho.analysis.coord %>%
  mutate(context = `MCA_bladelet.morpho$SU`, "dim1" = coord.Dim.1, "dim2" = coord.Dim.2) %>%
  dplyr::select(context, dim1, dim2)

# Multi-comparison bootstrap
boot_dim1.bladelet.morpho <- WRS2::mcppb20(dim1 ~ context, data = MCA_bladelet.morpho.analysis.coord, nboot=5000)
boot_dim2.bladelet.morpho <- WRS2::mcppb20(dim2 ~ context, data = MCA_bladelet.morpho.analysis.coord, nboot=5000)

# Convert results to data frames
multi_comp_dim1.bladelet.morpho <- as.data.frame(boot_dim1.bladelet.morpho$comp)
multi_comp_dim2.bladelet.morpho <- as.data.frame(boot_dim2.bladelet.morpho$comp)

# Add contexts names and transform to column "context"
multi_comp_dim1.bladelet.morpho$context <- comparison.MCA
multi_comp_dim1.bladelet.morpho <- multi_comp_dim1.bladelet.morpho[ !duplicated(names(multi_comp_dim1.bladelet.morpho)) ]
multi_comp_dim2.bladelet.morpho$context <- comparison.MCA
multi_comp_dim2.bladelet.morpho <- multi_comp_dim2.bladelet.morpho[ !duplicated(names(multi_comp_dim2.bladelet.morpho)) ]

multi_comp.bladelet.morpho <- full_join(multi_comp_dim1.bladelet.morpho, multi_comp_dim2.bladelet.morpho, by = "context", suffix = c("dim1", "dim2"))

multi_comp.bladelet.morpho <- select(multi_comp.bladelet.morpho, context, `p-valuedim1`, `p-valuedim2`)

# Build list with pairwise p-values
labels_dim1.bladelet.morpho <- setNames(multi_comp.bladelet.morpho$`p-valuedim1`, multi_comp.bladelet.morpho$context)
labels_dim2.bladelet.morpho <- setNames(multi_comp.bladelet.morpho$`p-valuedim2`, multi_comp.bladelet.morpho$context)

# Attribute letters to each grouping based on alpha = .05
grouping_dim1.bladelet.morpho <- data.frame(multcompView::multcompLetters(labels_dim1.bladelet.morpho)["Letters"])
grouping_dim1.bladelet.morpho <- rownames_to_column(grouping_dim1.bladelet.morpho, "context")

grouping_dim2.bladelet.morpho <- data.frame(multcompView::multcompLetters(labels_dim2.bladelet.morpho)["Letters"])
grouping_dim2.bladelet.morpho <- rownames_to_column(grouping_dim2.bladelet.morpho, "context")

grouping.bladelet.morpho <- full_join(grouping_dim1.bladelet.morpho, grouping_dim2.bladelet.morpho, by = "context", suffix = c("dim1", "dim2"))

grouping.bladelet.morpho <- mutate(grouping.bladelet.morpho, Lettersdim1 = toupper(Lettersdim1), Lettersdim2 = toupper(Lettersdim2))


# Plot trimmed mean and CI #
mca.bladelet.morpho.summary <- MCA_bladelet.morpho.analysis.coord %>% 
  dplyr::group_by(context) %>% 
  dplyr::summarise(meandim1 = mean(dim1, trim = 0.2), meandim2 = mean(dim2, trim = 0.2), 
                   LCIdim1 = MeanCI(dim1, trim = 0.2)[2], MCIdim1 = MeanCI(dim1, trim = 0.2)[3],
                   LCIdim2 = MeanCI(dim2, trim = 0.2)[2], MCIdim2 = MeanCI(dim2, trim = 0.2)[3]) 

m1.bladelet.morpho <- max(mca.bladelet.morpho.summary$MCIdim1)
m2.bladelet.morpho <- max(mca.bladelet.morpho.summary$MCIdim2)


dim1.bladelet.morpho <- mca.bladelet.morpho.summary %>%
  mutate(context = fct_relevel(context, "A2int", "A2ext", "A1ext", "A1east", "D3b alpha")) %>%
  ggplot() +
  geom_point(stat = "identity", aes(x = context, y = meandim1), size = 3) +
  geom_errorbar(stat = "identity", aes(x = context, ymin = LCIdim1, ymax = MCIdim1), width = 0) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed", show.legend = FALSE) +
  theme_minimal() +
  geom_text(data = grouping.bladelet.morpho, aes(x = context, y = m1.bladelet.morpho + 0.08, label = Lettersdim1, size = 15), show.legend = FALSE) +
  xlab("") +
  ylab("Dimension 1") +
  scale_colour_discrete(l = 40) +
  theme(text = element_text(size = 14)) +
  coord_flip()

dim1.bladelet.morpho

dim2.bladelet.morpho <- mca.bladelet.morpho.summary %>%
  mutate(context = fct_relevel(context, "A2int", "A2ext", "A1ext", "A1east", "D3b alpha")) %>%
  ggplot() +
  geom_point(stat = "identity", aes(x = context, y = meandim2), size = 3) +
  geom_errorbar(stat = "identity", aes(x = context, ymin = LCIdim2, ymax = MCIdim2), width = 0) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed", show.legend = FALSE) +
  theme_minimal() +
  geom_text(data = grouping.bladelet.morpho, aes(x = context, y = m2.bladelet.morpho + 0.08, label = Lettersdim2, size = 15), show.legend = FALSE) +
  xlab("") +
  ylab("Dimension 2") +
  scale_colour_discrete(l = 40) +
  theme(text = element_text(size = 14)) +
  coord_flip()

dim2.bladelet.morpho

dim1.dim2.bladelet.morpho <- gridExtra::grid.arrange(dim1.bladelet.morpho, dim2.bladelet.morpho, nrow = 1)

var.dim1.dim2.bladelet.morpho <- ggarrange (MCA_bladelet.morpho.variable.category, dim1.dim2.bladelet.morpho, labels = c("a", "b"), nrow = 2)

ggsave("./output/figures/Figure_12.pdf", plot = var.dim1.dim2.bladelet.morpho, width = 8, height = 8, units = "in", dpi = 300)
