
rm(list = ls())

#library(limma)
library(readxl)
library(compareGroups)
library(ggbiplot)
library(stringr)
# library(LumiR)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(broom)
library(data.table)
library(purrr)  # For map and reduce functions
library(ggsignif)  # Load the ggsignif package
library(ggpubr)
library(rstatix)
library("RColorBrewer")
### Define the directory where your data files are located

setwd("/home/vivischuch/Documents/luminex_data/data/")

options(stringsAsFactors = F)

### Load necessary libraries
# You might need to install these libraries using install.packages() if you don't have them installed already
library(DataExplorer)
### Read the XLSX file
data <- read.csv("Results-Luminex.csv")
colnames(data)

data <- data %>%
  mutate(across(`GM.CSF`:`TNF.α`, ~as.numeric(gsub(",", "", .x))))

data <- filter(data, Rep != "AVERAGE")

str(data)

#Exploratory Data Analysis
create_report(data)

### Impute data, if necessary, using pmm method
imp <- mice::mice(data = data[, -c(1:2)], method = "pmm")
data[, -c(1:2)] <- mice::complete(imp)


#Exploratory Data Analysis
create_report(data)


# Add an identifier for each pair of replicates
data <- data %>%
  mutate(PairID = rep(1:(n()/2), each = 2))

# Calculate the average for each pair of replicates
average_data <- data %>%
  group_by(Sample, PairID) %>%
  summarise(across(`GM.CSF`:`TNF.α`, ~mean(as.numeric(.x), na.rm = TRUE)), .groups = 'drop')

# Remove PairID if it's no longer needed, keeping the averages per sample
average_data <- select(average_data, -PairID)


write.csv(data, "imputeddata_11mar24.csv", row.names = F)

metadata <- read.csv("lumi_metadata.csv", header = TRUE)
colnames(metadata)[1] <- "Sample"
data_raw_plot1 <- merge(average_data, metadata)
data_raw_plot1$groups <- paste0(data_raw_plot1$Treatment, "_", data_raw_plot1$Gestation)
data_raw_plot2 <- data_raw_plot1[, c(1, 17:20, 2:16)]
names(data_raw_plot2) <- gsub(".", "_", names(data_raw_plot2), fixed = TRUE)

data_raw_plot2 <- data_raw_plot2 %>%
  filter(!str_detect(Treatment, "rigi ag"))


### Change variables and create tables to compare de groups of interest
str(data_raw_plot2)


# Continue from the previous step where you've created Treatment_Gestation
data_raw_plot2 <- data_raw_plot2 %>%
  mutate(Treatment_Gestation = paste(Treatment, Gestation, sep="."))

data_raw_plot2$Treatment_Gestation <- as.character(data_raw_plot2$Treatment_Gestation)



library(dplyr)
library(tibble)

All_cytokine_data <- data_raw_plot2
All_cytokine_data <- as.data.frame(All_cytokine_data)
# rownames(All_cytokine_data) <- All_cytokine_data$...1
# All_cytokine_data <- All_cytokine_data[,-1]

# PCA with variables and different colors for groups (split by treatments)
library(ggfortify)
df <- All_cytokine_data[,c(6:20)]

png("PCA_biplot.png", width = 600, height = 400)
autoplot(prcomp(df, center = T, scale. = T), data = All_cytokine_data, loadings = T,
         loadings.label = T, loadings.label.repel = T, loadings.label.size = 3,
         loadings.label.alpha = 0.8, loadings.colour = "grey") +
  geom_point(size = 2,aes(color =  Gestation)) +
  theme_minimal()
dev.off()



# ##############################

head(data_raw_plot2)


names(data_raw_plot2) <- gsub(".", "_", names(data_raw_plot2), fixed = TRUE)

data_raw_plot2 <- merge(data_raw_plot2, metadata)

head(data_raw_plot2)

original_df <- data_raw_plot2
original_df$Sample <- NULL
original_df$Donor <- NULL
original_df$groups <- NULL
original_df$Treatment_Gestation <- NULL

# Convert the dataframe to a data.table
setDT(original_df)

head(original_df)
library(ggpubr)

# Function to split data.table into a list of dataframes
split_dataframe <- function(dt, col_indices) {
  num_cols <- ncol(dt)
  # Create a list of dataframes, each containing the selected columns
  df_list <- lapply(seq(3, num_cols, by = 1), function(i) {
    selected_cols <- c(1, 2, i)  # Select first two columns and the current index
    dt[, ..selected_cols, with = FALSE]
  })
  
  return(df_list)
}

# apply the function split_dataframe
df_list <- split_dataframe(original_df)
names(df_list) <- colnames(original_df[,c(3:17)])


# Function to do the boxplots
boxplots <- sapply(df_list, simplify = F, USE.NAMES = T, function(x){
  df <- x
  df <- df_list[[10]]
  title <- colnames(df)[3]
  colnames(df)[3] <- "value"
  colnames(df)[1] <- "treat"
  colnames(df)[2] <- "group"
  df <- as.data.frame(df)
  df<-df[!(df$treat=="rigi ag" | df$treat=="PIC" | df$treat=="LPS"| df$treat=="ZIKV" | df$treat=="IFNa"| df$treat=="IFNB" | df$treat=="IFNL"),]
  
  df$treat <- factor(df$treat, levels = c("NT", "CMV", "HIV"))
  
  stat.test <- df %>%
    group_by(treat) %>%
    t_test(value ~ group) %>%
    add_significance("p")
  stat.test
  
  my_palette <- c("coral2", "steelblue", "gold", "darkgreen", "purple", "orange", "pink")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    df, x = names(df)[1], 
    y = names(df)[3], 
    fill = names(df)[2], 
    palette = my_palette,
    size = 0.2,
    legend = "top",
    legend.title = "",
    bxp.errorbar = TRUE,
    bxp.errorbar.width = 0.15, 
    outlier.shape = 1,
    xlab = "", 
    ylab = "pg/ml",
    title = title
  ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  bxp
  stat.test <- stat.test %>% add_xy_position(x = "treat", dodge = 0.7, fun = "max", step.increase = 0.12)
  bxp.complex <- bxp + stat_pvalue_manual(
    stat.test,  label = "{p.signif}", hide.ns = TRUE, tip.length = 0.01
  ) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
  bxp.complex

  return(bxp.complex)
})


plots <- do.call("ggarrange", c(boxplots, ncol=5, nrow =3, common.legend = TRUE, legend="top"))

png("boxplots_all.png", width = 1920, height = 1080)
plots
dev.off()


lapply(names(boxplots), 
       function(x)ggsave(filename=paste(x,".png",sep=""), plot=boxplots[[x]]))


# #################

library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(data.table)
library(ggpubr)
library(tidyr)

# Convert the dataframe to a data.table
setDT(original_df)

# Create a 'Treatment_Gestation' interaction column
original_df[, Interaction := paste(Treatment, Gestation, sep = ".")]

# Specifying treatments to be removed
treatments_to_remove <- c("IFNB", "IFNL", "LPS", "PIC", "IFNa", "ZIKV")

# Filtering original_df to remove specified treatments
original_df <- original_df[!original_df$Treatment %in% treatments_to_remove, ]
original_df$Interaction <- as.factor(original_df$Interaction)

# Convert back to dataframe for easier manipulation with dplyr and ggplot2
df <- as.data.frame(original_df)

# Long format for Kruskal-Wallis test
df_long <- df %>%
  pivot_longer(
    cols = !c(Treatment, Gestation, Interaction), # Exclude these columns
    names_to = "Cytokine", 
    values_to = "Level"
  )

# grouped_data <- df_long %>%
#   group_by(Treatment, Gestation)

# Display the structure of the long format dataframe
head(df_long)
str(df_long)

df_long <- df_long %>%
  mutate(
    Interaction = gsub("CMV\\.Early Gestation", "CMV EG", Interaction),
    Interaction = gsub("CMV\\.Term", "CMV TG", Interaction),
    Interaction = gsub("HIV\\.Early Gestation", "HIV EG", Interaction),
    Interaction = gsub("HIV\\.Term", "HIV TG", Interaction),
    Interaction = gsub("NT\\.Early Gestation", "NT EG", Interaction),
    Interaction = gsub("NT\\.Term", "NT TG", Interaction),
    Cytokine = gsub("_", "-", Cytokine)  # Remove underline from Cytokine column
  ) %>%
  mutate(Interaction = factor(Interaction, levels = c("NT EG", "HIV EG", "CMV EG", "NT TG", "HIV TG", "CMV TG")))

kruskal_results <- df_long %>%
  group_by(Cytokine) %>%
  kruskal_test(Level ~ Interaction) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>% # Adjust for multiple comparisons using Benjamini-Hochberg method
  add_significance("p") # Add significance column based on adjusted p-values


# Group the data by Treatment, Gestation, and Cytokine, and apply Shapiro-Wilk test
normality_test_results <- df_long %>%
  group_by(Treatment, Gestation, Cytokine) %>%
  summarise(p_value_shapiro = shapiro.test(Level)$p.value) %>%
  ungroup()

# View the results
print(normality_test_results)

# View the results
print(kruskal_results)

library(ggplot2)
library(ggsignif)
library(dunn.test)
library(dplyr)

plot_cytokine_with_significance <- function(cytokine_name, df_long, all_significant_comparisons) {
  df_cytokine <- filter(df_long, Cytokine == cytokine_name)
  
  df_cytokine$Interaction <- factor(df_cytokine$Interaction, 
                                    levels = c("HIV EG", "CMV EG", "NT EG", "HIV TG", "CMV TG", "NT TG"))
  
  post_hoc_results <- dunn.test(x = df_cytokine$Level, g = df_cytokine$Interaction, method = "bonferroni")
  
  significant_comparisons <- all_significant_comparisons[all_significant_comparisons$Cytokine == cytokine_name, ]
  
  p <- ggplot(df_cytokine, aes(x = Interaction, y = Level, fill = Gestation)) +
    geom_boxplot() +
    theme_classic() +  # Change to a minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 18)) +
    labs(title = cytokine_name, y = "pg/ml", x = "") +
    scale_fill_manual(values = c("#F8766D", "#00BFC4"))
  # scale_fill_brewer(palette = "Set2")
  
  if (nrow(significant_comparisons) > 0) {
    max_y <- max(df_cytokine$Level, na.rm = TRUE)
    y_increment <- max_y * 0.1  # Adjust as needed
    
    for (i in 1:nrow(significant_comparisons)) {
      comp <- strsplit(as.character(significant_comparisons$comparisons[i]), " - ")[[1]]
      if (length(comp) == 2) {  # Ensure we have a pair to compare
        p_value <- as.numeric(significant_comparisons[i, "P.adjusted"])
        significance <- ifelse(p_value > 0.05 & p_value < 0.07, "*", 
                               ifelse(p_value > 0.01 & p_value <= 0.05, "**", 
                                      ifelse(p_value <= 0.01, "***", "")))
        p <- p + geom_signif(comparison = list(comp),
                             map_signif_level = TRUE,
                             textsize = 6,
                             vjust = 0.5,
                             y_position = max_y + i * y_increment,
                             annotations = significance)
      }
    }
  }
  
  return(p)
}

unique_cytokines <- unique(df_long$Cytokine)

all_significant_comparisons <- data.frame()

for (cytokine in unique_cytokines) {
  df_cytokine <- filter(df_long, Cytokine == cytokine)
  post_hoc_results <- dunn.test(x = df_cytokine$Level, g = df_cytokine$Interaction, method = "bonferroni")
  significant_comparisons <- data.frame(Cytokine = cytokine,
                                        comparisons = post_hoc_results$comparisons,
                                        P.adjusted = post_hoc_results$P.adjusted)
  significant_comparisons <- significant_comparisons[significant_comparisons$P.adjusted < 0.07, ]
  all_significant_comparisons <- rbind(all_significant_comparisons, significant_comparisons)
}

# Print the table of all significant comparisons
print(all_significant_comparisons)



# Apply the function to all proteins
plots <- list()


for (cytokine in unique_cytokines) {
  plots[[cytokine]] <- plot_cytokine_with_significance(cytokine, df_long, all_significant_comparisons)
}

# Example: display or save the plot for a specific protein
# print(plots[["IL_5"]])  # Replace "GM_CSF" with the actual name of the protein

# Optionally, save plots to files
for (cytokine in names(plots)) {
  ggsave(paste0(cytokine, ".svg"), plots[[cytokine]], width = 10, height = 6)
}

library(gridExtra)

# List to store plots
plot_list <- list()

# Create and store plots
for (cytokine in names(plots)) {
  plot_list[[cytokine]] <- plots[[cytokine]] + theme(legend.position = "none")  # Remove legend from individual plots
}

# Arrange plots in a grid
grid_arranged <- do.call(grid.arrange, c(plot_list, ncol = 5))


# Save the grid of plots with a single legend
ggsave("grid_of_plots_with_legend.png", grid_arranged, width = 16, height = 12)

ggsave("grid_of_plots_with_legend.svg", grid_arranged, width = 16, height = 12)

