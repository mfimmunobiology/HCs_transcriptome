
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
data <- read.csv("Results-Table 1.csv")
colnames(data)

data <- data %>%
  mutate(across(`GM.CSF`:`TNF.α`, ~as.numeric(gsub(",", "", .x))))

data <- filter(data, Rep != "AVERAGE")

str(data)

#Exploratory Data Analysis
# create_report(data)



### Evaluatew the batch effect (more than 1 plate)- if applicable
# To do this, you need to remove the categorical variables
# data2 = data[,-c(1:2)]
# data2[is.na(data2)] = 0
# input_data_pc  = prcomp(data2,
#                         center = T,
#                         scale. = T)
# ggbiplot(input_data_pc,
#          obs.scale = 1,
#          var.scale = 1,
#          groups = data2$Plate,
#          ellipse = T,
#          var.axes = F,
#          circle = T) +
#   theme(legend.direction = 'vertical',
#         legend.position = 'right')+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position="bottom",
#         aspect.ratio = 1)+
#   ggtitle("")



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


str(data_raw_plot2)


library(ggplot2)
library(dplyr)

# Assuming 'data_raw_plot2' is your dataframe and it includes the measure in a column named 'IL_6'.
# Creating a new column combining 'Treatment' and 'Gestation' if not already done
data_raw_plot2 <- data_raw_plot2 %>%
  mutate(Treatment_Gestation = paste(Treatment, Gestation, sep="."))

# Plotting IL-6 levels across these combined groups, with color based on Treatment
ggplot(data_raw_plot2, aes(x=Treatment_Gestation, y=IL_6, fill=Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title="IL-6 Levels by Treatment and Gestation", x="Treatment and Gestation", y="IL-6 (pg/ml)") +
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.title=element_blank()) +
  scale_fill_brewer(palette="Paired") # Colors the boxplots based on Treatment



# Continue from the previous step where you've created Treatment_Gestation
data_raw_plot2 <- data_raw_plot2 %>%
  mutate(Treatment_Gestation = paste(Treatment, Gestation, sep="."))

data_raw_plot2$Treatment_Gestation <- as.character(data_raw_plot2$Treatment_Gestation)

# Then list the unique values to verify your group names
unique(data_raw_plot2$Treatment_Gestation)
# Step 1: Extract unique treatments and gestations
unique_treatments <- unique(gsub("\\..*", "", data_raw_plot2$Treatment_Gestation))

# Step 2: Identify pairs for comparison
comparison_list <- lapply(unique_treatments, function(treatment) {
  early <- paste(treatment, "Early Gestation", sep=".")
  term <- paste(treatment, "Term", sep=".")
  if(early %in% data_raw_plot2$Treatment_Gestation & term %in% data_raw_plot2$Treatment_Gestation){
    return(c(early, term))
  }
})

# Remove NULL elements
comparison_list <- Filter(Negate(is.null), comparison_list)


ggplot(data_raw_plot2, aes(x=Treatment_Gestation, y=TNF_α, fill=Treatment)) +
  geom_boxplot() +
  geom_signif(comparisons = comparison_list, map_signif_level = TRUE) +
  theme_minimal() +
  labs(title="TNF_α Levels by Treatment and Gestation", x="Treatment and Gestation", y="TNF_α (pg/ml)") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.title=element_blank()) +
  scale_fill_brewer(palette="Paired")

ggsave(filename = "luminex_analysis_plot_TNF_α.svg", plot = last_plot())



#PCA ANALYSIS



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


# data_raw_plot <- as.data.frame(t(data_raw_plot2))
names(data_raw_plot2) <- gsub(".", "_", names(data_raw_plot2), fixed = TRUE)
# data_raw_plot$groups <- metadata$groups
# data_raw_plot2$Sample <- rownames(data_raw_plot2)
data_raw_plot2 <- merge(data_raw_plot2, metadata)

head(data_raw_plot2)

original_df <- data_raw_plot2
original_df$Sample <- NULL
original_df$Donor <- NULL
original_df$groups <- NULL
original_df$Treatment_Gestation <- NULL

# original_df[, 4:18] <- log(original_df[3:17], 10)

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

# Adjusted function to perform post-hoc tests and generate plots for all proteins
plot_cytokine_with_significance <- function(cytokine_name, df_long, all_significant_comparisons) {
  df_cytokine <- filter(df_long, Cytokine == cytokine_name)
  
  post_hoc_results <- dunn.test(x = df_cytokine$Level, g = df_cytokine$Interaction, method = "bonferroni")
  
  significant_comparisons <- all_significant_comparisons[all_significant_comparisons$Cytokine == cytokine_name, ]
  
  p <- ggplot(df_cytokine, aes(x = Interaction, y = Level, fill = Treatment)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = cytokine_name, y = "pg/ml", x = "") +
    scale_fill_brewer(palette = "Pastel1")
  
  if (nrow(significant_comparisons) > 0) {
    max_y <- max(df_cytokine$Level, na.rm = TRUE)
    y_increment <- max_y * 0.1  # Adjust as needed
    
    for (i in 1:nrow(significant_comparisons)) {
      comp <- strsplit(as.character(significant_comparisons$comparisons[i]), " - ")[[1]]
      if (length(comp) == 2) {  # Ensure we have a pair to compare
        p_value <- as.numeric(significant_comparisons[i, "P.adjusted"])
        significance <- ifelse(p_value > 0.07, "*", ifelse(p_value > 0.05, "**", ifelse(p_value > 0.03, "***", "***")))
        p <- p + geom_signif(comparison = list(comp),
                             map_signif_level = TRUE,
                             textsize = 4,
                             vjust = 0.5,
                             y_position = max_y + i * y_increment,
                             annotations = significance)
      }
    }
  }
  
  return(p)
}


unique_cytokines <- unique(df_long$Cytokine)

# Combine significant comparisons for all cytokines into one data frame
all_significant_comparisons <- data.frame()

for (cytokine in unique_cytokines) {
  df_cytokine <- filter(df_long, Cytokine == cytokine)
  post_hoc_results <- dunn.test(x = df_cytokine$Level, g = df_cytokine$Interaction, method = "bonferroni")
  significant_comparisons <- data.frame(Cytokine = cytokine,
                                        comparisons = post_hoc_results$comparisons,
                                        P.adjusted = post_hoc_results$P.adjusted)
  significant_comparisons <- significant_comparisons[significant_comparisons$P.adjusted < 0.1, ]
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
print(plots[["IL_5"]])  # Replace "GM_CSF" with the actual name of the protein

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
ggsave("grid_of_plots_with_legend.svg", grid_arranged, width = 15, height = 10)


