# Title     : TODO
# Objective : TODO
# Created by: chema
# Created on: 2023-10-26

# Load the UpSetR package
library(UpSetR)

# Read the CSV file into a data frame
gene_data <- read.csv("gene_data.csv")

# Create the UpSet plot
upset(gene_data, sets.bar.color="lightblue", order.by = "freq")
