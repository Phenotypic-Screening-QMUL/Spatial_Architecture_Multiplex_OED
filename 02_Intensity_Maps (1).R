# Useful Packages
library("gplots")       # Graphical plots
library("RColorBrewer") # Color palettes
library("matrixStats")  # Efficient matrix computations
library("plyr")        # Data manipulation
library("dplyr")       # Data manipulation
library("data.table")  # High-performance data manipulation ``
library("stringr")     # String manipulation
library("ggplot2")     # Elegant plots
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(scales)
library(FNN)
library(tidyr)

Data_Intensity <- fread("Tidy_Intensity_Data.csv", encoding = "Latin-1")

Data_Intensity <- Data_Intensity %>%
  mutate(Slide = paste(Slide, `Region`, sep = "_"))

# Select classification columns
intensity_cols <- names(Data_Intensity)[c(3:7,9:12)]  # Extract column names, not data

#Example Slides
Slide_Examples <- Data_Intensity %>% filter(Slide %in% c("92865_R0" , "92882_R1"))

for (Slide in unique(Slide_Examples$Slide)) {
  
  # Ensure filtering is correct
  Slide_data <- Slide_Examples %>%
    filter(Slide == !!Slide)  # Use !! to correctly reference the loop variable
  
  # Loop through each intensity column
  for (col in intensity_cols) {
    
    # Generate plot for the current intensity column
    TMA_plot_colored <- ggplot(Slide_data, aes(x = XCentroid, y = YCentroid, color = .data[[col]])) +
      geom_point(size = 0.1) +
      scale_color_viridis(trans = "log10") +  
      labs(title = paste("Cell Map -", Slide, "Colored by", col),
           x = "X Centroid", y = "Y Centroid", color = col) +  
      theme_minimal() +  
      theme(plot.background = element_rect(fill = "white"),  
            panel.background = element_rect(fill = "white", color = NA),  
            legend.background = element_rect(fill = "white", color = NA),  
            legend.key = element_rect(fill = "white", color = NA),  
            axis.text = element_blank(),  
            axis.title = element_text(size = 12),  
            axis.ticks = element_blank()) +  
      scale_x_continuous(labels = NULL) +  
      scale_y_continuous(labels = NULL) 
    
    # Save plot
    ggsave(paste0("Intensity Cell Maps/", Slide, "_Cell_Map_", gsub(" ", "_", col), ".png"), TMA_plot_colored)
  }
  
  # Define custom colors for Analysis Region
  region_colors <- c("Epi" = "#008000", "Stroma" = "#E4D00A", "Sub001" = "#FF5733")
  
  # Modify the plot to use the custom colors
  TMA_plot_category <- ggplot(Slide_data, aes(x = XCentroid, y = YCentroid, color = `Analysis Region`)) +
    geom_point(size = 0.1) +
    scale_color_manual(values = region_colors) +  # Use custom colors
    labs(title = paste("Cell Map -", Slide, "Colored by Category"),
         x = "X Centroid", y = "Y Centroid", color = "Category") +  # Set axis titles
    theme_minimal() +  # Use minimal theme
    theme(plot.background = element_rect(fill = "white"),  
          panel.background = element_rect(fill = "white", color = NA),  
          legend.background = element_rect(fill = "white", color = NA),  
          legend.key = element_rect(fill = "white", color = NA),  
          axis.text = element_blank(),  
          axis.title = element_text(size = 12),  
          axis.ticks = element_blank()) +  
    scale_x_continuous(labels = NULL) +  
    scale_y_continuous(labels = NULL)
  
  # Save the plot
  ggsave(paste0("Intensity Cell Maps/", Slide, "_Cell_Map_Colored_by_Region.png"), TMA_plot_category)
  
}

# Scale Intensity Columns - This will be used in later clustering
Data_Intensity <- Data_Intensity %>%
  mutate(across(all_of(intensity_cols), ~ as.vector(scale(.))))

write.csv(Data_Intensity, "Data_Intensity_Scaled.csv")
