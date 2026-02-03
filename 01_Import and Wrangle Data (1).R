# Load necessary libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)    # For pivot_longer function
library(gridExtra) # For arranging plots
library(stringr)


# Load the dataset
Data_1 <- fread("B1_ST.csv", encoding = "Latin-1")
Data_2 <- fread("B2_ST.csv", encoding = "Latin-1")
Data_3 <- fread("B3_ST.csv", encoding = "Latin-1")


# Bind rows, keeping only common columns
Data <- bind_rows(Data_1, Data_2,Data_3)

# Remove columns that are all NA
Data <- Data %>% select(where(~ !all(is.na(.))))

Data <- Data %>%
  mutate(
    Slide = str_extract(`Image Location`, "(?<=CellDIVE_)\\d+"),
    Region = str_extract(`Image Location`, "R\\d+")
  )

#Remove Cytoplasm
Data <- Data %>% select(-contains("cytoplasm"))

#Remove Nuclei
Data <- Data %>% select(-contains("nuclei"))

#Meta Data 
Meta_Data <- Data[,c(2,18,33:38)]

#Classifier Areas
Classifier_Areas <- Data %>% select(contains("mm²")) %>% bind_cols(Meta_Data)

#Remove Area Vars
Data <- Data %>% select(-contains("mm²"))

#Intensity Data
Intensity_Data <- Data %>% select(contains("intensity")) %>% bind_cols(Meta_Data)

#Remove Intensity Vars
Data <- Data %>% select(-contains("intensity")) 

# Select Morphology Data including Area, Roundness, and Perimeter, then combine with Meta Data
Morphology_Data <- Data %>%
  select(contains("Area") | contains("Roundness") | contains("Perimeter")) %>%
  bind_cols(Meta_Data) 

# Remove Morphology Variables (Area, Roundness, and Perimeter)
Data <- Data %>%
  select(-contains("Area"), -contains("Roundness"), -contains("Perimeter"))

#Positive Cell Data
Positive_Cell_Data <- Data %>% select(contains("positive")) %>% bind_cols(Meta_Data)

#Remove Positive  Vars
Data <- Data %>% select(-contains("positive"))


# Create the combined table with all entries
clinical_class_lookup <- data.table(
  Slide = c("92865", "92866", "92867", "92868", "92869", "92870", 
            "92871", "92872", "92873", "92874", "92875", "92876", 
            "92877", "92878", "92879", "92880", "92881", "92882", 
            "92883", "95459", "95460", "95461", "95462", "95463", 
            "95464", "95465",  "95468", "95469", 
            "95470",  "95472", "98087", "98088", "98089", 
            "98090", "98091", "98092", "98093", "98094", "98095", 
            "98096", "98097", "98098", "98099", "98100", "98101", 
            "98102", "98103", "98104"),
  Clinical_Class = c("Non-dysplastic (Control)", "Non-dysplastic (Control)", "Mild OED", "Mild OED", "Mild OED", "Mild OED", 
                     "Moderate OED", "Moderate OED", "Moderate OED", "Moderate OED", 
                     "Moderate OED", "Moderate OED", "Severe OED", "Severe OED", 
                     "Severe OED", "Severe OED", "Severe OED", "OSCC", "OSCC", 
                     "Veruccous OED", "Veruccous OED", "Veruccous OED", "HPV OED", 
                     "HPV OED", "HPV OED", "HPV OED", "Mild OED", "Mild OED", 
                     "Mild OED", "Moderate OED", 
                     "Moderate OED", "Moderate OED", "Moderate OED", "Moderate OED", 
                     "Moderate OED", "Severe OED", "Severe OED", "Severe OED", 
                     "Severe OED", "Severe OED", "Severe OED","Severe OED", "Hyperkeratosis", 
                     "Hyperkeratosis", "Hyperkeratosis", 
                     "Mild OED", "Mild OED", "Mild OED"),
  Transformation = c("Normal Control","Normal Control", 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
                     0, 0, "Cancer Control", "Cancer Control", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 1, 1, 0, 0, 0, 1, 1, 1,0, 0, 0, 0, 
                     0, 0, 0,0,0)
)

# Function to merge clinical_class_lookup with a data frame
add_clinical_class <- function(df) {
  merged_df <- merge(df, clinical_class_lookup, by = "Slide", all.x = TRUE)
  return(merged_df)
}

# Apply the function to all data frames - if needed
Classifier_Areas <- add_clinical_class(Classifier_Areas)
Intensity_Data <- add_clinical_class(Intensity_Data)
Morphology_Data <- add_clinical_class(Morphology_Data)
Positive_Cell_Data <- add_clinical_class(Positive_Cell_Data)
Data <- add_clinical_class(Data)  

# Clan Intensity Data
Intensity_Data <- Intensity_Data [,-c(2,8)] #remove unwanted cols

# Categories to exclude - don't remove if focusing on transformation status
#Intensity_Data <-  Intensity_Data %>% filter (!Clinical_Class %in% c("Hyperkeratosis", "HPV OED", "Veruccous OED"))

# Calculate the centroid of each cell
Intensity_Data <- Intensity_Data %>%
  mutate(
    XCentroid = (XMin + XMax) / 2,
    YCentroid = -((YMin + YMax)) / 2
  )

Intensity_Data <- Intensity_Data %>%
  filter(
    (`Analysis Region` == "Epi" & `Classifier Label` == "Epi") |
      (`Analysis Region` %in% c("Stroma", "Sub001") & `Classifier Label` == "Stroma")
  )

Positive_Cell_Data <- Positive_Cell_Data %>%
  filter(
    (`Analysis Region` == "Epi" & `Classifier Label` == "Epi") |
      (`Analysis Region` %in% c("Stroma", "Sub001") & `Classifier Label` == "Stroma")
  )


Meta_Data <- Meta_Data %>%
  filter(
    (`Analysis Region` == "Epi" & `Classifier Label` == "Epi") |
      (`Analysis Region` %in% c("Stroma", "Sub001") & `Classifier Label` == "Stroma")
  )


write.csv(Intensity_Data, "Tidy_Intensity_Data.csv")

library(tibble)  # for add_row()

# Count unique slides per Clinical_Class
Clinical_Class_Counts <- Intensity_Data %>%
  distinct(Slide, Clinical_Class) %>%
  count(Clinical_Class)

# Add total row
Clinical_Class_Counts <- Clinical_Class_Counts %>%
  add_row(Clinical_Class = "Total", n = sum(.$n))

write.csv(Clinical_Class_Counts, "Clinical_Class_Counts.csv", row.names = FALSE)

# Count unique slides per Transformation
Transformation_Counts <- Intensity_Data %>%
  distinct(Slide, Transformation) %>%
  count(Transformation)

# Add total row
Transformation_Counts <- Transformation_Counts %>%
  add_row(Transformation = "Total", n = sum(.$n))

write.csv(Transformation_Counts, "Transformation_Counts.csv", row.names = FALSE)


# Save processed datasets for later use
write.csv(Classifier_Areas, "Classifier_Areas.csv", row.names = FALSE)
write.csv(Intensity_Data, "Intensity_Data.csv", row.names = FALSE)
write.csv(Morphology_Data, "Morphology_Data.csv", row.names = FALSE)
write.csv(Positive_Cell_Data, "Positive_Cell_Data.csv", row.names = FALSE)
write.csv(Data, "Data.csv", row.names = FALSE)
