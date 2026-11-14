   #UPLOADING DATABASE ONTO MY WORKSPACE 

#Exploratory analysis of Antimicrobial Resistance Patterns in Bacterial Isolates Using R 
#packages for the study
#install.packages("tidyverse" , "janitor")
#library(tidyverse)
#install.packages("janitor")
#library(janitor)
#getwd()

# this is the data retrival stage
#amr <- read.csv("C:/Users/bisma/OneDrive/Documents/CPY_SAVES/AMR WITH R/simulated_amr_data.csv")
#view(amr)
#amr <- read.csv("C:/Users/bisma/OneDrive/Documents/CPY_SAVES/AMR WITH R/simulated_amr_data.csv",na = C("","NA"))|> clean_names()

#Lets look at the data 
head(amr)
glimpse(amr)
summary(amr)

#Count the number of the isolates that belong to each species
count(amr ,species)
amr |> count(species)

# an inspection summary of the  first part of the work 
writeLines(
  c(
    paste("number of isolates:",nrow(amr)),
    paste("number of columns:",ncol(amr)),
    paste("species:", paste(unique(amr$species), collapse = ";"))
  ),
  "results/data_overview.txt"
  ) 


#CLEANING AND PREPARING OF THE UPLOADED DATA
  #TO CHECK DATA QUALITY AND CONVERTING THE DATA FROM WIDE TO LONG FORMAT 
#library(tidyverse)
#library(janitor)
#COUNTING MISSING VALUES IN COLUMNS

missing_summary <- amr |> 
  summarise(across(everything(),~ sum(is.na(.)))) |> 
  pivot_longer(everything(), names_to = "variable", values_to = "missing_values")
print(missing_summary)
write_csv(missing_summary , "results/missing_summary.csv")

#CONVERT ANTIBIOTIC COLUMN INTO TWO COLUMNS: ANTIBIOTIC AND RESULTS

amr_long <- amr |> 
  pivot_longer(
    cols = all_of(c("AMP","CIP","GEN","CRO","TET","SXT")),
    names_to = "antibiotics",
    values_to = "results"
  ) |> 
  mutate(
    results =str_to_upper(str_trim(results)),
    results =na_if(results, "")
  ) |> 
  filter(!is.na(results) , results %in% c('S',"I","R"))


#COUNT THE DATASET WE JUST CREATED
head(amr_long)
summarise(amr_long)
glimpse(amr_long)
count(amr_long)

#RESISTANCE ANALYSIS
 #REISITANCE RATES AND MULTIDRUG RESISTANCE ANALYSIS 
# Resistance rates = number of R results / number of valid results * 100
resistance_summary <- amr_long |> 
  group_by(species, antibiotics ) |> 
  summarise(
    tested = n(),
    resistant = sum(results == "R"),
    resistance_percent = 100 * resistant /tested,
    .groups ="drop"
    ) |> 
  arrange(species, desc(resistance_percent))

print(resistance_summary)
write.csv(resistance_summary, "results/resistance_by_species_antibiotic.csv")

#MDR = multiple drug resistance 
isolate_mdr <- amr_long |>
  group_by(isolate_id, species) |> 
  summarise(
    antibiotics_tested = n(),
    antibiotics_resistant = sum(results== "R"),
    mdr = antibiotics_resistant >= 3,
    .groups = "drop"
     )

mdr_summary <-  isolate_mdr |> 
  group_by(species) |> 
  summarise(
    isolates = n(),
    mdr_isolates = sum(mdr),
    mdr_percent = 100 * mdr_isolates /isolates,
    .groups = "drop"
    )


write.csv(isolate_mdr, "results/mdr_by_isolates.csv")
write.csv(mdr_summary, "results/mdr_summary_by_species.csv")


# figures to visualize my results 
  #Figure 1 : Heatmap of resistance percentages by species and antibiotics 
heatmap <- ggplot(
  resistance_summary, 
  aes(x = antibiotics,y= species, fil = resistance_percent)
 ) +
  geom_tile(colour = "white") + 
  geom_text(aes(label = paste0(round(resistance_percent), "%")), size = 3) + 
  scale_fill_gradient(low = "white", high = "firebrick",
                      name = "Resistance (%)") +
  labs(
    title = "Antibiotic resistance patterns by bacterial species",
    x = "Antibiotic code ",
    y = "Bacterial species"
    ) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

eggsave ("figures/resistance_heatmap.png", heatmap,
         width = 9, height = 5, dpi = 300)


# Figure 2: Multidrug resistance percentage by species
mdr_plot <- ggplot(
  mdr_summary,
  aes(x = reorder(species, mdr_percent), y = mdr_percent)
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Proportion of isolates classified as MDR",
    x = "Bacterial species",
    y = "MDR isolates (%)"
  ) +
  theme_minimal()

ggsave("figures/mdr_by_species.png", mdr_plot,
       width = 8, height = 5, dpi = 300)




#Figure 3: Number of resistant results per antibiotic

amr_long <- amr |>
  pivot_longer(
    cols = all_of(c("AMP", "CIP", "CRO", "GEN", "TET", "SXT")),
    names_to = "antibiotic",
    values_to = "result"
  ) |>
  mutate(result = str_to_upper(str_trim(result))) |>
  filter(!is.na(result), result %in% c("S", "I", "R"))

antibiotic_plot_data <- amr_long |>
  group_by(antibiotic) |>
  summarise(resistant_results = sum(result == "R"), .groups = "drop")


antibiotic_plot <- ggplot(
  antibiotic_plot_data,
  aes(x = reorder(antibiotic, resistant_results), y = resistant_results)
) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(
    title = "Count of resistant test results by antibiotic",
    x = "Antibiotic code",
    y = "Number of resistant results"
  ) +
  theme_minimal()

ggsave("figures/resistant_results_by_antibiotic.png", antibiotic_plot,
       width = 7, height = 5, dpi = 300)

message("Figures saved in the figures/ folder.")