# cleaning the birdbase dataset into something that can be used 

library(readxl)
library(dplyr)
library(stringr)

# 1) Read the raw sheet WITHOUT treating the first row as header

raw <- read_excel(
  "BIRDBASE v2025.1 Sekercioglu et al. Final.xlsx",
  sheet = "Data",
  col_names = FALSE
)

#2) 
# Row 1 = big section headings (Species ID / Taxonomy / ...)
# Row 2 = actual variable names (IOC 15.1, English Name, Latin, ...)
header <- as.character(unlist(raw[2, ]))

# Drop first two rows from the data
bird_raw <- raw[-c(1, 2), ]

# Assign names from row 2
names(bird_raw) <- header





# 3) Rename key columns to names that make sense
bird_traits <- bird_raw %>%
  rename(
    species_ioc   = `IOC 15.1`,
    name_en       = `English Name (BirdLife > IOC > Clements>AviList)`,
    name_latin    = `Latin (BirdLife > IOC > Clements>AviList)`,
    order         = Order,
    family_ioc    = `Family IOC 15.1`,
    family_clem   = `Family Clements v2024b`,
    genus         = Genus,
    species_epithet = Species,
    iucn2024      = `2024 IUCN Red List category`,
    female_min_mass   = `Female MinMass`,
    female_max_mass   = `Female MaxMass`,
    male_min_mass     = `Male MinMass`,
    male_max_mass     = `Male MaxMass`,
    unsexed_min_mass  = `Unsexed MinMass`,
    unsexed_max_mass  = `Unsexed MaxMass`,
    avg_mass          = `Average Mass`,
    elev_xmin         = Xmin,
    elev_norm_min     = NormMin,
    elev_range        = `Elevational Range`,
    elev_norm_max     = NormMax,
    elev_xmax         = Xmax,
    primary_habitat   = `Primary Habitat`,
    habitat_breadth   = HB,
    primary_diet      = `Primary Diet`,
    diet_breadth      = DB,
    esi               = ESI,
    flightlessness    = Flightlessness,
    mig               = Mig,
    alt_move          = Alt,
    irreg_move        = Irreg,
    dispersive        = Disp,
    sedentary         = Sed
  ) %>%
  # optional: use Latin name as a clean species label
  mutate(
    species_label = if_else(
      !is.na(name_en),
      paste0(name_en, " (", name_latin, ")"),
      name_latin
    )
  )

# 4) Convert numeric-looking columns to numeric
num_cols <- c(
  "female_min_mass", "female_max_mass",
  "male_min_mass", "male_max_mass",
  "unsexed_min_mass", "unsexed_max_mass",
  "avg_mass",
  "elev_xmin", "elev_norm_min", "elev_range",
  "elev_norm_max", "elev_xmax",
  "habitat_breadth", "diet_breadth", "esi"
)

bird_traits <- bird_traits %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

write.csv(bird_traits, "bird_traits.csv")


