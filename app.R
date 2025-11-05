# app.R ---------------------------------------------------------------

library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(cluster)   # daisy (Gower)
library(ape)       # pcoa
library(DT)
library(tidyr)

## ---- Load bird_traits from CSV -------------------------------------

bird_traits <- read.csv("bird_traits.csv", stringsAsFactors = FALSE)

# Make sure key numeric traits are numeric (just in case CSV did something goofy)
num_cols <- c(
  "female_min_mass", "female_max_mass",
  "male_min_mass",   "male_max_mass",
  "unsexed_min_mass","unsexed_max_mass",
  "avg_mass",
  "elev_xmin", "elev_norm_min", "elev_range",
  "elev_norm_max", "elev_xmax",
  "habitat_breadth", "diet_breadth", "esi"
)

bird_traits <- bird_traits %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

## ---- Shiny UI ------------------------------------------------------

ui <- fluidPage(
  titlePanel("BIRDBASE trait-space explorer (PCoA on Gower)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Filters"),
      
      selectInput(
        "order_filter",
        "Limit to order:",
        choices = c("All", sort(unique(na.omit(bird_traits$order)))),
        selected = "All"
      ),
      
      selectInput(
        "family_filter",
        "Limit to family (Clements):",
        choices = c("All", sort(unique(na.omit(bird_traits$family_clem)))),
        selected = "All"
      ),
      
      selectInput(
        "iucn_filter",
        "Limit to IUCN category:",
        choices = c("All", sort(unique(na.omit(bird_traits$iucn2024)))),
        selected = "All"
      ),
      
      pickerInput(
        "species_select",
        "Select focal species (optional):",
        choices = sort(unique(bird_traits$species_label)),
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE
        )
      ),
      
      hr(),
      h4("Traits for distance"),
      
      pickerInput(
        "trait_select",
        "Traits to use in Gower distance:",
        choices = c(
          # continuous-ish
          "avg_mass",
          "elev_range",
          "habitat_breadth",
          "diet_breadth",
          "esi",
          # categorical / codes
          "primary_habitat",
          "primary_diet",
          "iucn2024",
          "flightlessness",
          "mig",
          "alt_move",
          "irreg_move",
          "dispersive",
          "sedentary",
          "order",
          "family_clem"
        ),
        selected = c("avg_mass", "elev_range", "habitat_breadth", "diet_breadth", "esi"),
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE
        )
      ),
      
      numericInput(
        "min_species",
        "Minimum species to run PCoA:",
        value = 20,
        min = 5,
        max = 500
      ),
      
      numericInput(
        "max_species",
        "Maximum species (random subset if exceeded):",
        value = 500,
        min = 50,
        max = 3000
      ),
      
      selectInput(
        "colour_by",
        "Colour points by:",
        choices = c(
          "primary_habitat",
          "primary_diet",
          "iucn2024",
          "order",
          "family_clem",
          "mig"
        ),
        selected = "primary_habitat"
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Ordination",
                 br(),
                 textOutput("summary_text"),
                 plotOutput("pcoa_plot", height = "600px",
                            hover = "plot_hover"),
                 verbatimTextOutput("hover_info")
        ),
        tabPanel("Trait distributions",
                 br(),
                 plotOutput("trait_hist", height = "450px")
        ),
        tabPanel("Data",
                 br(),
                 DTOutput("trait_table")
        ),
        tabPanel("About / Citation",
                 br(),
                 h3("BIRDBASE dataset"),
                 p("This app uses the BIRDBASE dataset as described in:"),
                 tags$p(
                   tags$em(
                     "Şekercioğlu, Ç.H., Kittelberger, K.D., Mota, F.M.M. et al. ",
                     "BIRDBASE: A Global Dataset of Avian Biogeography, ",
                     "Conservation, Ecology and Life History Traits. ",
                     "Sci Data 12, 1558 (2025). https://doi.org/10.1038/s41597-025-05615-3"
                   )
                 ),
                 p("Please cite this paper if you use insights from this app in any publication or report.")
        )
      )
    )
  )
)

## ---- Shiny server --------------------------------------------------

server <- function(input, output, session) {
  
  # 1) Filter species pool based on order/family/IUCN and optional species selection
  filtered_species <- reactive({
    dat <- bird_traits
    
    if (input$order_filter != "All") {
      dat <- dat %>% filter(order == input$order_filter)
    }
    if (input$family_filter != "All") {
      dat <- dat %>% filter(family_clem == input$family_filter)
    }
    if (input$iucn_filter != "All") {
      dat <- dat %>% filter(iucn2024 == input$iucn_filter)
    }
    if (length(input$species_select) > 0) {
      dat <- dat %>% filter(species_label %in% input$species_select)
    }
    
    # Random subset if too many species (performance guardrail)
    if (nrow(dat) > input$max_species) {
      set.seed(123)
      dat <- dat %>% slice_sample(n = input$max_species)
    }
    
    dat
  })
  
  # 2) Prepared trait data (drop NAs in selected traits)
  trait_data <- reactive({
    req(input$trait_select)
    dat <- filtered_species()
    
    cols_needed <- unique(c(
      "species_label",
      "order", "family_clem",
      "primary_habitat", "primary_diet",
      "iucn2024", "mig",
      input$trait_select
    ))
    
    dat <- dat %>% select(all_of(cols_needed))
    
    cc <- complete.cases(dat[, input$trait_select, drop = FALSE])
    dat <- dat[cc, ]
    
    dat
  })
  
  # 3) Run Gower + PCoA
  pcoa_res <- reactive({
    dat <- trait_data()
    
    if (nrow(dat) < input$min_species) {
      return(NULL)
    }
    
    # Trait-only matrix for distances
    trait_mat <- dat[, input$trait_select, drop = FALSE]
    
    # IMPORTANT: convert character columns to factors for Gower
    trait_mat[] <- lapply(trait_mat, function(col) {
      if (is.character(col)) factor(col) else col
    })
    
    # Now Gower distance (handles numeric + factor/ordered/logical)
    gower_dist <- daisy(trait_mat, metric = "gower")
    
    # PCoA
    pcoa_out <- ape::pcoa(gower_dist)
    
    scores <- as.data.frame(pcoa_out$vectors[, 1:2])
    colnames(scores) <- c("PCoA1", "PCoA2")
    
    scores$species_label   <- dat$species_label
    scores$order           <- dat$order
    scores$family_clem     <- dat$family_clem
    scores$primary_habitat <- dat$primary_habitat
    scores$primary_diet    <- dat$primary_diet
    scores$iucn2024        <- dat$iucn2024
    scores$mig             <- dat$mig
    
    list(scores = scores, eig = pcoa_out$values)
  })
  
  
  # 4) Summary text
  output$summary_text <- renderText({
    dat_total <- filtered_species()
    dat_used  <- trait_data()
    
    if (is.null(pcoa_res())) {
      paste0(
        "Species after filters: ", nrow(dat_total),
        " | Species with complete selected traits: ", nrow(dat_used),
        " (need at least ", input$min_species, " for PCoA)."
      )
    } else {
      eig <- pcoa_res()$eig
      var1 <- round(100 * eig$Relative_eig[1], 1)
      var2 <- round(100 * eig$Relative_eig[2], 1)
      paste0(
        "Species after filters: ", nrow(dat_total),
        " | Species in PCoA (complete traits): ", nrow(dat_used),
        " | PCoA1 explains ", var1, "%; PCoA2 explains ", var2,
        "% of variance in Gower distances."
      )
    }
  })
  
  # 5) Ordination plot
  output$pcoa_plot <- renderPlot({
    res <- pcoa_res()
    validate(
      need(!is.null(res),
           "Not enough species with complete trait data to run PCoA. Increase species pool or relax filters.")
    )
    
    scores <- res$scores
    colour_var <- input$colour_by
    scores$colour_value <- scores[[colour_var]]
    
    ggplot(scores, aes(x = PCoA1, y = PCoA2, colour = colour_value)) +
      geom_point(alpha = 0.8) +
      theme_bw() +
      labs(
        colour = colour_var,
        x = "PCoA1",
        y = "PCoA2"
      )
  })
  
  # 6) Hover info (nearest species)
  output$hover_info <- renderPrint({
    res <- pcoa_res()
    if (is.null(res) || is.null(input$plot_hover)) return(NULL)
    
    scores <- res$scores
    d2 <- (scores$PCoA1 - input$plot_hover$x)^2 +
      (scores$PCoA2 - input$plot_hover$y)^2
    i <- which.min(d2)
    
    scores[i, c("species_label", "order", "family_clem",
                "primary_habitat", "primary_diet", "iucn2024", "mig")]
  })
  
  # 7) Trait histograms (for numeric traits only)
  output$trait_hist <- renderPlot({
    dat <- trait_data()
    req(nrow(dat) > 0, length(input$trait_select) > 0)
    
    cont_traits <- input$trait_select[
      sapply(dat[, input$trait_select, drop = FALSE], is.numeric)
    ]
    
    validate(
      need(length(cont_traits) > 0,
           "No numeric traits selected; choose at least one numeric trait for histograms.")
    )
    
    traits_long <- dat %>%
      select(species_label, all_of(cont_traits)) %>%
      pivot_longer(-species_label,
                   names_to = "trait",
                   values_to = "value")
    
    ggplot(traits_long, aes(x = value)) +
      geom_histogram(bins = 30) +
      facet_wrap(~ trait, scales = "free") +
      theme_bw()
  })
  
  # 8) Data table
  output$trait_table <- renderDT({
    trait_data()
  })
}

shinyApp(ui, server)
