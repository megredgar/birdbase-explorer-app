# app.R ---------------------------------------------------------------
# BIRDBASE Trait-space Explorer 🐦 (PCoA on Gower, interactive)

library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(cluster)   # daisy (Gower)
library(ape)       # pcoa
library(DT)
library(tidyr)
library(plotly)
library(scales)    # for hue_pal()

## ---- Load bird_traits from CSV -------------------------------------

bird_traits <- read.csv("bird_traits.csv", stringsAsFactors = FALSE)

# Make sure key numeric traits are numeric
num_cols <- c(
  "female_min_mass", "female_max_mass",
  "male_min_mass",   "male_max_mass",
  "unsexed_min_mass","unsexed_max_mass",
  "avg_mass",
  "elev_xmin", "elev_norm_min", "elev_range",
  "elev_norm_max", "elev_xmax",
  "habitat_breadth", "diet_breadth", "esi",
  "Clutch_Min", "Clutch_Max"
)

bird_traits <- bird_traits %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

## ---- Shiny UI ------------------------------------------------------

ui <- fluidPage(
  # Bubblegum pink header
  tags$div(
    style = "background-color:#FF69B4; padding: 15px 20px; margin-bottom: 20px;",
    h2(
      "BIRDBASE Trait-space Explorer 🐦",
      style = "color:white; margin:0;"
    )
  ),
  
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
        "Select focal species (optional subset):",
        choices = sort(unique(bird_traits$species_label)),
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE
        )
      ),
      
      textInput(
        "search_species",
        "Highlight species by name (string match):",
        value = ""
      ),
      
      hr(),
      h4("Traits for distance"),
      
      selectInput(
        "trait_preset",
        "Trait preset:",
        choices = c(
          "Custom" = "custom",
          "Body size & elevation" = "size_elev",
          "Habitat & diet specialization" = "hab_diet",
          "Breeding ecology" = "breeding",
          "Movement & residency" = "movement",
          "Social behaviour" = "social"
        ),
        selected = "custom"
      ),
      
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
          "Clutch_Min",
          "Clutch_Max",
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
          "family_clem",
          "Nest_Type",
          "Nest_SBS"
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
      ),
      
      checkboxInput(
        "show_hulls",
        "Show group hulls (by colour variable)",
        value = FALSE
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        # 1) Ordination
        tabPanel("Ordination",
                 br(),
                 textOutput("summary_text"),
                 br(),
                 fluidRow(
                   column(4, downloadButton("download_scores", "Download PCoA scores")),
                   column(4, downloadButton("download_traits", "Download filtered traits"))
                 ),
                 br(),
                 plotlyOutput("pcoa_plot", height = "600px"),
                 br(),
                 h4("Trait–axis associations (correlation with PCoA axes)"),
                 DTOutput("axis_traits_table")
        ),
        
        # 2) Trait distributions
        tabPanel("Trait distributions",
                 br(),
                 plotlyOutput("trait_hist", height = "450px")
        ),
        # 3) Data
        tabPanel("Data",
                 br(),
                 DTOutput("trait_table")
        ),
        # 4) Nearest neighbours
        tabPanel("Nearest neighbours",
                 br(),
                 fluidRow(
                   column(
                     6,
                     selectizeInput(
                       "nn_species",
                       "Focal species:",
                       choices = NULL,
                       options = list(
                         placeholder = "Type to search species...",
                         maxOptions = 200
                       )
                     )
                   ),
                   column(
                     2,
                     numericInput(
                       "nn_k",
                       "Number of nearest neighbours:",
                       value = 10,
                       min = 1,
                       max = 100
                     )
                   )
                 ),
                 br(),
                 DTOutput("nn_table")
        ),
        # 5) Functional diversity
        tabPanel("Functional diversity",
                 br(),
                 selectInput(
                   "fd_group",
                   "Group by:",
                   choices = c(
                     "order",
                     "family_clem",
                     "primary_habitat",
                     "primary_diet",
                     "iucn2024",
                     "mig"
                   ),
                   selected = "order"
                 ),
                 br(),
                 DTOutput("fd_table")
        ),
        # 6) Compare subsets
        tabPanel("Compare subsets",
                 br(),
                 selectInput(
                   "cmp_group_var",
                   "Variable to compare:",
                   choices = c(
                     "primary_habitat",
                     "primary_diet",
                     "iucn2024",
                     "order",
                     "family_clem",
                     "flightlessness",
                     "mig"
                   ),
                   selected = "primary_habitat"
                 ),
                 uiOutput("cmp_level_ui"),
                 br(),
                 plotlyOutput("compare_plot", height = "600px"),
                 br(),
                 verbatimTextOutput("compare_summary")
        ),
        # 7) About / Citation
        tabPanel("About / Citation",
                 br(),
                 h3("What this app does"),
                 p("This app lets you explore the global trait space of birds using the BIRDBASE dataset, 
                   via a Gower distance matrix and Principal Coordinates Analysis (PCoA). 
                   You can filter species, choose different trait sets, and then visualize, 
                   compare, and export trait-space information to support exploratory analyses 
                   and idea generation for projects and papers."),
                 
                 h3("How to use the app"),
                 tags$ol(
                   tags$li(
                     tags$b("Start with the Ordination tab:"),
                     " use the filters in the left panel to limit the species pool 
                     (e.g., by order, family, or IUCN category). 
                     Choose a trait preset or custom trait set. 
                     The ordination will update to show species in a 2D trait space, 
                     coloured by the variable you chose (e.g., primary habitat or IUCN status)."
                   ),
                   tags$li(
                     tags$b("Highlight or search specific species:"),
                     " use the 'Highlight species by name' box to bump up the size/alpha 
                     of matching species in the ordination (e.g., type 'loon' or 'gull')."
                   ),
                   tags$li(
                     tags$b("Interpret the axes:"),
                     " scroll down in the Ordination tab to the 'Trait–axis associations' table. 
                     This shows correlations between numeric traits and PCoA1/PCoA2, 
                     helping you interpret what each axis represents 
                     (e.g., body size vs elevational niche vs specialization)."
                   ),
                   tags$li(
                     tags$b("Explore trait distributions:"),
                     " on the Trait distributions tab, inspect histograms of numeric traits 
                     used in the ordination to see how the focal species pool is distributed 
                     in body size, elevation, specialization, etc."
                   ),
                   tags$li(
                     tags$b("Look at the actual rows:"),
                     " the Data tab shows the cleaned trait table used for the current 
                     filters and trait selection, which you can sort, search, and export."
                   ),
                   tags$li(
                     tags$b("Find nearest neighbours:"),
                     " on the Nearest neighbours tab, choose a focal species and get the 
                     closest species in trait space (based on PCoA coordinates). 
                     This is useful for identifying functional analogues or potential 
                     'backup' species if one species is lost."
                   ),
                   tags$li(
                     tags$b("Summarize functional diversity:"),
                     " on the Functional diversity tab, pick a grouping variable 
                     (e.g., order, primary habitat, IUCN category) and see simple 
                     metrics like number of species, centroid in trait space, 
                     mean pairwise distance, mean nearest-neighbour distance, 
                     and 2D hull area. This gives a sense of how much trait 
                     space each group occupies."
                   ),
                   tags$li(
                     tags$b("Compare two subsets:"),
                     " on the Compare subsets tab, choose a variable and two levels 
                     (e.g., LC vs non-LC, migratory vs sedentary, open vs forest habitats). 
                     The ordination will highlight these two subsets plus 'other' species, 
                     and a summary will report the distance between their centroids in trait space."
                   ),
                   tags$li(
                     tags$b("Export for further analysis:"),
                     " in the Ordination tab, use the download buttons to save 
                     the PCoA scores and the filtered trait data, which you can 
                     feed into more formal models or analyses (e.g., GLMs, 
                     phylogenetic models, SDMs, etc.)."
                   )
                 ),
                 
                 h3("Examples of questions you can explore"),
                 tags$ul(
                   tags$li(
                     tags$b("Global context of a focal clade:"),
                     " filter to an order or family (e.g., loons, waterfowl, passerines) 
                     and ask where they sit in global trait space. 
                     Are they concentrated in a narrow region? 
                     Do they span a wide range of body sizes, habitats, or diets?"
                   ),
                   tags$li(
                     tags$b("IUCN status and trait space:"),
                     " colour by IUCN category and turn on hulls. 
                     Are threatened species (VU/EN/CR) clustered in certain parts 
                     of trait space (e.g., large-bodied, specialists), 
                     or scattered throughout?"
                   ),
                   tags$li(
                     tags$b("Functional redundancy:"),
                     " use the Nearest neighbours tab to see whether a threatened 
                     or focal species has many close trait neighbours or is 
                     functionally isolated. This can suggest whether communities 
                     have 'backup' species for certain roles."
                   ),
                   tags$li(
                     tags$b("Movement ecology and traits:"),
                     " filter or colour by migration / movement variables 
                     (e.g., mig, alt_move, sedentary) and check whether 
                     migratory vs sedentary species differ systematically 
                     along body size or specialization axes."
                   ),
                   tags$li(
                     tags$b("Breeding ecology patterns:"),
                     " use a breeding trait preset (clutch size, nest type, 
                     cooperative breeding) and colour by habitat or IUCN status 
                     to see if certain reproductive strategies cluster 
                     in particular environments or risk categories."
                   ),
                   tags$li(
                     tags$b("Linking to regional monitoring data:"),
                     " once you have a list of species occurring in a region 
                     (e.g., boreal wetlands, Prairie Region point-count data), 
                     you can filter to those species here and see how that community 
                     sits in global trait space: Do boreal birds occupy 
                     a narrow slice of global avian traits, or are they 
                     functionally diverse compared to the global pool?"
                   ),
                   tags$li(
                     tags$b("Trait predictors of model performance:"),
                     " if you have SDM or ISDM performance metrics per species, 
                     you can export the trait table and PCoA scores and then test 
                     whether species with poor model performance share 
                     particular traits (e.g., extreme specialists, rare, 
                     unusual movement strategies)."
                   )
                 ),
                 
                 h3("Citation"),
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

## ---- Helper: hull area ---------------------------------------------

hull_area <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  h <- chull(x, y)
  xh <- x[h]
  yh <- y[h]
  area <- 0.5 * abs(sum(xh * c(yh[-1], yh[1]) - yh * c(xh[-1], xh[1])))
  return(area)
}

## ---- Shiny server --------------------------------------------------

server <- function(input, output, session) {
  
  # Trait preset behaviour
  observeEvent(input$trait_preset, {
    if (input$trait_preset == "custom") return(NULL)
    
    preset_traits <- switch(
      input$trait_preset,
      "size_elev" = c("avg_mass", "elev_range", "elev_norm_min", "elev_norm_max"),
      "hab_diet"  = c("primary_habitat", "habitat_breadth",
                      "primary_diet", "diet_breadth", "esi"),
      "breeding"  = c("Clutch_Min", "Clutch_Max", "Nest_Type", "Nest_SBS",
                      "Mono", "Poly", "Coop"),
      "movement"  = c("mig", "alt_move", "irreg_move",
                      "dispersive", "sedentary"),
      "social"    = c("Social_1", "Social_2", "Social_3",
                      "Social_4", "Social_5", "Social_6"),
      character(0)
    )
    
    preset_traits <- preset_traits[preset_traits %in% colnames(bird_traits)]
    
    updatePickerInput(
      session,
      "trait_select",
      selected = preset_traits
    )
  })
  
  # 1) Filter species pool
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
    
    if (nrow(dat) > input$max_species) {
      set.seed(123)
      dat <- dat %>% slice_sample(n = input$max_species)
    }
    
    dat
  })
  
  # 2) Prepared trait data
  trait_data <- reactive({
    req(input$trait_select)
    dat <- filtered_species()
    
    cols_needed <- unique(c(
      "species_label",
      "order", "family_clem",
      "primary_habitat", "primary_diet",
      "iucn2024", "mig",
      "flightlessness",
      "Nest_Type", "Nest_SBS",
      "Mono", "Poly", "Coop",
      paste0("Social_", 1:6),
      input$trait_select
    ))
    
    cols_existing <- cols_needed[cols_needed %in% colnames(dat)]
    dat <- dat %>% select(all_of(cols_existing))
    
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
    
    trait_mat <- dat[, input$trait_select, drop = FALSE]
    
    trait_mat[] <- lapply(trait_mat, function(col) {
      if (is.character(col)) factor(col) else col
    })
    
    gower_dist <- daisy(trait_mat, metric = "gower")
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
  
  # 5) Ordination plot (interactive) with richer colour palette
  output$pcoa_plot <- plotly::renderPlotly({
    res <- pcoa_res()
    validate(
      need(!is.null(res),
           "Not enough species with complete trait data to run PCoA. Increase species pool or relax filters.")
    )
    
    scores <- res$scores
    colour_var <- input$colour_by
    scores$colour_value <- as.factor(scores[[colour_var]])
    
    # dynamic hue palette for arbitrary number of levels
    colour_levels <- levels(scores$colour_value)
    pal_cols <- hue_pal()(length(colour_levels))
    names(pal_cols) <- colour_levels
    
    if (!is.null(input$search_species) && nzchar(input$search_species)) {
      scores$highlight <- grepl(input$search_species,
                                scores$species_label,
                                ignore.case = TRUE)
    } else {
      scores$highlight <- FALSE
    }
    
    base_aes <- aes(
      x = PCoA1, y = PCoA2,
      colour = colour_value,
      text = paste0(
        "<b>", species_label, "</b><br>",
        "Order: ", order, "<br>",
        "Family: ", family_clem, "<br>",
        "Habitat: ", primary_habitat, "<br>",
        "Diet: ", primary_diet, "<br>",
        "IUCN: ", iucn2024, "<br>",
        "Migration: ", mig
      )
    )
    
    p <- ggplot(scores, base_aes) +
      geom_point(
        aes(
          size = ifelse(highlight, 3, 1.5),
          alpha = ifelse(highlight, 1, 0.5)
        )
      ) +
      theme_minimal(base_size = 13) +
      scale_colour_manual(values = pal_cols, na.translate = FALSE) +
      scale_size_identity(guide = "none") +
      scale_alpha_identity(guide = "none") +
      labs(
        colour = colour_var,
        x = "PCoA1",
        y = "PCoA2"
      )
    
    if (isTRUE(input$show_hulls)) {
      hulls <- scores %>%
        filter(!is.na(colour_value)) %>%
        group_by(colour_value) %>%
        filter(n() >= 3) %>%
        slice(chull(PCoA1, PCoA2))
      
      if (nrow(hulls) > 0) {
        p <- p +
          geom_polygon(
            data = hulls,
            aes(x = PCoA1, y = PCoA2,
                fill = colour_value,
                group = colour_value),
            alpha = 0.15,
            colour = NA,
            inherit.aes = FALSE
          ) +
          scale_fill_manual(values = pal_cols, guide = "none")
      }
    }
    
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(legend = list(orientation = "h", y = -0.2))
  })
  
  # 6) Trait histograms (interactive)
  output$trait_hist <- plotly::renderPlotly({
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
    
    p_hist <- ggplot(traits_long, aes(
      x = value,
      fill = trait,
      text = paste0(
        "<b>", species_label, "</b><br>",
        "Trait: ", trait, "<br>",
        "Value: ", signif(value, 3)
      )
    )) +
      geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
      facet_wrap(~ trait, scales = "free") +
      theme_minimal(base_size = 13) +
      scale_fill_brewer(palette = "Set3") +
      labs(x = "Trait value", y = "Count", fill = "Trait")
    
    plotly::ggplotly(p_hist, tooltip = "text")
  })
  
  # 7) Data table
  output$trait_table <- renderDT({
    trait_data()
  })
  
  # 8) Nearest neighbours – server-side selectize for species list
  observeEvent(trait_data(), {
    dat <- trait_data()
    species_choices <- sort(unique(dat$species_label))
    updateSelectizeInput(
      session,
      "nn_species",
      choices = species_choices,
      server = TRUE
    )
  })
  
  nn_table <- reactive({
    res <- pcoa_res()
    validate(need(!is.null(res), "PCoA must be available to compute nearest neighbours."))
    scores <- res$scores
    req(input$nn_species, input$nn_k)
    
    focal_row <- scores %>% filter(species_label == input$nn_species)
    if (nrow(focal_row) == 0) return(NULL)
    
    coords <- as.matrix(scores[, c("PCoA1", "PCoA2")])
    focal_coords <- as.numeric(focal_row[1, c("PCoA1", "PCoA2")])
    
    dists <- sqrt((coords[, 1] - focal_coords[1])^2 +
                    (coords[, 2] - focal_coords[2])^2)
    
    scores$distance <- dists
    scores <- scores %>% filter(species_label != input$nn_species)
    
    scores %>%
      arrange(distance) %>%
      head(input$nn_k) %>%
      select(species_label, order, family_clem,
             primary_habitat, primary_diet, iucn2024, mig,
             PCoA1, PCoA2, distance)
  })
  
  output$nn_table <- renderDT({
    nn <- nn_table()
    validate(need(!is.null(nn), "No nearest neighbours available."))
    datatable(nn, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 9) Functional diversity metrics
  fd_metrics <- reactive({
    res <- pcoa_res()
    validate(need(!is.null(res), "PCoA must be available to compute functional diversity."))
    scores <- res$scores
    group_var <- input$fd_group
    scores$group <- as.factor(scores[[group_var]])
    
    scores <- scores %>% filter(!is.na(group))
    
    out <- scores %>%
      group_by(group) %>%
      group_modify(~{
        df <- .x
        n <- nrow(df)
        if (n < 2) {
          tibble(
            n_species = n,
            centroid_PCoA1 = mean(df$PCoA1),
            centroid_PCoA2 = mean(df$PCoA2),
            mean_pairwise_dist = NA_real_,
            mean_nn_dist = NA_real_,
            hull_area_2D = hull_area(df$PCoA1, df$PCoA2)
          )
        } else {
          coords <- as.matrix(df[, c("PCoA1", "PCoA2")])
          d <- as.matrix(dist(coords))
          mean_pairwise <- mean(d[upper.tri(d)], na.rm = TRUE)
          nn <- apply(d + diag(Inf, n), 1, min, na.rm = TRUE)
          tibble(
            n_species = n,
            centroid_PCoA1 = mean(df$PCoA1),
            centroid_PCoA2 = mean(df$PCoA2),
            mean_pairwise_dist = mean_pairwise,
            mean_nn_dist = mean(nn, na.rm = TRUE),
            hull_area_2D = hull_area(df$PCoA1, df$PCoA2)
          )
        }
      }) %>%
      ungroup()
    
    out
  })
  
  output$fd_table <- renderDT({
    datatable(fd_metrics(), options = list(pageLength = 15, scrollX = TRUE))
  })
  
  # 10) Compare subsets
  output$cmp_level_ui <- renderUI({
    res <- pcoa_res()
    if (is.null(res)) {
      return(helpText("Run PCoA (enough species) to enable comparison."))
    }
    scores <- res$scores
    gv <- input$cmp_group_var
    vals <- sort(unique(na.omit(scores[[gv]])))
    if (length(vals) < 2) {
      return(helpText("Not enough distinct levels in this variable to compare."))
    }
    tagList(
      selectInput("cmp_level_A", "Subset A:", choices = vals, selected = vals[1]),
      selectInput("cmp_level_B", "Subset B:", choices = vals, selected = vals[2])
    )
  })
  
  output$compare_plot <- plotly::renderPlotly({
    res <- pcoa_res()
    validate(need(!is.null(res), "PCoA must be available to compare subsets."))
    scores <- res$scores
    gv <- input$cmp_group_var
    scores$group_val <- scores[[gv]]
    req(input$cmp_level_A, input$cmp_level_B)
    
    scores$cmp_group <- "Other"
    scores$cmp_group[scores$group_val == input$cmp_level_A] <- "A"
    scores$cmp_group[scores$group_val == input$cmp_level_B] <- "B"
    scores$cmp_group <- factor(scores$cmp_group, levels = c("A", "B", "Other"))
    
    p <- ggplot(scores, aes(
      x = PCoA1, y = PCoA2,
      colour = cmp_group,
      text = paste0(
        "<b>", species_label, "</b><br>",
        gv, ": ", group_val
      )
    )) +
      geom_point(alpha = 0.8) +
      theme_minimal(base_size = 13) +
      scale_colour_manual(
        values = c("A" = "#1b9e77", "B" = "#d95f02", "Other" = "grey80")
      ) +
      labs(colour = "Subset")
    
    plotly::ggplotly(p, tooltip = "text")
  })
  
  output$compare_summary <- renderPrint({
    res <- pcoa_res()
    if (is.null(res)) return("PCoA not available.")
    scores <- res$scores
    gv <- input$cmp_group_var
    scores$group_val <- scores[[gv]]
    req(input$cmp_level_A, input$cmp_level_B)
    
    A <- scores %>% filter(group_val == input$cmp_level_A)
    B <- scores %>% filter(group_val == input$cmp_level_B)
    
    if (nrow(A) == 0 || nrow(B) == 0) {
      return("One of the subsets has zero species under current filters.")
    }
    
    centroid_A <- colMeans(A[, c("PCoA1", "PCoA2")])
    centroid_B <- colMeans(B[, c("PCoA1", "PCoA2")])
    dist_centroids <- sqrt(sum((centroid_A - centroid_B)^2))
    
    cat("Subset A:", input$cmp_level_A, " (n =", nrow(A), ")\n")
    cat("Subset B:", input$cmp_level_B, " (n =", nrow(B), ")\n\n")
    cat("Distance between centroids in PCoA space:", round(dist_centroids, 3), "\n")
  })
  
  # 11) Trait–axis correlations
  axis_traits <- reactive({
    res <- pcoa_res()
    dat <- trait_data()
    if (is.null(res) || nrow(dat) == 0) return(NULL)
    
    scores <- res$scores
    numerics <- input$trait_select[
      sapply(dat[, input$trait_select, drop = FALSE], is.numeric)
    ]
    if (length(numerics) == 0) return(NULL)
    
    out <- lapply(numerics, function(tr) {
      x <- dat[[tr]]
      c1 <- suppressWarnings(cor(x, scores$PCoA1, use = "pairwise.complete.obs"))
      c2 <- suppressWarnings(cor(x, scores$PCoA2, use = "pairwise.complete.obs"))
      data.frame(
        trait = tr,
        cor_PCoA1 = c1,
        cor_PCoA2 = c2,
        abs_cor_PCoA1 = abs(c1),
        abs_cor_PCoA2 = abs(c2)
      )
    })
    
    do.call(rbind, out) %>%
      arrange(desc(abs_cor_PCoA1))
  })
  
  output$axis_traits_table <- renderDT({
    at <- axis_traits()
    validate(need(!is.null(at), "No numeric traits to correlate with PCoA axes."))
    datatable(at, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 12) Downloads
  output$download_scores <- downloadHandler(
    filename = function() {
      paste0("pcoa_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- pcoa_res()
      validate(need(!is.null(res), "PCoA not available to download scores."))
      write.csv(res$scores, file, row.names = FALSE)
    }
  )
  
  output$download_traits <- downloadHandler(
    filename = function() {
      paste0("traits_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      dat <- trait_data()
      write.csv(dat, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
