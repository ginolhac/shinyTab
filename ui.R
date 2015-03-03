library("shiny")

shinyUI(fluidPage(theme = "bootstrap.css",
  
  # Application title
  title = "Mass spectrometry",
  p("replicate 7 for IL27 is discarded"),
  # add link to the rationale page
  tags$div(class="pull-right",
    tags$a(href = "rationale.html", "Rationale")
  ),
  # first row for input controls
  fluidRow(
    column(3,
           numericInput("MIN_REP", label = h6("protein in at least x / 9 replicates [1-9]"), value = 5
           )
    ),
    column(2,
           checkboxInput("IMPUTE", label = h6("Imputation of missing values"), value = TRUE
          )
    ),
    column(3,
           radioButtons("cmp", label = h5("Control versus"),
                        choices = list("IL27" = "IL27", 
                                       "calprotectin" = "calprotectin",
                                       "IL27 + calprotectin" = "IL27calprotectin"), 
                        selected = "IL27"
           )
    )
  ) ,
  # then results row
  # don't plot thing if replicates is not meaningful
  conditionalPanel("input.MIN_REP < 10 & input.MIN_REP > 0",
                   verbatimTextOutput("retained"),
                   # start navigation by tabs
                   tabsetPanel(
                     tabPanel("Variance stabilization", plotOutput(outputId = "massPlot"),
                              downloadButton('downloadVSN', 'Download')), 
                     tabPanel("MDS", plotOutput(outputId = "MDS"),
                              downloadButton('downloadMDS', 'Download')),
                     tabPanel("PCA", plotOutput(outputId = "PCA"),
                              downloadButton('downloadPCA', 'Download')),
                     tabPanel("H. Clusters", plotOutput(outputId = "hclust"),
                              downloadButton('downloadHCLUST', 'Download')),
                     tabPanel("Top hits", dataTableOutput(outputId = "topTable"),
                              downloadButton('downloadTable', 'Download'))
                   )
  )
))
