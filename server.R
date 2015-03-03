library("ggplot2")
library("reshape2")
library("tidyr")
library("dplyr")
library("gridExtra")
# following meeting with Enrico Glaab
library("vsn")
library("limma")
library("impute")
library("shiny")
source("utils.R") # own function

# 140201 project for Susann Mudrich

shinyServer(function(input, output) {

  # run once when the app is started
  mass <- read.table("mass.tsv", sep = "\t", header = TRUE) %>%
          tbl_df()  
  # reactive, main computation done for each changes
  # and common to all outputs
  massk <- reactive({
    # Filter for proteins that are found 
    # in all conditions, and at least 4 replicates
    mass.cast <- mass %>%
      group_by(Accession, Time, Condition) %>%
      mutate(n = n()) %>% 
      # at least in 5 replicates
      filter(n >= input$MIN_REP) %>% 
      group_by(Accession, Time) %>%
      mutate(count = n(), Score.med = median(Score)) %>%
      # at least in all experiments 
      filter(count >= (input$MIN_REP * 4)) %>%
      ## try to keep for imputation
      # first return to wide format
      dcast(Accession ~ Condition + Replicate, value.var = "Score")
    # convert to matrix after removing acc numbers and bad replicate
    mass.mat <- mass.cast %>%
      select(-c(Accession, IL27_7)) %>% as.matrix
    # try an imputation from RNA seq package
    if(input$IMPUTE) {
      mass.mat.knn <- impute.knn(mass.mat)$data
    }else {
      mass.mat.knn <- mass.mat
    }
    # do variance stabilization 
    vsndat <- vsn2(mass.mat.knn)
    vsndat <- exprs(vsndat)
    # return a list to get all objects
    return(list(cast = mass.cast,
                raw = mass.mat.knn,
                vsn = vsndat))
  })
  
  # plot variance stabization before and after
  output$massPlot <- renderPlot({
    # from Enrico Glaab, 
    # normalization should be first a variance stabilizing transformation
    # see issue:
    par(mfrow = c(1, 2))
    meanSdPlot(massk()$raw, main = "before vsn2")
    # correct using vsn2
    vsndat = vsn2(massk()$raw)
    vsndat = exprs(vsndat)
    # plot again
    meanSdPlot(massk()$vsn, main = "after vsn2")
    # repeat for saving plots
    pdf("vsn.pdf", height = 5, width = 9)
    par(mfrow = c(1, 2))
    meanSdPlot(massk()$raw, main = "before vsn2")
    meanSdPlot(massk()$vsn, main = "after vsn2")
    dev.off()
  })
  output$downloadVSN <- downloadHandler(
    filename = function() {
      "vsn.pdf"
    },
    content = function(file) {
      file.copy("vsn.pdf", file, overwrite = TRUE)
    }
  )
    
  # MDS for quality control
  output$MDS <- renderPlot({
    p <- plot.mds(massk()$vsn)
    plot(p)
    ggsave("MDS.pdf", height = 7, width = 9, p)
  })
  output$downloadMDS <- downloadHandler(
    filename = function() {
      "MDS.pdf"
    },
    content = function(file) {
      file.copy("MDS.pdf", file, overwrite = TRUE)
    }
  )
  # PCA for quality control
  output$PCA <- renderPlot({
    p <- plot.PCA(massk()$vsn)
    plot(p)
    ggsave("PCA.pdf", height = 7, width = 9, p)
  })
  output$downloadPCA <- downloadHandler(
    filename = function() {
      "PCA.pdf"
    },
    content = function(file) {
      file.copy("PCA.pdf", file, overwrite = TRUE)
    }
  )
  
  
  # Hierarchical clustering for quality controls
  output$hclust <- renderPlot({
    hcl = hclust(dist(t(massk()$vsn)))
    plot(hcl)
    pdf("hclust.pdf")
    plot(hcl)
    dev.off()
  })
  output$downloadHCLUST <- downloadHandler(
    filename = function() {
      "hclust.pdf"
    },
    content = function(file) {
      file.copy("hclust.pdf", file, overwrite = TRUE)
    }
  )
  
  # inform user on 
  output$retained <- renderText({
    paste("Retained proteins:", nrow(massk()$vsn), "in at least", input$MIN_REP,
          "replicates for the 4 conditions (", length(unique(mass$Accession)), "initially )")
  })
  
  output$topTable <- renderDataTable({
    if(input$cmp == "IL27") {
      sampleID <- c(1:9,1:8)
    } else{
      sampleID <- c(1:9,1:9)
    }
    # constrat using limma
    groups <- lapply(strsplit(colnames(massk()$vsn), "_"), function(x) x[1]) %>% as.character
    fit <- contrasts_limma(massk()$vsn, groups, sampleID, "control", input$cmp)
    # format as table and retrieve protein ID
    res <- topAccession(fit, massk()$cast, nrow(fit))
    output$downloadTable <- downloadHandler(
      filename = function() { paste("topTable", '.tsv', sep='') },
      content = function(file) {
        write.table(select(res, -Access), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    # create the html link
    res$Access <- paste("<a href='http://www.uniprot.org/uniprot/", 
                           res$Accession, "', target='_blank'>", res$Accession, "</a>", sep = "")
    res2 <- select(res, Access, logFC:B)
    res2
    # tweak number of entries to show by default
  },options = list(lengthMenu = c(10, 20, 50, 100), 
                   pageLength = 20,
                   # column sorted, highlighthed
                   orderClasses = TRUE),
  # escape html tags except for column 1: accession number
  escape = -1)
})




