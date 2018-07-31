#### cancerQC v2 single-file app

library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load test data
load("test_QC_data.RData")

### Helper functions for plotting

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Bigger plot text
bigger <- theme(legend.text=element_text(size=15), legend.title = element_text(size=15), axis.title = element_text(size=15), axis.text = element_text(size=15))
# Tilted x-axis labels
tiltedX <- theme(axis.text.x=element_text(angle=45,hjust=1))

# Show GMCs present in the test data
#gmc_choices <- paste0(GMCs$GMC, " - ", GMCs$CODE)
gmc_choices <- levels(QC_tumor$CENTER_annon)
                      
ui <- fluidPage(
  
  # Application title
  titlePanel("Cancer Sample QC"),
  h4("Visualising the quality of whole genome sequencing data"),
  a(href="https://github.com/martimij/CancerQC_v2", "Documentation (GitHub)"),
  br(),
  br(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      checkboxGroupInput(inputId = "gmc",
                         label = "Select Genomic Medicine Center(s):",
                         choices = gmc_choices),
      
      checkboxInput(inputId = 'all_none', 
                    label = 'All/None',
                    value = TRUE),
      
      dateRangeInput(inputId = "dateRange",
                     label = "Sample collection date:",
                     start = "2015-06-23",
                     end = "2017-01-06",
                     min = "2015-06-23",
                     max = "2017-01-06",
                     startview = "month"),
      
      radioButtons(inputId = "group_by",
                   label = "Group by:",
                   choices = list("TUMOUR TYPE" = "by_tumour_type", "GMC" = "by_gmc"),
                   selected = "by_gmc"),
      
      checkboxInput(inputId = "by_sample_type",
                    label = "Split by sample type (FF/FFPE)",
                    value = TRUE),
      # p(em("FF = 'Fresh Frozen'")),
      # p(em("FFPE = 'Formalin-Fixed Paraffin Embedded'")),
      
      fileInput(inputId = "QC_file",
                label = "Upload QC table:")
      
    ),
    
    mainPanel(
      tabsetPanel(

          tabPanel("SUMMARY",
                    br(),
                    sidebarPanel(
                      tableOutput("DistributionTable")
                      ),
                   mainPanel(
                     plotOutput("TumourSampleTypePlot")
                   )
          ),
          
          tabPanel("AT DROPOUT",
                   br(),
                   plotOutput("ATdropPlot")
          ),
            
          tabPanel("UNEVENNESS",
                   br(),
                   plotOutput("UnCoveragePlot")
          ),

          tabPanel("MAPPING RATE",
                   br(),
                   plotOutput("MappingPlot")
                     
          ),

          tabPanel("CHIMERAS",
                   br(),
                   plotOutput("ChimericPlot")

          ),

          tabPanel("DEAMINATION",
                   br(),
                   plotOutput("DeaminationPlot")
          ),

          tabPanel("GENOME COVERAGE",
                   br(),
                   plotOutput("CoveragePlot")

          ),

          tabPanel("COSMIC COVERAGE",
                   br(),
                   plotOutput("CosCoveragePlot")
                     
          ),

          tabPanel("DUPLICATION",
                   br(),
                   plotOutput("DuplicationPlot")

          )
      )
    )
  )
) 



server <- function(input, output, session) {

  # Reactivity for the All/None GMC checkbox
  observe({
    updateCheckboxGroupInput(
      session = session,
      inputId = 'gmc',
      choices = gmc_choices,
      selected = if (input$all_none) gmc_choices
    )
  })
  
  output$DistributionTable <- renderTable({
    
    req(input$gmc)
    #center  <- GMCs[GMCs$GMC %in% input$gmc,]$CODE
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
 
     if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      FreqData <- QC_tumor %>% 
        filter(CENTER_annon %in% center) %>% 
        group_by(TUMOUR_TYPE, GROUP) %>% 
        summarise(COUNT = n()) %>% 
        spread(GROUP, COUNT)
      # Add totals
      FreqData <- rbind(as.data.frame(FreqData), data.frame(TUMOUR_TYPE = "TOTAL", t(colSums(FreqData[,-1], na.rm = T))))
      names(FreqData)[1] <- "TUMOUR TYPE"
      # Fix cases when FF is missing
      if (!is.null(FreqData$FF)) {
        FreqData$FF <- as.integer(FreqData$FF)
      }
      # Fix cases when FFPE is missing
      if (!is.null(FreqData$FFPE)) { 
        FreqData$FFPE <- as.integer(FreqData$FFPE)
      }
     FreqData
    }
    
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      FreqData <- QC_tumor %>% 
        filter(CENTER_annon %in% center) %>% 
        group_by(GROUP, CENTER_annon) %>% 
        summarise(COUNT = n()) %>% 
        spread(GROUP, COUNT)
      # Add totals
      FreqData <- rbind(as.data.frame(FreqData), data.frame(CENTER_annon = "TOTAL", t(colSums(FreqData[,-1], na.rm = T))))
      names(FreqData)[1] <- "GMC"
      # Fix cases when FF is missing
      if (!is.null(FreqData$FF)) {
        FreqData$FF <- as.integer(FreqData$FF)
      }
      # Fix cases when FFPE is missing
      if (!is.null(FreqData$FFPE)) { 
        FreqData$FFPE <- as.integer(FreqData$FFPE)
      }
      FreqData

    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      FreqData <- QC_tumor %>% 
        filter(CENTER_annon %in% center) %>% 
        group_by(TUMOUR_TYPE) %>% 
        summarise(COUNT = n())
      # Add totals
      FreqData <- rbind(as.data.frame(FreqData), data.frame(TUMOUR_TYPE = "TOTAL", t(colSums(FreqData[,-1], na.rm = T))))
      names(FreqData)[1] <- "TUMOUR TYPE"
      FreqData$COUNT <- as.integer(FreqData$COUNT)
      FreqData
    }
    else {
      FreqData <- QC_tumor %>% 
        filter(CENTER_annon %in% center) %>% 
        group_by(CENTER_annon) %>% 
        summarise(COUNT = n())
      # Add totals
      FreqData <- rbind(as.data.frame(FreqData), data.frame(CENTER_annon = "TOTAL", t(colSums(FreqData[,-1], na.rm = T))))
      names(FreqData)[1] <- "GMC"
      FreqData$COUNT <- as.integer(FreqData$COUNT)
      FreqData
    }
    
  })
  
  
  output$TumourSampleTypePlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
  
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER_annon %in% center,]$TUMOUR_TYPE, QC_tumor[QC_tumor$CENTER_annon %in% center,]$GROUP))
      names(FreqData) <- c("TUMOUR TYPE", "SAMPLE TYPE", "FREQ")
      ggplot(FreqData, aes(x=`TUMOUR TYPE`, y = FREQ, fill= `SAMPLE TYPE`)) +
        geom_bar(stat = "identity") +
        bigger +
        tiltedX
    }
    
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER_annon %in% center,]$CENTER_annon, QC_tumor[QC_tumor$CENTER_annon %in% center,]$GROUP))
      names(FreqData) <- c("GMC", "SAMPLE TYPE", "FREQ")
      ggplot(FreqData, aes(x=GMC, y = FREQ, fill= `SAMPLE TYPE`)) +
        geom_bar(stat = "identity") +
        bigger +
        tiltedX
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER_annon %in% center,]$TUMOUR_TYPE))
      names(FreqData) <- c("TUMOUR TYPE", "FREQ")
      ggplot(FreqData, aes(x=`TUMOUR TYPE`, y = FREQ)) +
        geom_bar(stat = "identity") +
        bigger +
        tiltedX
    }
    else {
      FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER_annon %in% center,]$CENTER_annon))
      names(FreqData) <- c("GMC", "FREQ")
      ggplot(FreqData, aes(x=GMC, y = FREQ)) +
        geom_bar(stat = "identity") +
        bigger +
        tiltedX
    }
  })
  
  output$ATdropPlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=AT_DROP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
    }

  })
  
  output$UnCoveragePlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=COVERAGE_HOMOGENEITY)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
    }
  })
  
  output$MappingPlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
   
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=MAPPING_RATE_PER)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
    }
  })
  
  output$ChimericPlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
 
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=CHIMERIC_PER)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
    }
  })
  
  output$DeaminationPlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=DEAMINATION_MISMATCHES_PER)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
    }
  })
  
  output$CoveragePlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=MEDIAN_COV)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
    }
  })
  
  output$DuplicationPlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by){
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=DUPL_RATE)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
    }
  })
  
  output$CosCoveragePlot <- renderPlot({
    
    req(input$gmc)
    center <- input$gmc
    start_date <- input$dateRange[1]
    end_date <- input$dateRange[2]
    
    if (!is.null(input$QC_file)){
      QC <- read.csv(as.character(input$QC_file$datapath))
      QC <- QC[!duplicated(QC),]  # remove exact duplicates
      QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
      QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
      QC_tumor <- QC %>% filter(GROUP %in% tumor)
      QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
      QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
      QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
      QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"
    }
    
    QC_tumor <- QC_tumor %>% filter((COLLECTING_DATE >= start_date) & (COLLECTING_DATE <= end_date))
    
    if (input$by_sample_type & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")   + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER_annon %in% center,], aes(x=CENTER_annon, y=COSMIC_COV_LT30X)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

