#### cancerQC v2 single-file app

library(shiny)
library(dplyr)
library(data.table)
library(ggplot2)

# Load test data
load("test_QC_data.RData")

### Helper functions for plotting

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Bigger plot text
bigger <- theme(legend.text=element_text(size=15), legend.title = element_text(size=15), axis.title = element_text(size=15), axis.text = element_text(size=15))
# Tilted x-axis labels
tiltedX <- theme(axis.text.x=element_text(angle=45,hjust=1))

gmc_choices <- paste0(GMCs$GMC, " - ", GMCs$CODE)

ui <- fluidPage(
  
  # Application title
  titlePanel("Cancer Sample QC"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      checkboxGroupInput(inputId = "gmc",
                         label = "Select Genomic Medicine Center(s):",
                         choices = gmc_choices),
      
      checkboxInput(inputId = 'all_none', 
                    label = 'All/None',
                    value = TRUE),
      
      dateRangeInput(inputId = "dateRange",
                     label = "Collection date:",
                     start = "2015-06-23",
                     end = "2017-01-06",
                     min = "2015-06-23",
                     max = "2017-01-06",
                     startview = "year"),
      
      radioButtons(inputId = "group_by",
                   label = "Group by:",
                   choices = list("Tumour type" = "by_tumour_type", "GMC" = "by_gmc"),
                   selected = "by_gmc"),
      
      checkboxInput(inputId = "by_sample_type",
                    label = "Split by sample type (FF/FFPE)",
                    value = TRUE),
      
      fileInput(inputId = "QC_file",
                label = "Upload QC table:")
      
    ),
    
    mainPanel(
      tabsetPanel(

          tabPanel("SUMMARY",
                    sidebarPanel(
                      tableOutput("DistributionTable")
                      ),
                   mainPanel(
                     plotOutput("TumourSampleTypePlot")
                   )
          ),
          
          tabPanel("AT DROPOUT",
                   plotOutput("ATdropPlot")
          ),
            
          tabPanel("COVERAGE UNEVENNESS",
                       plotOutput("UnCoveragePlot")
          ),

          tabPanel("MAPPING RATE",
                       plotOutput("MappingPlot")
                     
          ),

          tabPanel("CHIMERIC READS",
                       plotOutput("ChimericPlot")

          ),

          tabPanel("DEAMINATION",
                       plotOutput("DeaminationPlot")
          ),

          tabPanel("GENOME COVERAGE",
                       plotOutput("CoveragePlot")

          ),

          tabPanel("COSMIC COVERAGE",
                       plotOutput("CosCoveragePlot")
                     
          ),

          tabPanel("DUPLICATION",
                       plotOutput("DuplicationPlot")

          )
      )
    )
  )
) 



server <- function(input, output, session) {

  # Reactivity for the All/None GMC checkbox (doesn't work)
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
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
    table <- as.data.frame.matrix(table(QC_tumor[QC_tumor$CENTER %in% center,]$TUMOUR_TYPE, QC_tumor[QC_tumor$CENTER %in% center,]$GROUP))
    table <- tibble::rownames_to_column(table)
    names(table)[1] <- "Tumour Type"
    table
  })
  
  output$TumourSampleTypePlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
    FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER %in% center,]$TUMOUR_TYPE, QC_tumor[QC_tumor$CENTER %in% center,]$GROUP))
    names(FreqData) <- c("TUMOUR TYPE", "SAMPLE TYPE", "FREQ")
    ggplot(FreqData, aes(x=`TUMOUR TYPE`, y = FREQ, fill= `SAMPLE TYPE`)) +
      geom_bar(stat = "identity") +
      bigger +
      tiltedX
  })
  
  output$ATdropPlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=AT_DROP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
    }

  })
  
  output$UnCoveragePlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COVERAGE_HOMOGENEITY)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
    }
  })
  
  output$MappingPlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MAPPING_RATE_PER)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
    }
  })
  
  output$ChimericPlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=CHIMERIC_PER)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
    }
  })
  
  output$DeaminationPlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DEAMINATION_MISMATCHES_PER)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
    }
  })
  
  output$CoveragePlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MEDIAN_COV)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
    }
  })
  
  output$DuplicationPlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate")  + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DUPL_RATE)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
    }
  })
  
  output$CosCoveragePlot <- renderPlot({
    req(input$gmc)
    gmcs <- sapply(1:length(input$gmc), function(x){
      strsplit(input$gmc[x], split = " - ")[[1]][1]
    })
    center  <- GMCs[GMCs$GMC %in% gmcs,]$CODE
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
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")   + bigger + tiltedX
    } 
    else if (input$by_sample_type & !("by_tumour_type" %in% input$group_by)) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
    }
    else if (!(input$by_sample_type) & "by_tumour_type" %in% input$group_by) {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger + tiltedX
    }
    else {
      ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COSMIC_COV_LT30X)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
    }
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

