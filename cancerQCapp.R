# Shiny web application for viewing cancer QC data
# Martina Mijuskovic
# May 2017

library(shiny)
library(dplyr)
library(data.table)
library(ggplot2)

#### Load and prepare objects

# Tumor sample groups
tumor <- c("FF", "FFPE")
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)
# Bigger plot text
bigger <- theme(legend.text=element_text(size=15), legend.title = element_text(size=15), axis.title = element_text(size=15), axis.text = element_text(size=15))
# Tilted x-axis labels
tiltedX <- theme(axis.text.x=element_text(angle=45,hjust=1))

# Load and clean QC metrics data (NEEDS MANUAL UPDATE)
QC <- read.csv("./Data/ready_to_plot.05-02-17.csv")
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

# Load GMC codes
GMCs <- read.csv("./Data/GMCs.csv")
GMCs <- GMCs %>% select(CODE, GMC)
GMCs$CODE <- as.character(GMCs$CODE)


# Define UI (user interface) - HTML code
ui <- fluidPage(
  
  # Application title
  titlePanel("Cancer Sample QC"),
  
  navlistPanel(
    tabPanel("SETTINGS",
                 mainPanel(
                   checkboxGroupInput(inputId = "gmc",
                                      label = "Select GMC(s):",
                                      choices = paste0(GMCs$GMC, " - ", GMCs$CODE),
                                      selected = "GMC1 - RGT"
                                      ),
                   
                   dateRangeInput(inputId = "dateRange",
                                  label = "Collection date:",
                                  start = "2015-06-23",
                                  end = "2017-01-06",
                                  min = "2015-06-23",
                                  max = "2017-01-06",
                                  startview = "year"
                                  ),

                   fileInput(inputId = "QC_file",
                            label = "Upload QC table:")
                  )
              ),
                        
    tabPanel("TUMOUR TYPE SUMMARY", 
             mainPanel(
               plotOutput("TumourPlot")
             )),
    
    tabPanel("SAMPLE TYPE SUMMARY", 
             sidebarPanel(
               tableOutput("DistributionTable")
             ),
             mainPanel(
               plotOutput("TumourSampleTypePlot")
             )),
    tabPanel("AT DROPOUT",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
             )),
             mainPanel(
               plotOutput("ATdropPlot")
             )
      
    ),
    
    tabPanel("COVERAGE UNEVENNESS",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type2",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type2",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("UnCoveragePlot")
             )
             
    ),
    
    tabPanel("MAPPING RATE",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type3",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type3",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("MappingPlot")
             )
             
    ),
    
    tabPanel("CHIMERIC READS",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type4",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type4",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("ChimericPlot")
             )
             
    ),
    
    tabPanel("DEAMINATION",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type5",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type5",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("DeaminationPlot")
             )
             
    ),
    
    tabPanel("GENOME COVERAGE",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type6",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type6",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("CoveragePlot")
             )
             
    ),
    
    tabPanel("COSMIC COVERAGE",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type8",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type8",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("CosCoveragePlot")
             )
             
    ),
    
    tabPanel("DUPLICATION RATE",
             sidebarPanel(
               checkboxInput(inputId = "by_tumor_type7",
                             label = "By tumour type",
                             value = FALSE
               ),
               checkboxInput(inputId = "by_sample_type7",
                             label = "By sample type (FF/FFPE)",
                             value = TRUE
               )),
             mainPanel(
               plotOutput("DuplicationPlot")
             )
             
    )
    
    
  )        
)



# Define server logic
server <- function(input, output) {
   
   output$TumourPlot <- renderPlot({
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
     FreqData <- as.data.frame(table(QC_tumor[QC_tumor$CENTER %in% center,]$TUMOUR_TYPE))
     names(FreqData) <- c("TUMOUR TYPE", "FREQ")
     ggplot(FreqData, aes(x="", y=FREQ, fill = `TUMOUR TYPE`)) + 
       geom_bar(width = 1, stat = "identity") + 
       theme_void() + 
       coord_polar("y", start=0) +
       theme(legend.text=element_text(size=15), legend.title = element_text(size=15))
   })
   
   output$DistributionTable <- renderTable({
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
     if (input$by_tumor_type & input$by_sample_type){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
     } 
     else if (input$by_sample_type) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=AT_DROP, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
     }
     else if (input$by_tumor_type) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=AT_DROP, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=AT_DROP, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "A/T Dropout") + bigger
     }
   })
   
   output$UnCoveragePlot <- renderPlot({
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
     if (input$by_tumor_type2 & input$by_sample_type2){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage")  + bigger + tiltedX
     } 
     else if (input$by_sample_type2) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COVERAGE_HOMOGENEITY, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
     }
     else if (input$by_tumor_type2) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COVERAGE_HOMOGENEITY, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COVERAGE_HOMOGENEITY, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "Unevenness of coverage") + bigger
     }
   })
   
   output$MappingPlot <- renderPlot({
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
     if (input$by_tumor_type3 & input$by_sample_type3){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate")  + bigger + tiltedX
     } 
     else if (input$by_sample_type3) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MAPPING_RATE_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
     }
     else if (input$by_tumor_type3) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MAPPING_RATE_PER, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MAPPING_RATE_PER, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "Mapping rate") + bigger
     }
   })
   
   output$ChimericPlot <- renderPlot({
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
     if (input$by_tumor_type4 & input$by_sample_type4){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads")  + bigger + tiltedX
     } 
     else if (input$by_sample_type4) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=CHIMERIC_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
     }
     else if (input$by_tumor_type4) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=CHIMERIC_PER, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=CHIMERIC_PER, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "% chimeric reads") + bigger
     }
   })
   
   output$DeaminationPlot <- renderPlot({
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
     if (input$by_tumor_type5 & input$by_sample_type5){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination")  + bigger + tiltedX
     } 
     else if (input$by_sample_type5) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DEAMINATION_MISMATCHES_PER, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
     }
     else if (input$by_tumor_type5) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DEAMINATION_MISMATCHES_PER, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DEAMINATION_MISMATCHES_PER, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "% deamination") + bigger
     }
   })
   
   output$CoveragePlot <- renderPlot({
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
     if (input$by_tumor_type6 & input$by_sample_type6){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage")  + bigger + tiltedX
     } 
     else if (input$by_sample_type6) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MEDIAN_COV, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
     }
     else if (input$by_tumor_type6) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=MEDIAN_COV, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=MEDIAN_COV, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "Median coverage") + bigger
     }
   })
   
   output$DuplicationPlot <- renderPlot({
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
     if (input$by_tumor_type7 & input$by_sample_type7){
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate")  + bigger + tiltedX
     } 
     else if (input$by_sample_type7) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DUPL_RATE, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
     }
     else if (input$by_tumor_type7) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=DUPL_RATE, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=DUPL_RATE, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "Duplication rate") + bigger
     }
   })
   
   output$CosCoveragePlot <- renderPlot({
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
     if (input$by_tumor_type8 & input$by_sample_type8) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")   + bigger + tiltedX
     } 
     else if (input$by_sample_type8) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COSMIC_COV_LT30X, colour = GROUP)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
     }
     else if (input$by_tumor_type8) {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=TUMOUR_TYPE, y=COSMIC_COV_LT30X, colour = TUMOUR_TYPE)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger + tiltedX
     }
     else {
       ggplot(QC_tumor[QC_tumor$CENTER %in% center,], aes(x=CENTER, y=COSMIC_COV_LT30X, colour = CENTER)) + geom_boxplot() + labs(x = "", y = "% Cosmic regions < 30X")  + bigger
     }
   })
   

     
}

# Run the application 
shinyApp(ui = ui, server = server)

