rm(list = ls())
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(shinythemes)
library(thematic)
library(RColorBrewer)
library(heatmaply)
library(plotly)
data <- read.csv("data/data_long_format.csv")
heat_map_df <- read.csv("data/heatmapdf.csv")
sorted_data <- read.csv("data/data_sorted.csv")
long_heatmap_df<- read.csv("data/long_heatmap_df.csv")
data <- subset(data, select = - X)
rownames(long_heatmap_df) <- long_heatmap_df$X
long_heatmap_df <- subset(long_heatmap_df, select= -X)
sorted_data <- subset(sorted_data, select = -X)
rownames(heat_map_df) <- heat_map_df$X
heat_map_df <- subset(heat_map_df, select = - X)
ui <-page_navbar(
  theme = bs_theme(bootswatch = "cosmo",bg = "white", fg = "black", 
                   heading_font =font_google("Roboto Mono"),
                   base_font = font_google("Roboto Mono"),
                   enable_rounded =TRUE
                   ),
  tags$style(' .well  { background-color: white !important;}'),
  title = "Borrelia Transcriptomic",
  bg = "#0062cc",
  underline = TRUE,
  nav_panel(
    title = "Visualization Tools",
            tabsetPanel(
              tabPanel("Plot",sidebarLayout(
                sidebarPanel(
                             selectInput("var", label = "Select a gene of interest", 
                                         choices = unique(data$Gene), selected = 'BB_0001_-_hypothetical_protein'),
                             checkboxInput("mean_box", label = "Mean value", value = FALSE),
                             checkboxInput("error_bar_box", label = "Error Bar", value = FALSE),
                             checkboxInput("anova_box", label = "Run analysis of varience", value = FALSE)
                ),
                mainPanel( card(
                  height = 600,
                     plotOutput("scatter_plot", width = "100%", height = "1000px"),
                     DTOutput("anova_table")
                ))
              )),
              tabPanel("Heat map", 
                       sidebarLayout(
                         sidebarPanel(
                           selectInput("P_value_filter", label = "P-value less than:", choices = list('5e-2' = 5e-2, '1e-2' = 1e-2, '5e-3' = 5e-3), selected = '5e-2'),
                           
                         ),
                         
                         
                         mainPanel(
                           card( height = 800,
                           plotlyOutput("heat_map",width = "100%",height = "2000px")
                           )
                           )
                       )
                       
                       
                       
                       
                       
                       
                       
              )
            )),
  nav_panel(title = "Datasets", 
            
            tabsetPanel(
              tabPanel(
                "Wide Format Data",
                card( height = 800,
                  DTOutput("wide_data"),
                  downloadButton("downloadWideData", "Download")
                )
              ),
              tabPanel(
                "Long Format Data",card(
                  height = 800,
                  DTOutput("long_data"),
                  downloadButton("downloadLongData", "Download")
                )
              )
            )
            
            ),
  
  
  
  
  
  nav_panel(title = "About",card(card_body(
            h2("Weigang Qiu LAB"),
                  p("Gene expression data was obtained from the original paper", a(em("RNA-Seq of Borrelia burgdorferi in Multiple Phases of Growth Reveals Insights into the Dynamics of Gene Expression, Transcriptome Architecture, and Noncoding RNAs.") ,
                                                                                href= "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164165" )),
            p("The project was executed under the supervision of ", strong("Dr. Weigang, Professor in the Department of Biological Sciences at Hunter College, City University of New York.")),
            p("Bioinformatics Assistant: ", strong("Laziz Asamov, an undergraduate student of Biological Science at Hunter College, City University of New York.")),
            p("Vist", a(em("Weigang Qiu LAB: Evolutionary Bioinformatics"),href = 'https://wiki.genometracker.org/w/Main_Page'), "website for more projects and information." )
))))
server <- function(input, output) {
  
  
  output$scatter_plot <- renderPlot({
    df_mean_std <- data[data$Gene == input$var,] %>%
      group_by(stage) %>%
      summarise_at(vars(logCounts), list(mean=mean, sd=sd)) %>% 
      as.data.frame()
    
    
    if(input$error_bar_box == TRUE & input$mean_box == TRUE){
      
      ggplot(df_mean_std, aes(x = stage, y = mean))+
        geom_point(size = 5, shape = 24, fill="darkturquoise")+
        geom_point(data = data[data$Gene == input$var,], aes(x = stage, y = logCounts), size = 4, shape = 21, fill = 'deepskyblue3', colour = 'deepskyblue3')+
        geom_errorbar(data = df_mean_std, aes(x = stage, y = mean, ymin = mean-sd, ymax = mean+sd,  alpha = 1),size = 1, width = 0.1, linetype = 3, colour = 'red')+
        theme_bw()+
        theme(legend.position = "none", 
              text=element_text(size=16,face="bold", family='Times New Roman'), 
              axis.title=element_text(size=18, family = "Times New Roman"),
              axis.title.x = element_text(margin=margin(20,0,0,0)), 
              axis.title.y = element_text(margin=margin(0,20,0,0)))+
        
        xlab("Stages")+
        ylab("Log 10 of FPKM")}
    
    else if(input$mean_box == TRUE){
      ggplot(df_mean_std, aes(x = stage, y = mean))+
        geom_point(size = 5, shape = 24, fill="darkturquoise")+
        geom_point(data = data[data$Gene == input$var,], aes(x = stage, y = logCounts), size = 4, shape = 21, fill = 'deepskyblue3', colour = 'deepskyblue3')+
        theme_bw()+
        theme(legend.position = "none", 
              text=element_text(size=16,face="bold", family='Times New Roman'), 
              axis.title=element_text(size=18, family = "Times New Roman"),
              axis.title.x = element_text(margin=margin(20,0,0,0)), 
              axis.title.y = element_text(margin=margin(0,20,0,0)))+
        
        xlab("Stages")+
        ylab("Log 10 of RPKM")} 
    else if(input$error_bar_box == TRUE){
      ggplot(data[data$Gene == input$var,], aes(x = stage, y = logCounts))+
        geom_point(size = 4, shape = 21, fill = 'deepskyblue3', colour = 'deepskyblue3')+
        geom_errorbar(data = df_mean_std, aes(x = stage, y = mean, ymin = mean-sd, ymax = mean+sd,  alpha = 1),size = 1, width = 0.1, linetype = 3, colour = 'red')+
        theme_bw()+
        theme(legend.position = "none", 
              text=element_text(size=16,face="bold", family='Times New Roman'), 
              axis.title=element_text(size=18, family = "Times New Roman"),
              axis.title.x = element_text(margin=margin(20,0,0,0)), 
              axis.title.y = element_text(margin=margin(0,20,0,0)))+
        
        xlab("Stages")+
        ylab("Log 10 of FPKM")
    }
    
    
    else{
      ggplot(data[data$Gene == input$var,], aes(x = stage, y = logCounts))+
        geom_point(size = 4, shape = 21, fill = 'deepskyblue3', colour = 'deepskyblue3')+
        theme_bw()+
        theme(legend.position = "none", 
              text=element_text(size=16,face="bold", family='Times New Roman'), 
              axis.title=element_text(size=18, family = "Times New Roman"),
              axis.title.x = element_text(margin=margin(20,0,0,0)), 
              axis.title.y = element_text(margin=margin(0,20,0,0)))+
        
        xlab("Stages")+
        ylab("Log 10 of PKM")
    }
    
    
  })
  output$heat_map<- renderPlotly({

    
    p_heat_map <- long_heatmap_df%>%
      
      filter(long_heatmap_df$P.values < input$P_value_filter)
    
    print(heatmaply(p_heat_map[1:9], showticklabels = c(TRUE, FALSE), cluster_columns = TRUE, cluster_rows = TRUE, column_text_angle = 90,
                    colors = colorRampPalette(brewer.pal(3, "BuPu"))(256),  colorbar_len = 1, dendrogram = 'row'))}
    
  )
  
  
  output$anova_table <- renderDT({
    if (input$anova_box == TRUE){
      
      model <- lm(counts ~ stage, data = data[data$Gene == input$var,])
      anv <- anova(model)
      datatable(anv, options = list(dom = 't'))
      
      
    }
    
  })
  output$wide_data <- renderDT({
   datatable(long_heatmap_df)
  })
  
  output$long_data <- renderDT(
    datatable(data)
  )
  output$downloadWideData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(long_heatmap_df, file)
    }
  )
  output$downloadLongData <- downloadHandler(
    
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data, file)
    }
  )
  
  
  
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)