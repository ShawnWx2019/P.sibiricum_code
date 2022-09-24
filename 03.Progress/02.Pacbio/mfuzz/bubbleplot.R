##########################################################
#     prj:  Enrichment result Enrichment Visualization shiny localhost ver.
#     date: 14 Feb 2021
#     Author: Shawn Wang
##########################################################
## options and import packages
options(stringsAsFactors = F)
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
jscode <- "shinyjs.closeWindow = function() { window.close(); }" # using shinyjs to stop shiny app
if (!require('shiny')) install.packages('shiny')
library(shiny)
if (!require('htmltools')) install.packages('htmltools')
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('dplyr')) install.packages('dplyr')
if (!require('stringr')) install.packages('stringr')
if (!require('magrittr')) install.packages('magrittr')
if (!require('shinyjs')) install.packages('shinyjs')
if (!require('colourpicker')) install.packages('colourpicker')
library(shiny)
library(shinyjs)
library(colourpicker)
library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)
## functions
## 01.GO Bubble Plot
GOBubblePlot = function(expfile, num, name,colormin,colormax,Gtype){
  ## 01 Import raw results
  data.raw = expfile
  ## ordered by qvalue
  data.raw = data.raw[order(data.raw$corrected.p.value.BH.method.),]
  id = as.character(data.raw$GO_Name)
  ## 02 Data cleaning
  data.clean = data.frame(id = id,
                          ratio = data.raw$HitsGenesCountsInSelectedSet/data.raw$AllGenesCountsInSelectedSet,
                          Gene_Number = data.raw$HitsGenesCountsInSelectedSet,
                          Qvalue = data.raw$corrected.p.value.BH.method.)
  ## how many go terms you want to show
  if (nrow(data.clean) <= num){
    data.clean <- data.clean
  } else {
    data.clean <- data.clean[c(1:num),]
  }
  ## rearrangement of dot order by ratio
  data.clean = data.clean[order(data.clean$ratio),]
  data.clean$id = factor(data.clean$id, levels = data.clean$id)
  y <- as.character(data.clean[,1])
  ## set y as factor for ggplot
  y <- factor(y,levels = unique(y))
  x <- data.clean[,2]
  Gene_Number <- data.clean[,3]
  Qvalue <- data.clean[,4]
  title <- paste(name,"Top",num,"GO enrichment result",sep = " ")
  font.size = 12
  
  ## plot step
  if (Gtype == "Bubble") {
    p.base = ggplot(data.clean, aes(x=x, y=y))+
      geom_point(aes(color = Qvalue, size = Gene_Number))+
      scale_color_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      scale_size_continuous(range=c(0, 8))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Ratio")
  } else  {
    p.base = ggplot(data.clean, aes(x=Gene_Number, y=y,fill = Qvalue))+
      geom_bar(stat = "identity")+
      scale_fill_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Number")
  }
  out = list(y = y,
             p = p.base)
  return(p.base)
  
}
## 02.KEGG
KEGGBubblePlot = function(expfile, num, name,colormin,colormax,Gtype){
  ## 01 Import raw results
  data.raw = expfile
  ## ordered by qvalue
  data.raw = data.raw[order(data.raw$corrected.p.value.BH.method.),]
  id = as.character(data.raw$Term.Name)
  colnames(data.raw)
  ## 02 Data cleaning
  data.clean = data.frame(id = id,
                          ratio = data.raw[,3]/data.raw[,4],
                          Gene_Number = data.raw[,3],
                          Qvalue = data.raw[,10])
  ## how many go terms you want to show
  if (nrow(data.clean) <= num){
    data.clean <- data.clean
  } else {
    data.clean <- data.clean[c(1:num),]
  }
  ## rearrangement of dot order by ratio
  data.clean = data.clean[order(data.clean$ratio),]
  data.clean$id = factor(data.clean$id, levels = data.clean$id)
  y <- as.character(data.clean[,1])
  ## set the graph size
  ## set y as factor for ggplot
  y <- factor(y,levels = unique(y))
  x <- data.clean[,2]
  Gene_Number <- data.clean[,3]
  Qvalue <- data.clean[,4]
  title <- paste(name,"Top",num,"GO enrichment result",sep = " ")
  font.size = 12
  ## plot step
  if (Gtype == "Bubble") {
    p.base = ggplot(data.clean, aes(x=x, y=y))+
      geom_point(aes(color = Qvalue, size = Gene_Number))+
      scale_color_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      scale_size_continuous(range=c(0, 8))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Ratio")
  } else  {
    p.base = ggplot(data.clean, aes(x=Gene_Number, y=y,fill = Qvalue))+
      geom_bar(stat = "identity")+
      scale_fill_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Number")
  }
  
  out = list(y = y,
             p = p.base)
  return(p.base)
}
## 03.Customized
customBubblePlot = function(expfile, num, name,colormin,colormax,Gtype){
  ## 01 Import raw results
  data.raw = expfile
  ## ordered by qvalue
  data.raw = data.raw[order(data.raw[,4]),]
  id = as.character(data.raw[,1])
  colnames(data.raw)
  ## 02 Data cleaning
  data.clean = data.frame(id = id,
                          ratio = data.raw[,2],
                          Gene_Number = data.raw[,3],
                          Qvalue = data.raw[,4])
  ## how many go terms you want to show
  if (nrow(data.clean) <= num){
    data.clean <- data.clean
  } else {
    data.clean <- data.clean[c(1:num),]
  }
  ## rearrangement of dot order by ratio
  data.clean = data.clean[order(data.clean$ratio),]
  data.clean$id = factor(data.clean$id, levels = data.clean$id)
  y <- as.character(data.clean[,1])
  ## set the graph size
  ## set y as factor for ggplot
  y <- factor(y,levels = unique(y))
  x <- data.clean[,2]
  Gene_Number <- data.clean[,3]
  Qvalue <- data.clean[,4]
  title <- paste(name,"Top",num,"GO enrichment result",sep = " ")
  font.size = 12
  ## plot step
  if (Gtype == "Bubble") {
    p.base = ggplot(data.clean, aes(x=x, y=y))+
      geom_point(aes(color = Qvalue, size = Gene_Number))+
      scale_color_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      scale_size_continuous(range=c(0, 8))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Ratio")
  } else  {
    p.base = ggplot(data.clean, aes(x=Gene_Number, y=y,fill = Qvalue))+
      geom_bar(stat = "identity")+
      scale_fill_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "Gene_Number")
  }
  out = list(y = y,
             p = p.base)
  return(p.base)
}
## 04.GSEA
GSEAplot = function(expfile, num, name,colormin,colormax,Gtype){
  ## 01 Import raw results
  data.raw = expfile
  ## ordered by qvalue
  data.raw = data.raw[order(data.raw[,4]),]
  id = as.character(data.raw[,1])
  colnames(data.raw)
  ## 02 Data cleaning
  data.clean = data.frame(id = id,
                          ratio = data.raw[,2],
                          Gene_Number = data.raw[,3],
                          Qvalue = data.raw[,4])
  ## how many go terms you want to show
  if (nrow(data.clean) <= num){
    data.clean <- data.clean
  } else {
    data.clean <- data.clean[c(1:num),]
  }
  ## rearrangement of dot order by ratio
  data.clean = data.clean[order(data.clean$ratio),]
  data.clean$id = factor(data.clean$id, levels = data.clean$id)
  y <- as.character(data.clean[,1])
  ## set the graph size
  ## set y as factor for ggplot
  y <- factor(y,levels = unique(y))
  x <- data.clean[,2]
  Gene_Number <- data.clean[,3]
  Pvalue <- data.clean[,4]
  title <- paste(name,"Top",num,"GO enrichment result",sep = " ")
  font.size = 12
  ## plot step
  if (Gtype == "Bubble") {
    p.base = ggplot(data.clean, aes(x=x, y=y))+
      geom_point(aes(color = Pvalue, size = Gene_Number))+
      scale_color_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      scale_size_continuous(range=c(0, 8))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "NES")
  } else  {
    p.base = ggplot(data.clean, aes(x=x, y=y,fill = Pvalue))+
      geom_bar(stat = "identity")+
      scale_fill_continuous(low=colormin, high=colormax, guide=guide_colorbar(reverse=TRUE))+
      labs(y=NULL) +
      ggtitle(title)+
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))+ 
      labs(x = "NES")
  }
  out = list(y = y,
             p = p.base)
  return(p.base)
}
## ShinyUI setting
ui = shinyUI(fluidPage(
  titlePanel(h4("Enrichment Bubble Plot")),
  
  # stop shiny app via js code
  useShinyjs(),
  extendShinyjs(text = jscode, functions = c("closeWindow")),
  actionButton("close", "Stop the App", class = "btn-warning"),
  # sidebar setting
  sidebarLayout(
    sidebarPanel(
      # Input files
      fileInput(inputId = "file",
                label = "Upload enrichment result"),
      help("The results of TBtools enrichment can be uploaded directly, results generated from other softwares need to be input in the required format"),
      # Select the input file type.
      radioButtons(inputId = "DataF",
                   label = "Input File Type",
                   choices = c(TBtools_GO = "TBtools_GO",
                               TBtools_KEGG = "TBtools_KEGG",
                               GSEA = "GSEA",
                               customized = "customized"),
                   selected = "TBtools_GO"),
      # Term number, how many terms do you want to demonstrate
      sliderInput(inputId = "num",
                  label = "Select the value from Slider",
                  min = 10,
                  max = 150,
                  value = 30),
      # Min color
      colourpicker::colourInput(inputId = "colormin",
                  label = "Pvalue color Minimum",
                  value = "red"),
      # Max color
      colourpicker::colourInput(inputId = "colormax",
                  label = "Pvalue color Maxmum",
                  value = "blue"),
      # Graph type, bubble plot or bar plot
      radioButtons(inputId = "Gtype",
                   label = "Graph type",
                   choices = c(Bubble = "Bubble",
                               bar = "bar"),
                   selected = "Bubble"),
      # Project name
      textInput(inputId = "name",
                label = "Project Name",
                value = ""),
      # output graph size width
      sliderInput(inputId = "w",
                  label = "Graph size-width",
                  min = 1,
                  max = 20,
                  value = 10),
      # output graph size height
      sliderInput(inputId = "h",
                  label = "Graph size-height",
                  min = 1,
                  max = 20,
                  value = 10),
      p( "If there is a longer GO term or Pathway name, it is recommended to increase the width of the picture. If the number of terms you select is to large, it is recommended to increase the height of the picture.")
    ),
    # demonstration part
    mainPanel(
      uiOutput("tb")
      )
    )
  )
)
# Server setting
server <- shinyServer(function(input,output){
  ## set para
  DataF <- reactive({
    as.character(input$DataF)
  })
  num <- reactive({
    as.numeric(input$num)
  })
  name <- reactive({
    as.character(input$name)
  })
  colormin <- reactive({
    as.character(input$colormin)
  })
  colormax <- reactive({
    as.character(input$colormax)
  })
  Gtype <- reactive({
    as.character(input$Gtype)
  })
  w <- reactive({
    as.numeric(input$w)
  })
  h <- reactive({
    as.numeric(input$h)
  })
  # import file
  data <- reactive({
    file1 <- input$file
    if(is.null(file1)){return()}
    read.table(file = file1$datapath,
               sep="\t",
               header = T,
               stringsAsFactors = F)
  })
  ## show tables
  output$finalTable <- renderTable({
    if(is.null(data())){return()}
    data()
  })


  ## info
  # output$info <- renderText(
  #   {
  #     paste("Data formate:",DataF(),"\n",
  #           "Term number:",num(),"\n",
  #           "Project name:",name(),"\n",
  #           "Color min:",colormin(),"\n",
  #           "Color max:",colormax(),"\n",
  #           "Gtype:",Gtype(),sep = "")
  #   }
  # )
  ## bubble plot
  p <- reactive({
    if(is.null(data())){return()}
    if (DataF() == "TBtools_GO") {
      GOBubblePlot(expfile = data(),
                   num = num(),
                   name = name(),
                   colormin = colormin(),
                   colormax = colormax(),
                   Gtype = Gtype())
    } else if (DataF() == "TBtools_KEGG") {
      KEGGBubblePlot(expfile = data(),
                     num = num(),
                     name = name(),
                     colormin = colormin(),
                     colormax = colormax(),
                     Gtype = Gtype())
    } else if (DataF() == "customized") {
      customBubblePlot(expfile = data(),
                       num = num(),
                       name = name(),
                       colormin = colormin(),
                       colormax = colormax(),
                       Gtype = Gtype())
    } else if (DataF() == "GSEA") {
      GSEAplot(expfile = data(),
               num = num(),
               name = name(),
               colormin = colormin(),
               colormax = colormax(),
               Gtype = Gtype())
    }
  })
  ## exhibit plot
  output$finalplot <- renderPlot({
    p()
  })
  ## download
  output$down <- downloadHandler(
    filename = function(){
      paste0(name(),"Bubble_plot.pdf")
    },
    content = function(file){
      ggsave(file,plot = p(),width = w(),height = h(),dpi = 300)
    }
  )
  ## table sets
  output$tb <- renderUI({
    if(is.null(data()))
      h5("Please upload your enrichment result")
    else
      tabsetPanel(tabPanel(title = "Table",tableOutput("finalTable")),
#                  tabPanel(title = "Information",tableOutput("info")),
                  tabPanel(title = "plot",
                           plotOutput("finalplot"),
                           downloadButton("down","Download the plot in PDF format")))
})
  ## close window
  observeEvent(input$close, {
    js$closeWindow()
    stopApp()
  })
}) 
## run shiny app
shinyApp(ui = ui,server = server)