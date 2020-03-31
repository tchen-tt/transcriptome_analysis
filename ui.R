library(shiny)
tagList(
  shinythemes::themeSelector(),
  navbarPage(
    "Data Analysis",
    tabPanel("Expression",
             sidebarPanel(
               h4("INPUT FILE"),
               fileInput("file", "File input:"),
               selectInput("species","Species",
                           choices=c("Humn","Zebrafish", "Mouse"),
                           selected="Zebrafish"),
               selectInput("control","Control",
                           choices=c(1:20),
                           selected=c(2:4),
                           multiple=TRUE
               ),
               selectInput("treat","Treat",
                           choices=c(1:20),
                           selected=c(5:7),
                           multiple=TRUE),
               tags$small(paste0("Note: Control and Treat values",
                                " are column in your put data,",
                                " and DEseq2 analysis result is Treat/control.")),
               br(),
               br(),
               h4("DESeq2"),
               
               sliderInput("log2FoldChange","|log2FoldChange|:",
                           min=0.5, max=10,
                           value=1),
               selectizeInput("P.adj","P.adj:",
                              choices=c(0.001, 0.01, 0.05),
                              selected=0.05)
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("ReadCount",
                          DT::dataTableOutput('readcount'),
                          br(),
                          br(),
                          h4("Correlation and cluster"),
                          br(),
                          fluidRow(
                            column(6,
                                   plotOutput("xishu", width = "370px", height = "370px")),
                            column(6,
                                   plotlyOutput("pca",width = "410px", height = "360px"))
                                   
                          )#,
                          #plotlyOutput("pca", width = "400px", height = "400px")
                          
                 ),
                 tabPanel("DESeq2.normalzed",
                          DT::dataTableOutput('normalzed'),
                          plotlyOutput("huoshan", width = "100%", 
                                       height = "400px", inline = FALSE)
                 ),
                 tabPanel("FPKM", 
                          DT::dataTableOutput('fpkm')
                 )
               )
             )
            
    ),
    tabPanel("GO",
             sidebarPanel(
               radioButtons("annotation_term", label="Choose gene list",
                            choices = list("diff"="diff", "up"="up", "down"="down"),
                            selected="diff")
             ),
             mainPanel(
               tags$h4("Go analysis result"),
               DT::dataTableOutput("GOresult")
             )
             ),
    tabPanel("KEGG","This Panel is intetionally KEGG Analysis")
  )
)