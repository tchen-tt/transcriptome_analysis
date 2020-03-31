library(shiny)
function(input, output){
  #exprReadcount <- read.csv(input$file[1,4], header=T, stringsAsFactors=F)
  # data
 
  #if(!exists("input$file")){
  #  print("dongiyna")
  #  source("./formalue.R")
  #  readcount <- rrsult$readcount
  #  diff <- rrsult$diff %>% na.omit
  #  normalized <- rrsult$normalized
  #
  filePath <- reactive({
    if(!is.null(input$file)){
      input$file
    }
  })
  
  control_treat_col <- reactive({
    list(control=as.numeric(input$control),
         treat=as.numeric(input$treat),
         species=input$species)
  })
  geneExpr <- reactive({
    if(!is.null(filePath()$datapath)){
      readcount <- read.csv(filePath()$datapath,
                            header=TRUE, stringsAsFactors=FALSE)
      samples_name <- names(readcount)
      group <- list(control=samples_name[control_treat_col()$control],
                    treat=samples_name[control_treat_col()$treat])
      da <- list("readcount"=readcount[, c(1, control_treat_col()$control, control_treat_col()$treat)],
                 "group"=group)
      
    }
  })
  
  inputfile_pre <- reactive({
    if(!is.null(geneExpr()$readcount)){
      rrsult <- DESeq_s(geneExpr()$readcount, geneExpr()$group)
    } else {
      rrsult
    }
      
  })
  inputfile <- reactive({
    rrsult <- inputfile_pre()
    P.adjs = as.numeric(input$P.adj)
    log2FoldChanges = input$log2FoldChange
    diff = rrsult$diff
    rrsult$diff <- diff %>% na.omit %>%
      dplyr::mutate(dfff1=ifelse( (abs(log2FoldChange)>=log2FoldChanges)&(padj<P.adjs), 3, 0),
               dfff2=dfff1 + sign(log2FoldChange),
               Type=factor(dfff1*dfff2, labels = c("nodiff", "down", "up"))) %>% 
      dplyr::select(-dfff1, -dfff2)
    rrsult
  })
  

  
  output$readcount <- DT::renderDataTable(
    DT::datatable(inputfile()$readcount, options = list(pageLength = 5))
  )
  output$normalzed <- DT::renderDataTable(
    DT::datatable(inputfile()$diff, options = list(pageLength = 5))
  )
  output$fpkm <- DT::renderDataTable(
    DT::datatable(inputfile()$normalized, options = list(pageLength = 5))
  )
  
  
  # picture
  # one xianguanxishu picture
  output$xishu <- renderPlot({
   file = inputfile()$readcount
   file %>% dplyr::select(2:ncol(file)) %>%
     `rownames<-`(file[,1]) %>% 
     as.matrix %>% cor %>%
     heatmap(scale="none", margins = c(7,7))
   })
  
  # two pca plotly
  output$pca <- renderPlotly({
    file = inputfile()$readcount
    cts = file %>% 
      dplyr::select(2:ncol(file)) %>%
      `rownames<-`(file[,1]) %>% 
      as.matrix
    cts <- cts[rowMeans(cts)>1,]
    testing <- prcomp(cts, center = TRUE, retx = TRUE)
    locat <- testing$rotation %>% as.data.frame %>% 
      dplyr::mutate(samplename=rownames(testing$rotation),
                    class = rep(c("Control","Treat"), 
                                c(length(control_treat_col()$control),
                                  length(control_treat_col()$treat))))
    varison <- summary(testing)
    P <- locat %>% ggplot(aes(PC1, PC2, fill=class)) + 
      geom_point(aes(txt=samplename)) + theme_bw() + 
      theme(panel.grid=element_blank()) +
      xlab(paste0("PC1 ", "(", varison$importance[2,1]*100, "%)")) +
      ylab(paste0("PC2 ", "(", varison$importance[2,2]*100, "%)"))
    ggplotly(P)
  })
  
  
  
  # three huoshantu
  output$huoshan <- renderPlotly({
    #cat(filePath()$datapath)
#    P.adjs = as.numeric(input$P.adj)
#    log2FoldChanges = input$log2FoldChange
#    diff = inputfile()$diff
    
#    testing <- diff %>% na.omit %>%
#      mutate(dfff1=ifelse( (abs(log2FoldChange)>=log2FoldChanges)&(padj<P.adjs), 3, 0),
#             dfff2=dfff1 + sign(log2FoldChange),
#             Type=factor(dfff1*dfff2, labels = c("nodiff", "down", "up"))) %>% 
#      select(-dfff1, -dfff2)
    P <- ggplot(inputfile()$diff, aes(log2FoldChange, -log10(pvalue), colour=Type)) + geom_point(aes(text=GeneID), data=inputfile()$diff) + theme_classic()
    ggplotly(P, tooltip = c("text"))
  })
}


