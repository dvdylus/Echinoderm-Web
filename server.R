library(shiny)
source("helpers.R")
#brief <- read.csv("Data/Brief.csv",header=T)



load("Data/Brief4_27_v5.RData")
all_images <- dir("www/WMISH/")
tablemax <- 30000
seqmax   <- 100
plotmax  <- 1000

#print(all_images)

shinyServer(function(input, output) {
  
  # a large table, reative to input$show_vars
  output$mytable1 = renderDataTable({
    library(ggplot2)
    brief[, input$show_vars, drop = FALSE]
  })

  infotbl <- reactive({ 
    
    infotbl <- data.frame()
    
    
    
    # decide searchfield and searchterms
    searchfield <- input$searchfield
	
    if (searchfield == "ID") { # Searches by specific IDs
      searchterms <- splitlines(input$namelist)
      #exprterms <- NULL
	 }else if (searchfield == "Cluster_ID") { # Searches by Expr Cluster
      searchterms <- splitlines(input$clusterlist)
      #exprterms <- NULL
    }else if (searchfield == "Annotation") { # Searches by Annotation
      searchterms <- splitlines(input$annotlist)
      #exprterms <- NULL
    }else if (searchfield == "Class.L1") { # Searches by Function Class
      searchterms  <- brief$ID[which(brief$Class.L1==input$Class.L1)]
	  exprterms <- splitlines(input$exprlist)
      #searchfield <- "ID"
    }
    
    # decide to maybe only use orhtologs or homologs
    relate <- input$relate
    
    infotbl <- exp.search(brief, searchterms, searchfield, exprterms, relate)
    rownames(infotbl) <- NULL
    
    return(infotbl)
    
  })
  
  
  
  output$desctbl <- renderTable({
    infotbl <- infotbl()
    infotbl.n <- nrow(infotbl)
    if (infotbl.n<1) { 
      return() 
    }
    if (infotbl.n>tablemax) { 
      infotbl <- data.frame(matrix(ncol=1, nrow=1))
      infotbl[1,1] <- paste("Over limit:", tablemax)
      colnames(infotbl) <- c("message")
      return(infotbl)
    }
    
    desctbl <- subset(infotbl, select=c("ID", "Annotation", "Relation","Class.L1", "Cluster_ID","Expr_Cluster"))
    
    return(desctbl)
    
  }, sanitize.text.function = function(x){x}, include.rownames=F)
  
  output$datatbl <- renderTable({
    infotbl <- infotbl()
    infotbl.n <- nrow(infotbl)
    if (infotbl.n<1) { 
      return() 
    }
    if (infotbl.n>tablemax) { 
      infotbl <- data.frame(matrix(ncol=1, nrow=1))
      infotbl[1,1] <- paste("Over limit:", tablemax)
      colnames(infotbl) <- c("message")
      return(infotbl)
    }
    
    exprlist <- infotbl$Cluster_ID
    exprtbl <- retrieveexpr(expr,exprlist) # load dataframe expr to retrive expr.data
    tblcols    <- colnames(exprtbl)
    timepoints <- tblcols[grepl("exp_", tblcols)]
    tblcols    <- c("Cluster_ID", "Annotation", timepoints[1:4],"fc","Expr_Cluster")
    datatbl    <- subset(exprtbl, select=tblcols)
    tblcols    <- gsub("exp_", "", tblcols)
    colnames(datatbl) <- tblcols
    
    return(datatbl)
    
  }, include.rownames=F, digits=2, format.args=list(big.mark=","))
  
  outwidth <- reactive({
    if (input$plottype == 'plotline') {
      outwidth <- 750
      
    } else if (input$plottype == 'plotheat') {
      
      infotbl <- infotbl()
      if (nrow(infotbl)<1)       { return(0) }
      if (nrow(infotbl)>plotmax) { return(0) }
      infotbl$fullname <- paste(infotbl$Cluster_ID, infotbl$Annotation, sep=" = ")
      longest <- max(nchar(as.character(infotbl$fullname)))
      outwidth <- 500 + longest*8
    } else if (input$plottype == 'diffplot'){
      outwidth <- 750
    }
    return(outwidth)
    
  })
  
  outheight <- reactive({
    
    if (input$plottype == 'plotline') {
      outheight <- 400
      
    } else if (input$plottype == 'plotheat') {
      infotbl <- infotbl()
      infotbl.n <- nrow(infotbl)   
      if (infotbl.n<1)           { return(0) }
      if (nrow(infotbl)>plotmax) { return(0) }
      outheight  <- infotbl.n * 25
      if (outheight < 300) {
        outheight <- 300
      } else if (outheight>20000) {
        outheight <- 20000
      }
    }  else if (input$plottype == 'diffplot'){
      outheight <- 400
    }  
    return(outheight)
  })
  
  output$profilePlot <- renderPlot({
    
    infotbl <- infotbl()
    infotbl.n <- nrow(infotbl)
    if (infotbl.n<1) { 
      return() 
    }
    if (infotbl.n>plotmax) { 
      return()
    }
    
    if (input$plottype == 'plotline') {
      if ( input$facet=='single' & infotbl.n<=100 ) {
        combine <- 'single'
      } else if (input$facet=='cluster') {
        combine <- 'cluster'
      } else {
        combine <- 'all'
      }            
      
      legend <- input$legend
      if ( legend==TRUE & infotbl.n>10 ) {
        legend <- FALSE
      }   
      
      exprlist <- infotbl$Cluster_ID
      exprtbl <- retrieveexpr(expr,exprlist) # load dataframe expr to retrive expr.data
      myplot <- exp.plot.line(exprtbl, combine=combine, ytrans=input$ytrans, legend=legend)
      print(myplot)
      
    } else if (input$plottype == 'plotheat') {
      exprlist <- infotbl$Cluster_ID
      exprtbl <- retrieveexpr(expr,exprlist) # load dataframe expr to retrive expr.data
      myplot    <- exp.plot.heat(exprtbl, unit="tpm")                
      print(myplot)
      
    } else if (input$plottype == 'diffplot'){
      if ( input$restrict=='detect') {
        combine <- 'Detectable'
      } else if (input$restrict=='gfold') {
        combine <- 'Gfold'
      } else {
        combine <- 'all'
      }            
      
      
      exprlist <- infotbl$Cluster_ID
      exprtbl <- retrieveexpr(differ,exprlist) # load dataframe expr to retrive expr.data
      myplot <- exp.plot.diff(exprtbl,combine)
      print(myplot)
    }
    
  }, width=outwidth, height=outheight)
  
  output$wmishtbl <- renderTable({
    
    infotbl <- infotbl()
    #print(infotbl[1,])
    infotbl.n <- nrow(infotbl)
    if (infotbl.n<1) { 
      return() 
    }
    
    all_images <- dir("www/WMISH/")
    
    wmishtable <- create_wmish_table_vertical(infotbl, all_images)
    #print(wmishtable)
    
#     dat <- data.frame(
#       country = c('Bl', 'MBl','G'),
#       flag = c('<img src=WMISH/G-Dri.jpg></img>',
#                paste(im_folder,"MBl-Dri.tif",sep=""),
#                '<img src=WMISH/G-Dri.jpg></img>')
#     )
#     dat
  return(wmishtable)
    
  }, sanitize.text.function = function(x){x})

  #, sanitize.text.function = function(x){x}
  
  output$mrna <- renderText({
    infotbl <- infotbl()
    infotbl.n <- nrow(infotbl)
    if (infotbl.n<1) { return() }
    if (infotbl.n>seqmax) { return(paste("Over limit:", seqmax)) }
    
    idlist <- infotbl$ID
    retrieveseq(sequences,idlist, 'Sequence')
  })
  
  output$downloadData <- downloadHandler(
    filename = "GeneTable.txt",
    content = function(file) {
      write.table(infotbl(), file, row.names=F, sep="\t")
    })
  
})