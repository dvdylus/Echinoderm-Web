library(shiny)
library(ggplot2)  # for the diamonds dataset

menuitem <- c(  "ID"               = "ID",
				"Cluster_ID"			= "Cluster_ID",
                "Annotation"         = "Annotation",
                "Function Class"     = "Class.L1"
)  

classitem <- c("Adhesion", "Apoptosis", "Biomineralization", "CalciumToolkit", "CellCycle", "Cytoskeleton", "DNAReplication", "Defensome", "EggActivation", "GPCRRhodopsin", "GTPase", "GermLineDeterminant", "Histone", "Immunity", "Kinase", "Metabolism", "Metalloprotease", "Nervous", "Oogenesis", "Phosphatase", "Signaling", "TF", "TranslationFactor", "ZNF")


shinyUI(pageWithSidebar(
  headerPanel("Transcriptome of Amphiura filiformis"),
    
  sidebarPanel(
    
    # construct HTML by 'tags'. design the width of sidebar panel
    tags$head(
      tags$style(type='text/css', ".span4 { max-width: 200px; }")
    ),
    
	# Left handside to select orthology
	wellPanel(

		selectInput("searchfield", "Genes", 
		choices = menuitem,
		selected = "ID"
		),

		conditionalPanel(
			condition = "input.searchfield == 'ID'",
			tags$textarea(id="namelist", rows=3, "AfiCDS.id1.tr5535") 
			),

		conditionalPanel(
			condition = "input.searchfield == 'Cluster_ID'",
			tags$textarea(id="clusterlist", rows=3, "Cluster-0.0") 
			),

		conditionalPanel(
			condition = "input.searchfield == 'Annotation'",
			tags$textarea(id="annotlist", rows=3, "Tgif") 
			),

		conditionalPanel(
			condition = "input.searchfield == 'Class.L1'",
			
			selectInput("Class.L1", "", 
				choices = classitem,
				selected = "TF"
				),
			h5("Expression Cluster"),
			selected="TF",  
			#textInput("exprlist", "Expression Cluster", "1") 
			tags$textarea(id="exprlist", rows=3, "1")
				)
		),

    # Left handside to select orthology
	wellPanel(
      
      #checkboxGroupInput("relate","Relation",c("Orhtolog"="ortho","Homolog"="homo"),inline=T)
      selectInput("relate", "Relation", 
                  choices = c("All"="all","Ortholog"="ortho","Homolog"="homo"),
                  selected = "all"
      )

      ),

    # Left handside PLOT panel taken from Qiang Tu
    wellPanel(     
      helpText("Plot"),
      radioButtons("plottype", "",
                   choices = c("Heatmap"  = "plotheat",
                               "Lineplot" = "plotline",
                               "Diffplot" = "diffplot"
                   ),
                   "plotheat"         
      ),
      
      conditionalPanel(
        condition = "input.plottype == 'plotline'",
        wellPanel(     
          checkboxInput("legend", "Show legend (<10)",  value=FALSE),
          radioButtons("ytrans", "Transform Y:",
                       choices = c("None"    = "",
                                   "Log10"   = "log",
                                   "Percent" = "percent"
                       ),
                       ""
          ),
          radioButtons("facet", "Group profiles by:",
                       choices = c("All"     = "all",
                                   "Cluster" = "cluster",
                                   "Gene (<100)"  = "single"),
                       "all"                    
          )
        )
      ),
      conditionalPanel(
        condition = "input.plottype == 'diffplot'",
        wellPanel(     
          radioButtons("restrict", "Restrict:",
                       choices = c("None"    = "",
                                   "Detectable"   = "detect",
                                   "Gfold" = "gfold"
                       ),
                       ""
          )
        )
      )
    ),
    # Left handside DOWNLOAD panel
    wellPanel(
      helpText("Table"),
      checkboxInput("only2col", "Briefer",  value=FALSE),
      downloadButton('downloadData', 'Download Full Table')
    )
    
    
  ),
  
    
  mainPanel(
    tabsetPanel(
      tabPanel('Brief',tableOutput("desctbl")),
      tabPanel("Data", tableOutput("datatbl")),
      tabPanel("Plot",  plotOutput("profilePlot", width='auto', height='auto')),
      tabPanel("WMISH", tableOutput("wmishtbl")),
      tabPanel("CDS",  verbatimTextOutput ("mrna"))
    )
  )
))