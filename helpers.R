library(ggplot2)
library(scales)
library(reshape2)
library(stringr)
source("ggplot_pLIB.R")

all_images <- dir("www/WMISH/")

#---------------------------------
# Function that based on search field will create sub-data
#
#	expdata
#	termlist
#	field
#	exprlist
#	relation
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
exp.search <- function(expdata, termlist, field, exprlist=NULL, relation="all") {
	result <- data.frame()

	if (field=="ID") {
		result <- subset(expdata, ID %in% termlist)
	} else if (field=="Annotation"){
		k.select <- which(grepl(termlist[1],expdata$Annotation, ignore.case = T))
		
		if(length(termlist)>1){
			for(i in 2:length(termlist)){
				k.select<- c(k.select,which(grepl(termlist[i],expdata$Annotation, ignore.case = T)))
			}
		}
			
		termlist.new <- as.character(unique(expdata$Annotation[k.select]))
		#print(termlist.new)
		result <- subset(expdata, Annotation %in% termlist.new)
	} else if (field=="Cluster_ID"){
		result <- subset(expdata, Cluster_ID %in% termlist)
	} else if (field=="Class.L1"){
		tmp <- subset(expdata, ID %in% termlist)
		# This is a very ugly way to ensure that if textarea is empty all expr_clusters are shown or if any character is introduced that it switches to all -> SUPER UGLY but works
		test_numeric <- as.numeric(exprlist)
		#print(test_numeric)
		if(length(test_numeric)==0){ # test whether there is no entry in text field
			test <- NULL
			#print("1A")
		} else if (length(which(is.na(test_numeric)))>0){ # test whether entry in textfield is not a number
			test <- NULL
			#print("1A")
		}else{
			test <- exprlist
			#print("1B")
		}
		if (!is.null(test)){
			result <- subset(tmp, Expr_Cluster %in% exprlist)
			#print("2A")				
		} else{
			result <- tmp
			#print("2B")
		}
	}

			if (nrow(result)>0) {
				# remove possible duplicates from multiple queries
				result <- subset(result, !duplicated(ID))

				# remove 'Sp-' prefix. might try to remove 'none' also
				# result$name <- gsub("Sp-", "", result$name)
			}

			if (relation=="ortho"){
				result <- result[which(result$Relation=="Ortholog"),]
				} else if(relation=="homo"){
					result <- result[which(result$Relation=="Homolog"),]
				}

				return(result)
			}


#---------------------------------
# Function that splitline????????
#
#	gene_list: list of genes taken from annotation by organism x
#	fc: fclust object obtained through fuzzy clustering
#	cluster_to_id: data.frame containing a map from id to expression clusters (CORSET output)
#	annotation: original annotation file as csv output of titus eel_pond algorithm
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
splitlines <- function(intext) {
  inputlines <- strsplit(intext, "\n")[[1]]
  inputlines <- inputlines[!grepl("^\\s*$", inputlines, perl=T)]  # ignore blank lines
  inputlines <- gsub("^\\s+", "", inputlines, perl=T)             # remove leading spaces
  inputlines <- gsub("\\s+$", "", inputlines, perl=T)
  return(inputlines)
}

#---------------------------------
# Function that retrieves sequences based on specified IDs
#
#	seqdb:
#	idlist:
#	seqtype:
#
#	OUTPUT: 
#---------------------------------
retrieveseq <- function(seqdb,idlist, seqtype) {
  allseq     <- subset(seqdb, ID %in% idlist, select=c("ID","Annotation","Cluster_ID", seqtype))
  
  returnstr <- ""
  for(i in 1:nrow(allseq)) {
    seqdesc   <- paste(">", allseq[i,"ID"], "-", allseq[i,"Cluster_ID"], " ", allseq[i,"Annotation"], sep="")
    seqstr    <- gsub("(.{1,60})", "\\1\n", allseq[i,seqtype])
    returnstr <- paste(returnstr, seqdesc, "\n", seqstr, "\n", sep="")
  }    
  return(returnstr)
}

#---------------------------------
# Function that retrieves expression data based on specified expression clusters
#
#	seqdb:
#	clusterlist:
#
#	OUTPUT: 
#---------------------------------
retrieveexpr <- function(seqdb,clusterlist) {
  allexpr     <- subset(seqdb, Cluster_ID %in% clusterlist) 
  return(allexpr)
}

#---------------------------------
# Creates a lineplot for selected genes
#
#	indata:
#	combine:
#	ytrans:
#	scales:
#	title:
#	legend:
#
#	OUTPUT: depending of parameters it can produce different types of lineplots as output
#---------------------------------
exp.plot.line <- function(indata, ylab="", combine="all", ytrans="", scales="fixed", title="", legend=F) {
  
  # input: Cluster_ID x timepoints dataframe
  # ytrans: percent, log
  # optional columns: Cluster_ID (otherwise rownames), cluster (otherwise "all")
  # return: a plot obj
  
  expdata <- data.frame(indata)
  tblcols    <- colnames(expdata)
  timepoints <- tblcols[grepl("exp_", tblcols)]
  points <- timepoints[1:4]
  tblcols    <- c("Cluster_ID", "Annotation","Expr_Cluster", points)
  datatbl    <- subset(expdata, select=tblcols)
  plotdata <- datatbl
  #points  <- colnames(expdata)[grep("exp_", colnames(expdata))]
  #expdata$Cluster_ID <- expdata$fullname
  #plotdata <- expdata[,c("ID", "cluster", points)]
  
  if (ytrans=="percent") {
    maxvalue          <- apply(plotdata[,points], 1, max)
    plotdata[,points] <- plotdata[,points] / maxvalue
  }
  
  #colnames(plotdata) <- simpname(colnames(plotdata))   # hrxx
  
  # melt
  plotdata.m <- melt(plotdata, id=c("Cluster_ID", "Annotation","Expr_Cluster"), variable.name="time")
  plotdata.m$time <- gsub("exp_", "", plotdata.m$time)
  plotdata.m$time[which(plotdata.m$time=="9hr")]<-"09hr"
  #plotdata.m$Cluster_ID <- gsub("^\\d\\d\\d ", "", plotdata.m$Cluster_ID)
  
  plot <- ggplot(data=plotdata.m, aes(x=time, y=value, group=Cluster_ID)) +
    labs(x="Time", title=title) +
    theme(axis.text.x = element_text(color='black'),
          axis.text.y = element_text(color='black'))   
  
  plot <- plot + theme_Publication(base_size=14)
  # decide legend
  if (legend==F) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # decide combine type: all, single or cluster
  if (combine == "all") {
    plot <- plot + geom_line(aes(color=Cluster_ID))
    
  } else if (combine == "single") {
    plot <- plot + geom_line(aes(color=Cluster_ID)) + facet_wrap(~Cluster_ID, scales=scales)
    
  } else if (combine == "cluster") {
    plot <- plot + geom_line(aes(color=Expr_Cluster)) + facet_wrap(~Expr_Cluster)
    
  } else {
    print("unknown combine type")
    return()
  }
  
  # Y transform
  values <- plotdata.m$value
  ymax <- max(values)*1.1
  ymin <- min(values[values>0])*0.9
  if (ytrans=="log") {
    plot <- plot + scale_y_log10(limits=c(ymin, ymax)) + labs(y="Transcripts/embryo")
  } else if (ytrans=="percent") {
    plot <- plot + scale_y_continuous(labels = percent, limits=c(0,1)) + labs(y="Percent of max")
  } else {
    plot <- plot + scale_y_continuous(limits=c(0, ymax)) + labs(y="Transcripts/embryo")
  }
  
  return(plot)
  
}

#---------------------------------
# Creates heatplot of time-course data
#
#	indata:
#	title:
#	palette:
#	uint:
#
#	OUTPUT: heatplot
#---------------------------------
exp.plot.heat <- function (indata, title="", palette="YlOrRd", basesize=18, unit="tpm") {
  
  if (unit == "tpm") {
    mybins      <- c(0,     5,         10,         20,         40,         80,       10000)
    mybinlabels <- c("A. <5", "B. 5-10", "C. 10-20", "D. 20-40", "E. 40-80", "F. >80")
  } else {  # assume cpe
    mybins      <- c(0,      250,          500,           1000,          2000,           4000,       1000000)
    mybinlabels <- c("A. <250", "B. 250-500", "C. 0.5-1K", "D. 1-2K", "E. 2-4K", "F. >4K")
  }

  
  expdata <- data.frame(indata)
  tblcols    <- colnames(expdata)
  timepoints <- tblcols[grepl("exp_", tblcols)]
  points <- timepoints[1:4]
  tblcols    <- c("Cluster_ID","Annotation", points)
  datatbl    <- subset(expdata, select=tblcols)
  plotdata <- datatbl
  
  # order by reverse fullname as on webpage, the head is seen first;
  # then remove the leading cluster IDs
  if(length(plotdata$Cluster_ID)>1){
    plotdata          <- plotdata[order(plotdata$Cluster_ID, decreasing=T), ]   
  }
  
  plotdata$fullname<-paste(plotdata$Cluster_ID, plotdata$Annotation, sep=" = ")
  #plotdata$fullname[which(plotdata$Annotation=='')]<-plotdata$Cluster_ID[which(plotdata$Annotation=='')]
  
  # melt
  plotdata.m <- melt(plotdata, id=c("Cluster_ID", "Annotation", "fullname"), variable.name="time")
  plotdata.m$Level <- cut(plotdata.m$value, breaks=mybins, labels=mybinlabels, include.lowest=T)
  plotdata.m$time <- gsub("exp_", "", plotdata.m$time)
  plotdata.m$time[which(plotdata.m$time=="9hr")]<-"09hr"
  
  plot <- ggplot(plotdata.m, aes(time, fullname, fill=Level)) + geom_tile() +
    scale_fill_brewer(palette=palette, name="Transcripts/Embryo") + 
    scale_y_discrete(expand = c(0, 1)) +
    theme_grey(base_size = basesize) +
    theme(axis.text.x = element_text(color='black'),
          axis.text.y = element_text(color='black')) +
    labs(title=title, x="Time", y="")
  plot <- plot + theme_Publication(base_size=14)
  return(plot)
}

#---------------------------------
# Creates differential plot for selected expression clusters
#
#	indata
#	combine
#	palette
#	basesize
#
#	OUTPUT: bar plot with differential genes
#---------------------------------
exp.plot.diff <- function(indata, combine="all", palette="YlOrRd", basesize=18) {
  
  # exception plot handles
  # set all NAs to zero
  expdata <- data.frame(indata)
  tblcols    <- colnames(expdata)
  timepoints <- tblcols[grepl("fc_", tblcols)]
  points <- timepoints[1:2]
  tblcols    <- c("Cluster_ID","Annotation","Effect", points)
  datatbl    <- subset(expdata, select=tblcols)

  
  if(combine=="Detectable"){
    plotdata <- datatbl[which(datatbl$Effect=="Detectable" | datatbl$Effect=="Gfold"),]
#     k.annot <- which(plotdata$Annotation!="")
#     k.not.annot <- which(plotdata$Annotation=="")
#     plotdata$fullname<-""
#     plotdata$fullname[k.annot]<-paste(plotdata$Cluster_ID[k.annot], plotdata$Annotation[k.annot], sep=" = ")
#     plotdata$fullname[k.not.annot]<-as.character(plotdata$Cluster_ID[k.not.annot])
  }else if (combine=="Gfold") {
    plotdata <- datatbl[which(datatbl$Effect=="Gfold"),]
#     k.annot <- which(plotdata$Annotation!="")
#     k.not.annot <- which(plotdata$Annotation=="")
#     plotdata$fullname<-""
#     plotdata$fullname[k.annot]<-paste(plotdata$Cluster_ID[k.annot], plotdata$Annotation[k.annot], sep=" = ")
#     plotdata$fullname[k.not.annot]<-as.character(plotdata$Cluster_ID[k.not.annot])
  }else {
    plotdata <- datatbl
#     k.annot <- which(plotdata$Annotation!="")
#     k.not.annot <- which(plotdata$Annotation=="")
#     plotdata$fullname<-""
#     plotdata$fullname[k.annot]<-paste(plotdata$Cluster_ID[k.annot], plotdata$Annotation[k.annot], sep=" = ")
#     plotdata$fullname[k.not.annot]<-as.character(plotdata$Cluster_ID[k.not.annot])
  }
  
  # order by reverse fullname as on webpage, the head is seen first;
  # then remove the leading cluster IDs
  if(length(plotdata$Cluster_ID)>1){
    plotdata          <- plotdata[order(plotdata$Cluster_ID, decreasing=T), ]   
  }
  
  #plotdata$Annotation[which(plotdata$Annotation=='')]<-plotdata$Cluster_ID[which(plotdata$Annotation=='')]
  # melt
  plotdata.m <- melt(plotdata, id=c("Cluster_ID","Effect","Annotation"), variable.name="differ")
  
  if (length(which(plotdata.m$value==-Inf))>0){
    plotdata.m$value[which(plotdata.m$value==-Inf)]=-10
  }
  
  if (length(which(plotdata.m$value==Inf))>0){
    plotdata.m$value[which(plotdata.m$value==Inf)]=10
  }
  
  if (length(which(is.na(plotdata.m$value)))>0){
    plotdata.m$value[is.na(plotdata.m$value)]=0
  }
 
  
  plotdata.m$value <- -1*plotdata.m$value
  
  plotdata.m$text <- ""
  plotdata.m$text[which(plotdata.m$differ=="fc_dmso")] <- as.character(plotdata.m$Annotation[which(plotdata.m$differ=="fc_dmso")])
  
  plot <- ggplot(plotdata.m,aes(Cluster_ID,value,fill=differ)) + 
    geom_bar(stat="identity", position="dodge") +
    #geom_text(aes(y=10,label=text), position="dodge",hjust=1, vjust=0, size=5, angle= 90 ,colour="black")+
    ylim(-10,10)+
    scale_fill_brewer(palette=palette, name="Differential") + 
    theme_Publication(base_size=14) +
    theme(axis.text.x = element_text(color='black',angle = 90, hjust = 1),
          axis.text.y = element_text(color='black')) + 
    geom_hline(yintercept=1.6, colour="#990000", linetype="dashed")+ geom_hline(yintercept=-1.6, colour="#990000", linetype="dashed") +
    ylab("log2(FC)") +
    xlab("Cluster")

  
  return(plot)
  
  
}

#---------------------------------
# Checks whether images for selected gene exist in the right folder and selects folder
#
#	indata: 
#  	allimagefolder:
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
check_wmish <- function(indata,allimagefolder) {
	
	#just_gene_names <- unique(gsub(".*Sp-", "", indata$Annotation, perl=T)) # find the right gene
  just_gene_names <- unique(indata$ID) # in the case that all folder were called with their ID!!!!
	k.folder.tmp <- match(just_gene_names,allimagefolder) # check which of the selected genes is in a folder and return their idx
	k.folder <- k.folder.tmp[which(!is.na(k.folder.tmp))]
  k.folder.annot <- indata$Annotation[which(indata$ID %in% allimagefolder)]
  k.folder.annot <- gsub(".*-Sp-","",k.folder.annot)
  k.folder.annot <- gsub("-WHL.*","",k.folder.annot)
	if(length(k.folder)<1){
		k.folder <- NA
	}
  if(length(k.folder.annot)<1){
    k.folder.annot <- NA
  }
	#print(k.folder)
  k.df <- data.frame(k.folder=k.folder, k.folder.annot=k.folder.annot)
	return(k.df)

}


#---------------------------------
# Checks whether WMISH was performed on currently selected list of genes!
#
#  gene_list: list of genes taken from annotation by organism x
#  fc: fclust object obtained through fuzzy clustering
#  cluster_to_id: data.frame containing a map from id to expression clusters (CORSET output)
#	annotation: original annotation file as csv output of titus eel_pond algorithm
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
get_wmish_filename <- function(gene,stage) {
	allf <- dir("www",recursive=T)
	allf.sub <- allf[which(grepl(gene,allf,fixed=T)==T)]
	fname.tmp <- allf.sub[which(grepl(stage,allf.sub,fixed=T)==T)]
	#print(fname.tmp)
	if(length(fname.tmp)>1){
		fname<-fname.tmp[1]
		}else{
			fname <- fname.tmp
		}
		#print(fname)
		if(length(fname)!=0){
			out.fname.tmp <- paste('<img src=',fname,sep="")
			out.fname <- paste(out.fname.tmp,' width="150" height="150"></img>',sep="")
			}else{
				out.fname=""
			}
			#print(out.fname)
			return(out.fname)
		}


#---------------------------------
# Find gene list in fuzzy clustering, for example transcription factors
#
#  gene_list: list of genes taken from annotation by organism x
#	fc: fclust object obtained through fuzzy clustering
#	cluster_to_id: data.frame containing a map from id to expression clusters (CORSET output)
#	annotation: original annotation file as csv output of titus eel_pond algorithm
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
create_wmish_table_horizontal <- function(indata,allimagefolders) {

	stages <- c("Cl","EBl","Bl","Bl-VV","MBl","G","Pr","Pl")

  k.df <- check_wmish(indata,allimagefolders)
	k.folder <- k.df$k.folder
	#print(k.folder)
	if(is.na(k.folder)){
		out.df <- data.frame(Output = "No WMISH performed")
		return(out.df)
		} else {
			out.df <- data.frame(Gene=allimagefolders[k.folder],Annot=k.df$k.folder.annot, Cleavage="",Early_Blastula="",Blastula="",Blastula_VV="",MesBlastula="",Gastrula="",Prism="",Pluteus="",stringsAsFactors=F)

			for(i in 1:length(k.folder)){
				k <- k.folder[i]
				
				out.df$Cleavage[i] <- as.character(get_wmish_filename(allimagefolders[k],stages[1]))
				out.df$Early_Blastula[i] <- as.character(get_wmish_filename(allimagefolders[k],stages[2]))
				out.df$Blastula[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[3]))
				out.df$Blastula_VV[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[4]))
				out.df$MesBlastula[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[5]))
				out.df$Gastrula[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[6]))
				out.df$Prism[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[7]))
				out.df$Pluteus[i] <-as.character(get_wmish_filename(allimagefolders[k],stages[8]))
			}

			return(out.df)
		}


	}
	
#---------------------------------
# Find gene list in fuzzy clustering, for example transcription factors
#
#  gene_list: list of genes taken from annotation by organism x
#	fc: fclust object obtained through fuzzy clustering
#	cluster_to_id: data.frame containing a map from id to expression clusters (CORSET output)
#	annotation: original annotation file as csv output of titus eel_pond algorithm
#
#	OUTPUT: vector with number of occurences for each cluster
#---------------------------------
create_wmish_table_vertical <- function(indata,allimagefolders) {

	stages <- c("Name","Cl","EBl","Bl","Bl-VV","MBl","G","Pr","Pl")

  k.df <- check_wmish(indata,allimagefolders)
	k.folder <- k.df$k.folder
	#print(k.folder)
	if(length(which(is.na(k.folder)))>0){
		out.df <- data.frame(Output = "No WMISH performed")
		return(out.df)
	} else{
		out.df <- data.frame(Stages=stages,stringsAsFactors=F)
		
		print(length(k.folder))
		for(i in 1:length(k.folder)){
			tmp.df <- data.frame(Stages=stages,stringsAsFactors=F)
			k <- k.folder[i]
      tmp.df$Stages[1] <- as.character(k.df$k.folder.annot[i])
			tmp.df$Stages[2] <- as.character(get_wmish_filename(allimagefolders[k],stages[2]))
			tmp.df$Stages[3] <- as.character(get_wmish_filename(allimagefolders[k],stages[3]))
			tmp.df$Stages[4] <- as.character(get_wmish_filename(allimagefolders[k],stages[4]))
			tmp.df$Stages[5] <- as.character(get_wmish_filename(allimagefolders[k],stages[5]))
			tmp.df$Stages[6] <- as.character(get_wmish_filename(allimagefolders[k],stages[6]))
			tmp.df$Stages[7] <- as.character(get_wmish_filename(allimagefolders[k],stages[7]))
			tmp.df$Stages[8] <- as.character(get_wmish_filename(allimagefolders[k],stages[8]))
			tmp.df$Stages[9] <- as.character(get_wmish_filename(allimagefolders[k],stages[9]))
			names(tmp.df)<-as.character(allimagefolders[k])
			#print(tmp.df$Stages)
			out.df <- cbind(out.df,tmp.df)
			#print[out.df[,i]]
		}

		return(out.df)
		}


	}

