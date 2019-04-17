#############################################
### Code to create function "GetGeneList" ###
### Hanna and Riley                       ###
#############################################

GetGeneList = function(Species,savefiles=FALSE,destfile){
  assembly = NULL
  if(missing(Species)){
		stop("ERROR: No species specified")
	}
	if(missing(destfile)){
		stop("ERROR: No path was specified for saving temporary and permanent files.")
	}
	cat("Please be patient, this could take a few minutes.","\n")
	dest = paste(destfile,"feature_table.txt.gz",sep="")
	term = paste(Species,"[orgn] latest_refseq[filter]",sep = "")
	AssemblyInfo = entrez_search(db="assembly",term,use_history = TRUE)
	AssemblySum = entrez_summary("assembly",web_history = AssemblyInfo$web_history)
	URL = paste(extract_from_esummary(AssemblySum,"ftppath_refseq"),"/",extract_from_esummary(AssemblySum,"assemblyaccession"),"_",extract_from_esummary(AssemblySum,"assemblyname"),"_feature_table.txt.gz",sep="")
	cat("Reading in file...","\n")
	   	download.file(URL,dest,cacheOK=TRUE)
	    remove(term,AssemblyInfo,AssemblySum,URL)
	NCBIList<-read.delim(gzfile(dest),header=FALSE,fill=TRUE,skip=1)
	colnames(NCBIList) = c("feature","class","assembly","assembly_unit","seq_type","chromosome","genomic_accession","start","end","strand","product_accession","non-redundant_refseq","related_accession","name","symbol","GeneID","locus_tag","feature_interval_length","product_length","attributes")
	#Section to Clean up Feature Table As Desired
	Assembly = unique(NCBIList[,'assembly',drop=FALSE])
	rownames(Assembly) = c(1:dim(Assembly)[1])
	if(dim(Assembly)[1]>1){
	  cat("Duplicate gene information may be present due to multiple assemblies and feature types.","\n","The following assembly builds are present in this list:","\n")
	  print(Assembly)
	  y = readline("Please choose which ASSEMBLY that you want to prioritize feature information from \n (e.g. 1 for the first assembly listed, 2 for the second, etc.). \n Other duplicate gene information (if any) will be removed from the list. \n ")
	  y = as.numeric(y)
	  if(abs(y) > nrow(Assembly)){
	    stop("ERROR: You specified a number outside the range possible for the assemblies. Please start over.")
	  }
	  ListSubA = subset(NCBIList,assembly == Assembly[y,1])
	}else{
	  cat("The only assembly information in this file is:","\n")
	  print(array(Assembly))
	  ListSubA = NCBIList
	}
	Features = unique(NCBIList[,'feature',drop=FALSE])
	rownames(Features) = c(1:dim(Features)[1])
	if(dim(Features)[1]>1){
	  cat("Duplicate gene information may be present due to multiple feature types.","\n","The following feature types are present in this list:","\n")
	  print(Features)
	  x = readline("Do you want to keep multiple feature type information? y = yes, n = no \n ")
	  if(x != "y"){
	    if(x !="n"){
	      stop("ERROR: You did not answer if you wanted duplicate feature type information removed.","\n","Please start over and enter y for yes or n for no when prompted.")
	    }
    }
	  if(x == "n"){
	    z = readline("Please choose which FEATURE TYPE that you want to prioritize information from \n (e.g. 1 for the first feature listed, 2 for the second, etc.). \n Other duplicate information (if any) will be removed from the list. \n")
	    z = as.numeric(z)
	    if(abs(z) > nrow(Features)){
	      stop("ERROR: You specified a number outside the range possible for the feature types.")
	    }
	    GeneID = array(unique(ListSubA[,"GeneID"]))
	    ListSubAF = ListSubA[which(ListSubA[,"GeneID"]%in%GeneID & ListSubA[,"feature"]==Features[z,1]),]
	    ListSubAF2 = subset(ListSubA, !(ListSubA[,"GeneID"] %in% ListSubAF[,"GeneID"]))
	    ListSubAF = rbind(ListSubAF,ListSubAF2)
	    remove(ListSubA,ListSubAF2)
	  }else{
	    ListSubAF = ListSubA
	    remove(ListSubA)
	  }
	}else{
	  cat("The only feature information in this file is:","\n")
	  print(array(Features))
	  ListSubAF = ListSubA
	  remove(ListSubA)
	}
	Class = unique(ListSubAF[,"class",drop=FALSE])
	rownames(Class) = c(1:dim(Class)[1])
	if(dim(Class)[1]>1){
	  cat("Duplicate gene information may be present due to multiple class types of a given feature.","\n","The following class types are present in the subsetted list:","\n")
	  print(Class)
	  c = readline("Do you want to keep multiple class type information? y = yes, n = no \n ")
	  if(c != "y"){
	    if(c !="n"){
	      stop("ERROR: You did not answer if you wanted duplicate feature type information removed.","\n","Please start over and enter y for yes or n for no when prompted.")
	    }
	  }
	  if(c == "n"){
	    a = readline("Please choose which CLASS TYPE that you want to prioritize information from \n (e.g. 1 for the first feature listed, 2 for the second, etc.). \n Other duplicate information (if any) will be removed from the list. \n")
	    a = as.numeric(a)
	    if(abs(a) > nrow(Class)){
	      stop("ERROR: You specified a number outside the range possible for the class types.")
	    }
	    ListSubAFC = ListSubAF[which(ListSubAF[,"GeneID"]%in%GeneID & ListSubAF[,"class"]==Class[a,1]),]
	    ListSubAFC2 = subset(ListSubAF, !(ListSubAF[,"GeneID"] %in% ListSubAFC[,"GeneID"]))
	    ListSubAFC = rbind(ListSubAFC,ListSubAFC2)
	    remove(ListSubAF,ListSubAFC2)
	  }else{
	    ListSubAFC = ListSubAF
	    remove(ListSubAF)
	  }
	}else{
	  cat("The only class information in this file for the feature(s) specified is:","\n")
	  print(array(Class))
	  ListSubAFC = ListSubAF
	  remove(ListSubAF)
	}
	FeatureList=ListSubAFC[order(ListSubAFC[,"chromosome"],ListSubAFC[,"start"]),] ##Sorting by Chr & Start Position##
	if(savefiles == TRUE){
		write.table(NCBIList,paste(destfile,"full_feature_table.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	  write.table(FeatureList,paste(destfile,"subsetted_feature_table.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	  remove(ListSubAFC,NCBIList)
	}
	if(savefiles == FALSE){
		unlink(dest)
		remove(ListSubAFC,NCBIList)
	}
	cat("Finished processing. The list will now be returned to the user.","\n")
	return(FeatureList)
}
