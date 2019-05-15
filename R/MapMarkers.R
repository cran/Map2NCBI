#############################################
### Code to create function "MapMarkers"  ###
### Hanna and Riley                       ###
#############################################

MapMarkers = function(features,markers,nAut,other=c("X"),savefiles=TRUE,destfile){
    chromosome = NULL
    if(missing(features)){
			stop("ERROR: Did not specify list of features to use for mapping.")
		}
		if(missing(markers)){
			stop("ERROR: Did not specify list of markers to be mapped.")
		}
		if(missing(nAut)){
			stop("ERROR: Did not specify the number of autosomes present in the marker file.")
		}
		if(savefiles == TRUE){
			if(missing(destfile)){
				stop("ERROR: No path was specified for the folder to save the output file.")
			}else{
			dest = paste(destfile,"MappedMarkers.txt",sep="")
			}
		}
    if(array(other)[1] == FALSE){
      chr = matrix(1:nAut,ncol=1)
      nchr = nAut
    } else{
      Aut = matrix(1:nAut,ncol=1)
      nchr = nAut
      other = matrix(other,ncol=1)
      Oth = matrix(0,dim(other)[1],1)
      for(i in 1:nrow(other)){
        Oth[i,1] = nchr+1
        nchr = nchr+1
      }
      chr = rbind(Aut,Oth)
      remove(Aut,Oth)
    }
		nCol = dim(array(c(colnames(markers),colnames(features),"Distance","Inside?")))
		MarkMap = matrix(0,1,nCol,byrow=TRUE,dimnames=list(c(1),array(c(colnames(markers),colnames(features),"Distance","Inside?"))))
		NotMapped = matrix(0,1,dim(markers)[2],dimnames = list(c(1),colnames(markers)))
		for(i in 1:nchr){
		  if(i > nAut){
		    j = i-nAut
		    k = other[j,1]
		  } else { k = i }
		  Chr_Features = subset(features, chromosome==k)
			Chr_Markers = subset(markers, chromosome==i)
			if(nrow(Chr_Markers) > 0 & nrow(Chr_Features) > 0){
				rownames(Chr_Markers) = 1:nrow(Chr_Markers)
				rownames(Chr_Features) = 1:nrow(Chr_Features)
				for(locus in 1:nrow(Chr_Markers)){
					MapIt = NULL
				  Distance = matrix(c(array(Chr_Features[,"start"] - Chr_Markers[locus,"position"]),array(Chr_Features[,"end"] - Chr_Markers[locus,"position"])),nrow=2,byrow=TRUE)
				  if(dim(array(Distance[1,which(Distance[1,]<0)]))[1]>=1){
  				  Inside = which(Distance[1,] < 0 & Distance[2,] >=0,arr.ind = TRUE)
  				  if(dim(array(Inside))[1]>1){
  				    Small = which(Distance[,array(Inside)]==min(abs(Distance[,array(Inside)])) | Distance[,array(Inside)]==-1*min(abs(Distance[,array(Inside)])),arr.ind = TRUE)
  				    Inside = array(Inside)[Small[2]]
  				  }
  				  if(dim(array(Inside))[1]==1){
	  			    MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[Inside,],cbind(matrix(0,1,1,dimnames = list(1,"Distance")),matrix("Yes,_Inside_Gene",1,1,dimnames = list(1,"Inside?")))))
	  			    MarkMap = rbind(MarkMap,MapIt)
	  			  }
				  }
  				if(dim(MarkMap[which(MarkMap[,1]==Chr_Markers[locus,1]),])[1]==0){
    				minDisLoc = which(Distance == min(abs(Distance)) | Distance == -1*min(abs(Distance)),arr.ind = TRUE)
    				minDis = Distance[minDisLoc[1],minDisLoc[2]]
    				if(minDis > 0){
    				  loc = "Before"
    				}else{ loc = "After" }
    				if(abs(minDis)<=2500 & abs(minDis)>0){
  					  MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[minDisLoc[2],],cbind(matrix(abs(minDis),1,1,dimnames = list(1,"Distance")),matrix(paste("Marker_is_<=_2500_bp_",loc,"_Feature",sep=""),1,1,dimnames = list(1,"Inside?")))))
  					}
  					if(abs(minDis)>2500 & abs(minDis)<=5000){
  					  MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[minDisLoc[2],],cbind(matrix(abs(minDis),1,1,dimnames = list(1,"Distance")),matrix(paste("Marker_is_>_2500_bp_<=5000_bp_",loc,"_Feature",sep=""),1,1,dimnames = list(1,"Inside?")))))
  					}
  					if(abs(minDis)>5000 & abs(minDis)<=25000){
  					  MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[minDisLoc[2],],cbind(matrix(abs(minDis),1,1,dimnames = list(1,"Distance")),matrix(paste("Marker_is_>_5000_bp_<=25000_bp_",loc,"_Feature",sep=""),1,1,dimnames = list(1,"Inside?")))))
  					}
  					if(abs(minDis)>25000 & abs(minDis)<=1000000){
  					  MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[minDisLoc[2],],cbind(matrix(abs(minDis),1,1,dimnames = list(1,"Distance")),matrix(paste("Nearest_feature_is_>_25,000_bp_",loc,"_Feature",sep=""),1,1,dimnames = list(1,"Inside?")))))
  					}
  					if(abs(minDis)>1000000){
  					  MapIt = cbind(Chr_Markers[locus,],cbind(Chr_Features[minDisLoc[2],],cbind(matrix(abs(minDis),1,1,dimnames = list(1,"Distance")),matrix(paste("Nearest_feature_is_>_1_Mb_",loc,"_Feature",sep=""),1,1,dimnames = list(1,"Inside?")))))
  					}
    				MarkMap = rbind(MarkMap,MapIt)
  				}
				}
			}else{
			  NotMapped = rbind(NotMapped,Chr_Markers)
			}
		}
		MarkMapF = MarkMap[2:nrow(MarkMap),]
		if(dim(NotMapped)[1]>1) {
      NotMapped = NotMapped[2:nrow(NotMapped),]
		}
    if(savefiles == TRUE){
			write.table(MarkMapF,dest,quote=FALSE,sep="\t",row.names=FALSE)
		  if(dim(NotMapped)[1]>0){
		    write.table(NotMapped,paste(destfile,"MarkersNotMapped.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
		  }
		}
		return(MarkMapF)
}
