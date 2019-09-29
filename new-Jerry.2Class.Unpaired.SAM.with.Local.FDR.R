rm(list=ls())

################# collapseIDs ######################

collapseIDs<-function(x,allids=row.names(x),method="mean"){

	allids<-as.vector(allids)
	ids<- levels(as.factor(allids))
	x.col<- matrix(nrow=length(ids), ncol=dim(x)[2])
	
	if(length(ids)==dim(x)[1]){ 
			dimnames(x)[[1]]<-allids
			return(x) 
	}
	
	for(i in 1:length(ids)){
		if(sum(allids==ids[i])>1){
			indices <- allids==ids[i] 
			if(method=="mean"){
				vals<-apply(x[indices,],2,mean)
			}
			if(method=="median"){
				vals<-apply(x[indices,],2,median)
			}
			if(method=="stdev"){   
				temp<- x[indices,]
				stdevs<- apply(temp,1,sd)
				vals<- temp[match(max(stdevs),stdevs),]
			}
			x.col[i,] <- vals
		}else{
			x.col[i,] <- x[allids==ids[i],]
		}
	}

	dimnames(x.col)<- list(ids,dimnames(x)[[2]])
	return(x.col)
	
}



################# readarray ######################

readarray<-function(dataFile,hr=1,impute=T,method="mean"){
	headerRows <- hr
	x<-read.delim(dataFile,sep="\t",header=T,fill=T)
	if(headerRows==1){
			classes<-NULL
			ids<-as.vector(t(x[,1]))
			xd<-x[,-1]
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)	
	}else{
			classes<-x[1:(headerRows-1),]
			dimnames(classes)[[1]]<-classes[,1]
			classes<-classes[,-1]
			xd<-x[(-1:-(headerRows-1)),]
			ids<-as.vector(t(xd[,1]))
			xd<-xd[,-1]
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)
	}
	features<- dim(xd)[1]
	samples<- dim(xd)[2]
	sampleNames<- names(xd)
	geneNames<-row.names(xd)
	xd<-apply(xd,2,as.numeric)
	row.names(xd)<-geneNames
	if(sum(apply(xd,2,is.na))>0 & impute){
		library(impute)
		allAnn<-dimnames(xd)
############################################################################################################################################################
		#data.imputed<-impute.knn(as.matrix(xd))############################################
############################################################################################################################################################
		data.imputed<-as.matrix(xd)

		xd<-data.imputed[1:features,]
		dimnames(xd)<-allAnn
	}
	classes[classes==""]<-NA
	return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}

################# readarray ######################

library(samr)

################################################################################################################################

#setwd("C:\Users\usary1\Documents\Work-Files\R scripts\from Chris Fan")
setwd("C:/Users/usary1/Documents/Work-Files/R scripts/from Chris Fan")

x<-readarray("KPB1B__NSG__NT_vs_CT_10Days.txt",hr=2)

for(i in 1:1)
{
	data=list(x=x$xd,y=x$classes[i,], geneid=as.character(rownames(x$xd)), genenames=as.character(rownames(x$xd)), logged2=T)          
	samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=300)

	delta= 0.5
	
	tiff(paste("SAM.plot.",rownames(x$classes[i,]),".tiff",sep=""))
	par(cex=2, font.lab=4, mar=c(4,4,1,1), mfrow=c(1,1))
	samr.plot(samr.obj,delta)
	legend("topleft",bty="n",legend=rownames(x$classes[i,]))
	dev.off()

	delta.table <- samr.compute.delta.table(samr.obj)
	siggenes.table<-samr.compute.siggenes.table(samr.obj, 0, data, delta.table,,all.genes=T,compute.localfdr=T)

	temp <- rbind(siggenes.table$genes.up,siggenes.table$genes.lo)
	write.table(temp, paste(rownames(x$classes[i,]),".samOut.V2.txt",sep=""), sep="\t",row.names=F)
	
}













