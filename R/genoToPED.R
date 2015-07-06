genoToPED=function(saveAT,FamilyID=NULL,indID=NULL,PaternalID=NULL,MaternalID=NULL,Sex=NULL,pheno=NULL,geno,comment=NULL,pheno.Orimissing=NA,geno.Orimissing=NA,geno.Oricoding=c(0,1,2),geno.Pedcoding=c("AA","AB","BB"),chr=NULL,rs=NULL,Morgan=NULL,bp=NULL){
	#name for the ped and map file
	n=nrow(geno)
	nmar=ncol(geno)
	if(is.null(PaternalID)){PaternalID=rep(0,n)}
	if(is.null(MaternalID)){MaternalID=rep(0,n)}
	if(is.null(FamilyID)){FamilyID=rep(1,n)}
	if(is.null(indID)){indID=c(1:n)}
	if(is.null(Sex)){Sex=rep(3,n)}
	if(is.null(pheno)){
		pheno=rep(-9,n)
		}else{
	if(length(pheno)!=n) stop("number of phenotypes must equal to the number of rows in geno")	
	if(is.na(pheno.Orimissing)){
		if(any(is.na(pheno))){
			pheno[which(is.na(pheno))]=-9
			}
	}else{
	if(any(pheno==pheno.Orimissing)){pheno[which(pheno==pheno.Orimissing)]=-9}
	}
	}
	if(is.null(chr)){chr=rep(1,nmar)}else{if(length(chr)!=nmar){stop("chr must be a vector of size nmar")}}
	if(is.null(rs)){rs=paste(c(1:nmar),sep="")}else{if(length(rs)!=nmar)stop("rs must be a vector of size nmar")}
	if(!is.null(Morgan)){if(length(Morgan)!=nmar)stop("Morgan must be a vector of size nmar")}
	if(is.null(bp)){bp=c(1:nmar)}else{if(length(bp)!=nmar)stop("bp must be a vector of size nmar")}
	map=cbind(chr,rs,Morgan,bp)
	write.table(map,file=paste(saveAT,".map",sep=""),row.names=F,col.names=F,quote=F)
    saveAT.ped=paste(saveAT,".ped",sep="")
    
    if(is.na(geno.Orimissing)){
    	if(any(is.na(geno))){
    		geno[which(is.na(geno))]="00"
    	}
    }else{
    
    if(any(geno==geno.Orimissing)){
    	geno[which(geno==geno.Orimissing)]="00"
    	}
    }
	cat("#",comment,"\n",file=saveAT.ped)
	for(i in 1:nrow(geno)){
		genoi=geno[i,]
if(is.numeric(geno.Oricoding)){
	if(length(geno.Oricoding)!=length(geno.Pedcoding)){
		stop("geno.Oricoding must have the same number of elements as geno.Pedcoding")
		}
		for(j in 1:length(geno.Oricoding)){
	genoi[which(genoi==geno.Oricoding[j])]=geno.Pedcoding[j];
	}
	} 
	cat(FamilyID[i],indID[i],PaternalID[i],MaternalID[i],Sex[i],pheno[i],unlist(strsplit(genoi,split="")),"\n",sep=" ",append=T,file=saveAT.ped)
		
	}
	
}

getG=function(Pedgeno=NULL,geno,prune.r2=0.1,plinkpath="~/bin/plink",removePlinkFiles=T){
	#use plink to prune the genotypes
   if(prune.r2<1){
   if(is.null(Pedgeno)){
   saveAt=paste("prune")
   genoToPED(saveAT=saveAt,geno=wheat.X,geno.Oricoding=c(0,1),geno.Pedcoding=c("AA","BB"))
   }else{
   	saveAt=Pedgeno
   }
   system(paste(plinkpath, "--noweb --file", saveAt, "--indep-pairwise 100 10",prune.r2,sep=" "),show.output.on.console=F)
   prune.in=scan("plink.prune.in")
   geno=geno[,prune.in]
    if(removePlinkFiles==T){
      file.remove("plink.nosex")
      file.remove("plink.log")
   }
  
  }

  geno=scale(geno,center=T,scale=F)
  G=tcrossprod(geno)
  G=G/mean(diag(G))

  return(G)
}

