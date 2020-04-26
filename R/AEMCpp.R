## 11 April 2020 / Thinh:
# add a parameter "tau" to do cluster for heterogeneous samples
# rewrite AEM, EM, getSE in C++
## 27 Nov 2019 / Nghia:
# add standard error for the estimates
## 04 June 2019 / Nghia:
# add parameter: isoform.method=average/total to report the expression of the individual members of a paralog i) average (default) or ii) total from the paralog set
## 01 Apr 2019 / Wenjiang:
# add "merge.paralogs" parameter to turn on/off the paralog merging in XAEM. The default is off, which will generate the same set of isoforms between different projects. To turn it on, just add "merge.paralogs=TRUE"
# Example of command: Rscript buildCRP.R in=eqClass.txt isoform.out=X_matrix.RData merge.paralogs=TRUE

## Take the workdir and core arguments
workdir=NULL
core = 8 #default
merge.paralogs = TRUE ## default is to combine paralogs in the updated X to obtain the best performance
foutr="XAEM_paralog_expression.RData"
design.matrix="X_matrix.RData"
isoform.method="average" #  "average" or "total"
remove.ycount=TRUE
tau=0

saveSubset=FALSE #save singleton and paralogs 
noBiasResult=FALSE #export the results without bias correction
foutr_noBias="XAEM_noBiasCor.RData"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
  if (res[1]=="isoform.out") fout=as.character(res[2])
  if (res[1]=="paralog.out") foutr=as.character(res[2])
	if (res[1]=="merge.paralogs") merge.paralogs=as.logical(res[2])
  if (res[1]=="isoform.method") isoform.method=as.character(res[2])
  if (res[1]=="remove.ycount") remove.ycount=as.logical(res[2])
  if (res[1]=="tau") tau=as.numeric(res[2])
}

cat("\n AEM_update_X_beta.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n workdir: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
cat("\n isoform.out: ",fout)
cat("\n paralog.out: ",foutr)
cat("\n merge.paralogs: ",merge.paralogs)
cat("\n isoform.method: ",isoform.method)
cat("\n remove.ycount: ",remove.ycount)
cat("\n tau: ",tau)
cat("\n ----------------------------------------------------- ")

library(Rcpp)
library(RcppParallel)
library(parallel)

source("/path/to/R/Rsource.R")
source("/path/to/R/cluster.R")
sourceCpp("/path/to/R/Rsource.cpp")


options(stringsAsFactors=FALSE)
setwd(workdir)

#load input data
load(design.matrix)
load("Ycount.RData")

#set parallel
core = min(core, defaultNumThreads())
RcppParallel::setThreadOptions(numThreads = core)

##### start from here
###############################################################################
cat("\n...Estimation using AEM algorithm...\n")

X.y= Y

x_est = list()
y_est = list()
x_colnames = list()
for(i in 1:length(Y)){
  x.y =  X.y[[i]]
  xloc = which(colnames(x.y) != "sample1")
  y = x.y[,-xloc]
  x = matrix(x.y[,xloc], ncol=length(xloc))
  colnames(x) = colnames(x.y)[xloc]
  x_colnames[[i]] = colnames(x)
  x_est[[i]] = x
  y_est[[i]] = y
} 

AEMfun = function(i){
  res = AEMCpp(x_est[[i]], y_est[[i]], 5)
  colnames(res$X) = x_colnames[[i]]
  res
}

EST = mclapply(1:length(X.y), AEMfun, mc.cores=core)

###############################################################################
if(tau > 0){
  tau = 1 - tau
  cat("\n...Clustering for heterogeneous dataset...\n")
  res = exec_cluster(Y, EST, tau)
  original_Y = res[["original_Y"]]
  cluster_Y = res[["cluster_Y"]]
  original_EST = res[["original_EST"]]
  cluster_EST = mclapply(1:length(cluster_Y), function(i){
    y = cluster_Y[[i]]
    xloc = which(colnames(y) != "sample1")
    ys = y[,-xloc]
    xs = matrix(y[,xloc], ncol=length(xloc))
    res = AEMCpp(xs, ys)
    colnames(res$X) = colnames(y)[xloc]
    res
  }, mc.cores=core)

  EST = c(original_EST, cluster_EST)
  Y = c(original_Y, cluster_Y)
  y_est = lapply(Y, function(y){
    xloc = which(colnames(y) != "sample1")
    y[,-xloc]
  }) 
}

###############################################################################
X.y= Y

x.all = lapply(1:length(X.y), function(i){
  x = EST[[i]]$X
  x
})

# run add paralog step with X from EST result
if(merge.paralogs){
  cat("\n...Merging paralogs...\n")
  for(i in 1:length(X.y)){
    temp = try(ccrpfun(x.all[[i]]), silent = TRUE)
    if(class(temp) == "try-error"){
      temp = ccrpfun(x.all[[i]], clim=5)
    }
    x.all[[i]] = temp
  }
}

###############################################################################
beta.all = list()# use the new X to calculate new beta
se.all = list()# keep standard error

beta.all = lapply(1:length(X.y), function(i){
  res = parallelEMCpp(x.all[[i]], y_est[[i]])
  rownames(res) = colnames(x.all[[i]])
  t(res)
})
se.all = lapply(1:length(X.y), function(i){
  res = getSECpp(x.all[[i]], beta.all[[i]], y_est[[i]])
  colnames(res) = colnames(x.all[[i]])
  res
})

if (saveSubset) save(beta.all,file="Beta_final_paralog.Rdata")

###############################################################################
# process singletons

estfun = function(mat){
  CRP.y = mat$crpcount
  #sum(!(names(CRP.y)==names(CRP))) 
  ## cluster info
  npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
  table(npat)
  loc1 = which(npat==1)  # clusters with 1 pattern

  TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
   names(TC1)=names(CRP.y[loc1])
  ## single tx
  est.all =  c(TC1)
  return(est.all)
}

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
result_est=NULL
for(id in 1:length(flist)){ # call crpcount()
  load(flist[id])
  est1 = estfun(y)# estimation step
  result_est = cbind(result_est,est1)
}
colnames(result_est)=samplename1

if (saveSubset) save(result_est,file="Est_result_Singletons.RData")

###############################################################################
# expand beta matrix to have the same paralogs.
# split paralog to separeate isoforms for se matrix 
if(tau > 0){
  start = length(original_Y) + 1
  n_sample = nrow(original_EST[[1]]$BETA)

  #expand for beta matrix
  concat_beta = list()
  temp = NULL
  for(i in start:length(beta.all)){
      if(is.null(temp)){
        temp = beta.all[[i]]
        if(nrow(temp) == n_sample){
          concat_beta[[length(concat_beta) + 1]] = temp
          temp = NULL
        }
        next 
      }
      paralog_1 = colnames(temp)
      paralog_2 = colnames(beta.all[[i]])
      diff_1 = setdiff(paralog_2, paralog_1)
      diff_2 = setdiff(paralog_1, paralog_2)
      if(length(diff_1) > 0){
        for(j in 1:length(diff_1)){
          temp = cbind(temp, new=0)
        }
        colnames(temp) = c(paralog_1, diff_1)
      }
  
      new_beta = beta.all[[i]]
      if(length(diff_2) > 0){
        for(j in 1:length(diff_2)){
          new_beta = cbind(new_beta, new=0)
        }
        colnames(new_beta) = c(paralog_2, diff_2)
      }

      temp = as.matrix(rbind(as.data.frame(temp), as.data.frame(new_beta)))
      if(nrow(temp) == n_sample){
        concat_beta[[length(concat_beta) + 1]] = temp
        temp = NULL
      }
  }
  beta.all = c(beta.all[1:(start-1)], concat_beta)

  #split for se matrix
  temp = NULL
  concat_se = list()
  for(i in start:length(se.all)){
    se_isoforms = lapply(colnames(se.all[[i]]), function(s) unlist(strsplit(s, " ")))
    n_isoforms = sapply(se_isoforms, length)
    col_id = lapply(1:length(n_isoforms), function(k) rep(k, n_isoforms[k]))
    col_id = unlist(col_id)
    new_se = se.all[[i]][, col_id]
    colnames(new_se) = unlist(se_isoforms)

    temp = as.matrix(rbind(as.data.frame(temp), as.data.frame(new_se)))

    if(nrow(temp) == n_sample){
      concat_se[[length(concat_se) + 1]] = temp
      temp = NULL
    }
  }
  se.all = c(se.all[1:(start-1)], concat_se)
}

###############################################################################
seq1 = do.call(cbind,beta.all)
seq1 = t(seq1)
XAEM_count = rbind(result_est,seq1)

if(merge.paralogs){ #export standard error only if using merge.paralogs=TRUE
  se1=do.call(cbind,se.all);se1=t(se1)
  result_se=sqrt(result_est) #use sqrt(ycount) for singletons
  XAEM_se=rbind(result_se,se1)

  foutr_se=foutr
  foutr_se=gsub(".RData","",foutr_se)
  foutr_se=paste(foutr_se,".standard_error.RData",sep="")
  save(XAEM_se, file=foutr_se)
}
##### done with the estimation
###############################################################################

### collect the list of txnames from count data, more than 1 transcripts if the isoform is a paralog
txList = sapply(rownames(XAEM_count),function(x){return(unlist(strsplit(x," ")))})
txNum = sapply(txList,length)
### get median txlength for paralog
txlength_names = names(txlength)
txLen=mclapply(txList, function(x){
  pick = txlength_names %in% x
  return(median(txlength[pick]))
}, mc.cores=core)
txLen = unlist(txLen)
### compute TPM
#normalize count to length
isoform_lenNorm=apply(XAEM_count,2,function(x)return(x/txLen))
libsize_lenNorm=apply(isoform_lenNorm,2,sum)
XAEM_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
XAEM_tpm=t(XAEM_tpm)
#keep information of raw output of XAEM
save(XAEM_count, XAEM_tpm,file=foutr)

## expand isoforms can not separated from CRP
paralogID=which(txNum >1)
paralog_count=XAEM_count[paralogID,]
paralog_tpm=XAEM_tpm[paralogID,]
#get mapping ID
matchID=sapply(c(1:length(paralogID)), function(x) rep(paralogID[x],txNum[paralogID][x]))
matchID=unlist(matchID)
#get isoform names
matchNames=sapply(c(1:length(paralogID)), function(x) unlist(strsplit(names(paralogID[x])," ")))
matchNames=unlist(matchNames)
#update to data XAEM_count
expandDat=matrix(0,ncol = ncol(XAEM_count), nrow = sum(txNum[paralogID]))
expandDat=XAEM_count[matchID,] #if (isoform.method=="total"):the counts of isoform members are equal to the count of paralog
if (isoform.method=="average"){ #the counts of isoform members are equal to the average count of paralog
  expandDat_txnum=txNum[match(names(matchID),names(txNum))]
  expandDat2=apply(cbind(expandDat_txnum,expandDat),1,function(x) x[-1]/x[1])
  expandDat=t(expandDat2)
}
rownames(expandDat)=matchNames
isoform_count=XAEM_count
isoform_count=rbind(isoform_count,expandDat)
isoform_count=isoform_count[-paralogID,]
isoform_count=rowsum(isoform_count, group=rownames(isoform_count), na.rm=TRUE)
##recalculate TPM
txLen2=txlength[match(rownames(isoform_count),names(txlength))]
isoform_lenNorm=apply(isoform_count,2,function(x)return(x/txLen2))
libsize_lenNorm=apply(isoform_lenNorm,2,sum)
isoform_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
isoform_tpm=t(isoform_tpm)

#export to file
pick=names(txlength) %in% rownames(isoform_count)
pick=which(!pick)
if(length(pick)>0){
  newDat=matrix(0,ncol = ncol(isoform_count), nrow = length(pick))
  rownames(newDat)=names(txlength)[pick]
  colnames(newDat)=colnames(isoform_count)
  isoform_count=rbind(isoform_count,newDat)
  isoform_tpm=rbind(isoform_tpm,newDat)
}
isoform_count=isoform_count[names(txlength),]
isoform_tpm=isoform_tpm[names(txlength),]

save(isoform_count,isoform_tpm,file=fout)


###### No bias correction 
if (noBiasResult){
  cat("\n Get results with no bias correction")
  #paralogs
  X.y= Y
  beta.all = list()
  for(i in 1:length(X.y))
  {
   x2=NULL
   x.y =  X.y[[i]]
   #xloc = grep('N', colnames(x.y))
   xloc = which(colnames(x.y) != "sample1")
   y = x.y[,-xloc]
   x2 = x.y[,xloc,drop=FALSE] ## X is not updated using AEM algorithm
    beta1 = foreach(j=1:ncol(y)) %dopar% EM(x2,y[,j])
   beta2 = do.call(rbind,beta1)
   rownames(beta2) = NULL
  }
  #save(beta.all,file="Beta_final_paralog.RData")
  ## singletons
  estfun = function(mat){
    CRP.y = mat$crpcount
    sum(!(names(CRP.y)==names(CRP))) 
    ## cluster info
    npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
    table(npat)
    loc1 = which(npat==1)  # clusters with 1 pattern

    TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
     names(TC1)=names(CRP.y[loc1])
    ## single tx
    est.all =  c(TC1)
    return(est.all)
  }

  flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
  result_est=NULL
  for(id in 1:length(flist)){ # call crpcount()
    load(flist[id])
    est1 = estfun(y)# estimation step
    result_est = cbind(result_est,est1)
  }
  colnames(result_est)=samplename1
  seq1 = do.call(cbind,beta.all);seq1=t(seq1)
  XAEM_count = rbind(result_est,seq1)

  ### collect the list of txnames from count data, more than 1 transcripts if the isoform is a paralog
  txList = sapply(rownames(XAEM_count),function(x){return(unlist(strsplit(x," ")))})
  txNum=sapply(txList,length)
  ### get median txlength for paralog
  txLen=sapply(txList, function(x){
    pick=names(txlength) %in% x
    return(median(txlength[pick]))
  })
  ### compute TPM
  #normalize count to length
  isoform_lenNorm=apply(XAEM_count,2,function(x)return(x/txLen))
  libsize_lenNorm=apply(isoform_lenNorm,2,sum)
  XAEM_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
  XAEM_tpm=t(XAEM_tpm)
  #keep information of raw output of XAEM
  save(XAEM_count, XAEM_tpm,file=foutr_noBias)
} #done for no bias


### clean Ycount
# if (remove.ycount){
#   system("rm -rf Ycount")
#   system("rm -rf Ycount.RData")
# }

cat("\n...Done...\n")
