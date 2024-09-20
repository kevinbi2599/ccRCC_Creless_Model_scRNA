library(glmnet)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)
library(ComplexHeatmap)

library(Seurat)
library(useful)
library(harmony)

library(RColorBrewer)
library(Scillus)
library(ggrepel)
library(stringr)
library(dplyr)

set.seed(1234)

setwd("/broad/VanAllenLab/xchip/cga_home/kevinbi/StranskyLinehan_Mouse")

#Bring in functions from similarity.R.
#https://github.com/constantAmateur/scKidneyTumors/tree/master

getPopulationOffset = function(y){
  if(!is.factor(y))
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}


#' Do the OvR fit for every variable.  This just does a simple CV selection of regularisation amount.  Far from ideal, but should be good enough for the main conclusions.
multinomialFitCV = function(x,y,nParallel=1,...){
  fits = list()
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {
        tryCatch(
          cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,lambda=exp(seq(-10,-3,length.out=100)),...),
          error = function(e) {
            warning(sprintf("Could not fit model for variable %s",mark))
            return(NULL)
          })
      })
  }
  return(fits)
}

#' Load training data
loadTrainingData = function(dirs){
  mDat = NULL
  toc = NULL
  for(dir in dirs){
    metadata = readRDS(file.path(dir,'metadata.RDS'))
    ttoc = readRDS(file.path(dir,'tableOfCounts.RDS'))
    rownames(metadata$meta.data) = paste0(metadata$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(metadata$meta.data)))
    colnames(ttoc) = rownames(metadata$meta.data)
    if(is.null(mDat)){
      mDat = metadata$meta.data
    }else{
      mDat = rbind(mDat,metadata$meta.data)
    }
    if(is.null(toc)){
      toc = ttoc
    }else{
      toc = cbind(toc,ttoc)
    }
  }
  return(list(toc=toc,mDat=mDat))
}

###

#Read in normal kidney and tumor objects and prepare for training

seur_n = readRDS(file = "Manuscript_Objects/NonImmune_NormalKidney_AmbFilt_SeuratObject_Final.RDS")

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(seur_n, label = T)

seur_t = readRDS(file = "Manuscript_Objects/NonImmune_Tumor_AmbFilt_SeuratObject.RDS")

options(repr.plot.width=12, repr.plot.height=8)
DimPlot(seur_t, label = T)

#table(Idents(seur_n))
#table(Idents(seur_t))

n_omits = c("Immune_Contaminant","Unidentified","Unidentified_PT")
seur_n = subset(seur_n, idents = n_omits, invert = T)
seur_n$AnnotTrain = Idents(seur_n)

t_omits = c("Immune_Contaminant","Contaminant")
seur_t = subset(seur_t, idents = t_omits, invert = T)
table(Idents(seur_t))

seur_t$Annot_Granular = droplevels(seur_t$Annot_Granular)
table(seur_t$Annot_Granular)

meta = seur_t@meta.data
meta = mutate(meta, Annot_Test = ifelse(Annot_Granular %in% c("Tumor_Neat1","Tumor_Cpe","Tumor_Cycling_Mki67","Tumor_ISG_High_Ifit3"),paste(Annot_Granular, Sample, sep = "."),paste(Annot_Granular)))
table(meta$Annot_Test)

seur_t$Annot_Test = meta$Annot_Test
seur_t = seur_t
Idents(seur_t) = seur_t$Annot_Test
seur_t$AnnotTrain = Idents(seur_t)

seur_n_raw = seur_n@assays$RNA@counts
seur_t_raw = seur_t@assays$RNA@counts
shared_genes = intersect(row.names(seur_n_raw),row.names(seur_t_raw))

rtoc = cbind(seur_n_raw[shared_genes,], seur_t_raw[shared_genes,])
saveRDS(rtoc, file = "Manuscript_Objects/rtoc_NormalKidney_Tumor_Final.RDS")

#Read in results of feature selection DE
select = readRDS("Manuscript_Objects/NormalKidney_COI_FeatureSelection_Output_Final.RDS")
select = subset(select, p_val_adj < 0.001 & avg_log2FC > 0.5)
table(select$cluster)

model_genes = unique(as.character(select$gene))
length(model_genes)

#Original list of excluded features from Young et al, 2018
#hkGeneREGEX='^(Eif[0-9]|Rpl[0-9]|Rps[0-9]|Rpn1|Polr[0-9]|Snx[0-9]|Hsp[ab][0-9]|H1fx|H2af[vxyz]|Prka|Nduf[abcsv]|Psm[abcdefg][0-#9]|Uba[0-9]|Ube[0-9]|Usp[0-9]|Txn)'
#coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(rtoc),value=TRUE), #Poorly characterised
#                        grep('Malat1',rownames(rtoc),value=TRUE), #Contamination
#                        grep('^Hb',rownames(rtoc),value=TRUE), #Contamination
#                        grep('^mt-',rownames(rtoc),value=TRUE), #Mitochondria
#                        grep(hkGeneREGEX,rownames(rtoc),value=TRUE) #Housekeeping genes
#                        ))

#Custom list of excluded features
hkGeneREGEX='^(Fth1|Ftl1|Xist|Eif[0-9]|Rpl[0-9]|Rps[0-9]|Rpn1|Polr[0-9]|Snx[0-9]|Hsp[ab][0-9]|H1fx|H2af[vxyz]|Prka|Nduf[abcsv]|Psm[abcdefg][0-9]|Uba[0-9]|Ube[0-9]|Usp[0-9]|Txn)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',row.names(rtoc),value=TRUE), #Poorly characterised
                        grep('Malat1',row.names(rtoc),value=TRUE), #Contamination
                        grep('^Hb',row.names(rtoc),value=TRUE), #Contamination
                        grep('^mt-',row.names(rtoc),value=TRUE), #Mitochondria
                        grep(hkGeneREGEX,row.names(rtoc),value=TRUE) #Housekeeping genes
                        ))

###

#Train models for normal cell type identities

trainDat = list(toc = seur_n@assays$RNA@data, mDat = seur_n@meta.data)
trainDat$mDat$Trainer = trainDat$mDat$AnnotTrain
cells = c(rownames(trainDat$mDat))
corner(cells)
classes = c(trainDat$mDat$AnnotTrain)
corner(classes)

excludeGenes=coreExcludeGenes
includeGenes=c()

dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(Matrix::rowSums(dat>0)>10 & !(rownames(dat)%in%excludeGenes)),]
dat = dat[intersect(model_genes, row.names(dat)),]
dim(dat)
dat = t(dat)

fitNormalClusters = multinomialFitCV(dat, classes, nParallel = 2)

str(fitNormalClusters)

saveRDS(fitNormalClusters, file = "Manuscript_Objects/TrainedModels_NormalKidney_SelectedFeatures_Final.RDS")

