## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SUFA)
genedata=download_genedata()
# data(genedata)

## ----set_seed-----------------------------------------------------------------
set.seed(1)

## ----extraction---------------------------------------------------------------
bulk.data= genedata$bulk
array.data= genedata$array1
array.data2= genedata$array2

## ----extraction_types---------------------------------------------------------
Ident.bulk= genedata$bulk.types
Ident.array= genedata$array1.types
Ident.array2= genedata$array2.types

## ----install_genefilter-------------------------------------------------------
if (!require("genefilter", quietly = TRUE)){
  BiocManager::install("genefilter")
}
library(Biobase);library(genefilter)

## ----hvg_bulk-----------------------------------------------------------------
dat=ExpressionSet(assayData = bulk.data)
filter.dat=varFilter(dat, var.cutoff = 0.95)
hvgs_in_bulk=rownames(filter.dat)
rm(filter.dat,dat)

## ----hvg_array1---------------------------------------------------------------
dat=ExpressionSet(assayData = array.data)
filter.dat=varFilter(dat, var.cutoff = 0.95)
hvgs_in_array=rownames(filter.dat)
rm(filter.dat,dat)

## ----hvg_array2---------------------------------------------------------------
dat=ExpressionSet(assayData = array.data2)
filter.dat=varFilter(dat, var.cutoff = 0.95)
hvgs_in_array2=rownames(filter.dat)
rm(filter.dat,dat)

## -----------------------------------------------------------------------------
genes.use=Reduce(intersect,list(hvgs_in_array,hvgs_in_array2,hvgs_in_bulk))
length(genes.use)

## ----gene_heatmap, fig.dim = c(7, 5), fig.cap="Heatmaps of two microarray and a bulk RNAseq dataset on a common set of genes. The X and Y axes represent genes and cells across studies, respectively."----
Y=list(t(array.data[genes.use,]), t(array.data2[genes.use,]),t(bulk.data[genes.use,] ))
Ident.list=list(Ident.array,Ident.array2,Ident.bulk)

Y_cent=mapply(function(Idents,y){
  cell.types=unique(Idents)
  y_cent=lapply(cell.types, function(ident,Identt,Y){
    scale(Y[which(Identt==ident),],scale=F)
  } ,Y=y,Identt=Idents)
  Reduce(rbind,y_cent)
}, Ident.list, Y )

Y_ordered=mapply( function(y, name){
  d=dist(y)
  hc2=cluster::agnes(d)
  library(reshape2)
  dimnames(y)=NULL
  melted_cormat <- melt(t((y[hc2$order,]/max(abs(y)))))
  melted_cormat$Dataset=as.factor(name)
  melted_cormat
},Y_cent,c("Microarray 1", "Microarray 2", "Bulk RNAseq"),SIMPLIFY = F)

Y_ordered=Reduce(rbind,Y_ordered)
library(ggplot2)
ggplot(data = Y_ordered, aes(x=Var1, y=Var2, fill=value))+  geom_tile(aes(fill = value)) +
  theme (panel.grid.major = element_blank(), panel.border = element_blank(),
         panel.background = element_blank(), legend.position = "none",text = element_text(size=15))+ 
  scale_fill_viridis_c(option="B",alpha=1) +  facet_grid(rows=vars(Dataset),scales = "free")+
  xlab("Genes")+ylab("Cells")

## ----echo=FALSE---------------------------------------------------------------
rm(Y_ordered)

## ----scaling------------------------------------------------------------------
sd.median=median( sapply(Y_cent, function(y) apply(y,2,sd)))
Y_cent_scaled=lapply(Y_cent, "/",   sd.median)

## ----fit_SUFA, message=FALSE, results=FALSE, warning=FALSE, comment=FALSE-----
res.sufa<-fit_SUFA(Y_cent_scaled,qmax=25,nthreads = 6,nrun = 7.5e3,
                                     nleapfrog = 4, leapmax = 9, del_range = c(0,.01))

## ----set_burnin---------------------------------------------------------------
burn=500

## ----wbic---------------------------------------------------------------------
WBIC(res.sufa ,Y_cent_scaled, model="SUFA",burn=5e2,ncores=1)

## ----get_cormat---------------------------------------------------------------
phis=apply(res.sufa$Lambda[,,-(1:burn)],3,identity, simplify = F)
diags=apply(res.sufa$residuals[-(1:burn),],1,identity, simplify = F)

cormats=abind::abind( mapply(function(phi, ps) cov2cor( tcrossprod(phi)+diag(ps)) 
                             ,phis,diags, SIMPLIFY = F),along = 3 )
cormat_mean=apply(cormats, c(1, 2), mean)

alpha=.05
low = apply(cormats, c(1, 2), quantile, alpha/2)
high = apply(cormats, c(1, 2), quantile, 1 - alpha/2)
cormat_mean[(low < 0) & (high > 0)]=0

## ----plot_cormat ,echo=FALSE,fig.dim= c(5.5,4.5), fig.cap ="Shared correlation structre of the genes"----
C2=cormat_mean; 
d <- as.dist(1-abs(C2)) 
# Hierarchical clustering using Complete Linkage
hc2=cluster::agnes(d)
C2[cormat_mean==0]=NA
heatplot(C2[hc2$order,hc2$order])

rm(C2)

## ----subgenes-----------------------------------------------------------------
Z2=apply(cormats, c(1, 2), mean)
rownames(Z2)=colnames(Z2)=genes.use
Z2[abs(Z2)<.25]=0
indices.connected.genes=which(rowSums( Z2!=0)-1 >10 )
subgenes=genes.use[indices.connected.genes]

## ----circos_plot, echo=FALSE, message=FALSE, warning=FALSE, comment=FALSE, dpi=200,out.width="120%", fig.cap ="Circos plot of the dependency structure between genes, blue (red) connection implies a positive (negative), their opacities being proportional to the corresponding association strengths. Associations are only plotted for absolute correlation >= 0.25."----
library(circlize)
C2=Z2[subgenes,subgenes]; 
d <- as.dist(1-abs(C2)) 
# Hierarchical clustering using Complete Linkage
hc2=cluster::agnes(d)
indices.connected.genes.ordered=hc2$order
subgenes.connected.ordered=subgenes[indices.connected.genes.ordered]

est_partcor_plot_subgenes_connected_ordered=C2[subgenes.connected.ordered,subgenes.connected.ordered] 
col_red_blue = colorRamp2(c(-.5,0,.5), c("red","white", "blue"), transparency = 0.00)
diag(est_partcor_plot_subgenes_connected_ordered)=.0001
par(cex=.2)
chordDiagramFromMatrix(est_partcor_plot_subgenes_connected_ordered, order=c(subgenes.connected.ordered), symmetric=TRUE, transparency = 0.25, grid.col="blue", col=col_red_blue, annotationTrack=c("grid"), keep.diagonal=T, scale=T, reduce=-100, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(est_partcor_plot_subgenes_connected_ordered))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

## ----loadings, message=FALSE, warning=FALSE, comment=FALSE,out.width="120%", fig.height=7, dpi=300, fig.cap ="Plots the top six columns with highest sum-of-squares values from the shared loading matrix, and repeat for the top three columns of each study-specific loading matrix."----
library(dplyr);library(ggplot2);library(dplyr);library(forcats);library(viridis)
# studies=c("Microarray 1", "Microarray 2","BulkRNASeq")
studies=c("Array 1", "Array 2","Bulk")
loadings=lam.est.all(res.sufa,burn = burn)

loadings_shared=data.frame(Loads= loadings$Shared[,1:6],Genes=genes.use, Loading="Shared");
loadings_study=foreach(s =1:length(loadings$Study_specific),.combine = rbind,.packages = c("SUFA","abind")) %do%{
  data.frame(Loads= loadings$Study_specific[[s]][,1:3] ,Genes=genes.use, Loading=studies[s])
}
loadings.m=rbind(melt(loadings_shared,id.vars = c("Loading","Genes")),
melt(loadings_study,id.vars = c("Loading","Genes")) )
loadings.m=subset(loadings.m, subset=Genes%in% (subgenes.connected.ordered<-subgenes[hc2$order]) )
loadings.m%>%
  mutate(Genes = fct_relevel(Genes,
                             subgenes[hc2$order] ),
         Loading= fct_relevel(Loading,c("Shared","Bulk","Array 1","Aarray 2") ) ) %>%
  ggplot( aes(Genes, (value), fill=value)) + 
  facet_wrap(~Loading+ variable, nrow=1,scales="free_x") + #place the factors in separate facets
  geom_bar(stat="identity") + ylab("Loading Strength") +coord_flip() + 
  theme_bw(base_size=10)+
  theme(legend.position = "bottom",axis.text.y = element_text(hjust = 1,size=2),
        axis.text.x = element_text(size=2, angle = 90), strip.text.x = element_text(size = 2),
        legend.key.height=unit(.01,"in"),legend.key.width=unit(.3,"in"), 
        legend.title = element_text(size = 2.5), legend.text = element_text(size = 2.5), 
        axis.title.x = element_text(size = 3), axis.title.y = element_text(size = 3))+
  scale_fill_viridis("Loading", option = "A") ->pplot
pplot

