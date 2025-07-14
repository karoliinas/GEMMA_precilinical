setwd("/Users/ks423044/Desktop/GEMMA/mouse_model/Mouse_codes")
library(mixOmics)
library(BiocParallel)

# Read normalized and filtered data matrices - samples in rows and variables in columns
# I have named the mice using strain and the mouse id 
# e.g. C57BL/6J_34
# All the matrices need to have matching samples

# The two strains of mice were analyzed separately
#
#
#
# Metadata
group_data=read.table("Testdata/group_data.txt")
#
# C57BL/6J
qpcr=read.table("Testdata/Qpcr_black.txt")
tcell=read.table("Testdata/Tcell_black.txt")
metabolites=read.table("Testdata/Metabolites_black.txt")
microbiome_c=read.table("Testdata/Microbiome_genus_black.txt")
scfa=read.table("Testdata/Scfa_black.txt")
behaviour=read.table("Testdata/Behaviour_black.txt")

# vector with the hFMT status
group=group_data[rownames(qpcr),"Group"]

# BALB/C
#qpcr=read.table("Qpcr_black.txt")
#tcell=read.table("Tcell_black.txt")
#metabolites=read.table("Metabolites_black.txt")
#microbiome_c=read.table("Microbiome_genus_black.txt")
#scfa=read.table("Scfa_black.txt")
#behaviour=read.table("Behaviour_black.txt")

# vector with the hFMT status
#group=group_data[group_data$Strain=="BALB/c",]

# Select data blocks to integrate here

# Microbiome Caecum
xalb=list(qpcr,
          microbiome_c, 
          scfa,
          behaviour,
          metabolites)
names(xalb)=c("Qpcr", "MicrobiomeC","Scfa","Behaviour","Metabolites")


lapply(xalb, dim) # Check that the number of samples match

# Run PCA on the whole data
# Select number of components
pca_list=list()
for(i in 1:length(xalb)){
  p=pca(xalb[[i]], ncomp=4, scale=T)
  pca_list[[i]]=p
  
}
# Plot the variances by components
names(pca_list)=names(xalb)
par(mfrow=c(2,3))
for(i in 1:length(pca_list)){
  plot(pca_list[[i]], main=names(pca_list)[i])
}
# Run PCA on optimal number of components
pca_list=list()
for(i in 1:length(xalb)){
  p=pca(xalb[[i]], ncomp=2, scale=T) # change ncomp based on the variances
  pca_list[[i]]=p
  
}
# Plot the PCA results
par(mfrow=c(2,3))
for(i in 1:length(pca_list)){
  plotIndiv(pca_list[[i]], group=group,
            col = c('goldenrod','dodgerblue'),
            title=names(xalb)[[i]],
            size.title=rel(1),
            size.subtitle=rel(1),
            legend=F,
            ind.names=F,
            style="graphics",
            ellipse=T)
}


# Run basic sPLS on the pairs of datatypes to evaluate pairwise correlations
list.keepX=c(5,5)
list.keepY=c(5,5)

pls_list=list()
pls_names=c()
idx=combn(5,2)
for(i in 1:ncol(idx)){
  var1=names(xalb)[idx[1,i]]
  var2=names(xalb)[idx[2,i]]
  print(paste0(var1," vs ", var2))
  plsname=paste0("pls",var1,"_",var2)
  pls=pls(xalb[[var1]], xalb[[var2]])#, keepX=list.keepX, keepY=list.keepY)
  pls_list[[i]]=pls
  pls_names=append(pls_names, plsname)
}
names(pls_list)=pls_names

par(mfrow=c(2,2))
for(i in 1:length(pls_list)){
  var1=names(xalb)[idx[1,i]]
  var2=names(xalb)[idx[2,i]]
  plotVar(pls_list[[i]],title=names(pls_list)[[i]],
          legend=c(var1,var2), cutoff = 0.5, 
          var.names=F, style='graphics',
          pch = c(16, 17), cex = c(2,2), 
          col = c('thistle', 'peachpuff'))
}

# Plot the results
setwd("output_dir")
pdf("PLScorrelation_circles_microbiomeC.pdf", width=10, height=25)
par(mfrow=c(5,2))
for(i in 1:length(pls_list)){
  var1=names(xalb)[idx[1,i]]
  var2=names(xalb)[idx[2,i]]
  plotVar(pls_list[[i]],title=names(pls_list)[[i]],
          legend=c(var1,var2), cutoff = 0.5, 
          var.names=T, style='graphics',
          cex = c(0.7,0.7), 
          col = c('slategray', 'firebrick'))
}
dev.off()


# design matrix filled with 0.1s
design = matrix(0.1, ncol = length(xalb), nrow = length(xalb), 
                dimnames = list(names(xalb), names(xalb)))
diag(design) = 0 # set diagonal to 0s
yalb=as.character(group)
design

# form basic DIABLO model
basic.diablo.model = block.plsda(X = xalb, Y = yalb, ncomp = 4, design = design) 
plotDiablo(basic.diablo.model)

# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, ncomp=4, validation = 'Mfold', folds=3, nrepeat=100) 
# for small sample size use loo
#perf.diablo = perf(basic.diablo.model, ncomp=4, validation = 'loo') 

par(mfrow=c(1,1))
plot(perf.diablo) # plot output of tuning

# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

# Set the number of variables to test, change to suit each datatype
# Microbiome F
#test.keepX = list (Qpcr = c(1:5), 
#                   MicrobiomeF = c(1:4, seq(5, 12, 2)),
#                   Scfa = c(1:3),
#                   Behaviou. = c(1:4),
#                   Metabolites = c(1:4, seq(5, 12, 2)))

# Microbiome C
test.keepX = list (Qpcr = c(1:5), 
                   MicrobiomeC = c(1:4, seq(5, 12, 2)),
                   Scfa = c(1:3),
                   Behaviour = c(1:4),
                   Metabolites = c(1:4, seq(5, 12, 2)))

p=MulticoreParam(workers=parallel::detectCores()-2)
tuned=tune.block.splsda(X=xalb,Y=yalb, ncomp=2, dist="max.dist",test.keepX=test.keepX, validation="Mfold", folds=3, nrepeat=100,
                        design=design, BPPARAM = p) #Change ncomp and dist based on the results from perf.diablo, use validation=loo for small sample size


#Number of variables to keep
keepx=tuned$choice.keepX

#Default mode = PLS regression
final.diablo.model = block.splsda(X = xalb, Y = yalb, ncomp = 2, 
                                  design = design, keepX=keepx)


vars_selected=list()
# check the features selected
for(i in 1:length(xalb)){
  block=names(xalb)[i]
  vars=selectVar(final.diablo.model, block = block, comp = 1) #change comp = 2 for the second component
  vars_selected[[i]]=vars
}
names(vars_selected)=names(xalb)
#print the variables
for(i in 1:length(vars_selected)){
  print(paste0(names(vars_selected)[i],":"))
  print(vars_selected[[i]][[1]][[1]])
  print("")
}

# Plotting & saving results

group_col=c("goldenrod", "dodgerblue")
names(group_col)=c("ASD+GI","Siblings")

# PCA type of plot (1st component from each datatype against one another
# with correlation value between the datatypes)
par(mfrow=c(1,1))
pdf("black_mC_diablo_report_2comp.pdf",width=9, height=8)
plotDiablo(final.diablo.model,col.per.group=group_col,ncomp=1)
dev.off()

# PCA type of plot 1st and 2nd PLS component of each datatype
pdf("black_mC_diablo_sample_plot_2comp.pdf",width=9, height=8)
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots',
          legend.title="Group",
          group=group,
          size.title=rel(1),
          size.subtitle=rel(1),
          style="lattice",
          ellipse=T,
          comp=c(1,2),
          col = c("goldenrod","dodgerblue"))
dev.off()

# Arrow plot
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')


# Correlations in a circos plot
pdf("black_mC_circos_plot1:2_2comp.pdf", width=8, height=8)
circosPlot(final.diablo.model, comp=1:2,cutoff=0.8, group=as.factor(yalb),
           size.variables=0.7, size.labels=1, line=T,
           color.blocks= c("thistle",'tomato',"aquamarine4","lightblue","peachpuff"))
dev.off()

# Loadings plot
pdf("black_mC_Loadings_comp1_2comp.pdf", width=8, height=8)
plotLoadings(final.diablo.model, contrib = 'max', method = 'median', comp=1, size.name=1, size.title = rel(1),
             size.subtitle = rel(1))
dev.off()

pdf("black_mC_Loadings_comp2_2comp.pdf", width=8, height=8)
plotLoadings(final.diablo.model, contrib = 'max', method = 'median', comp=2, size.name=1, size.title = rel(1),
             size.subtitle = rel(1))

dev.off()
# Variates names with correlations on unit circle
par(mfrow=c(1,1))
pdf("black_mC_variates_2comp.pdf", width=10, height=10)
plotVar(final.diablo.model, var.names = T, 
        style = 'graphics', legend = TRUE)
dev.off()

# Heatmap of the selected variables
# 1st component
pdf("black_mC_heatmap_comp1_2comp.pdf", width=10, height=10)
cimDiablo(final.diablo.model, comp=1, size.legend = 0.8, color.Y=c("goldenrod","dodgerblue"), 
          color.blocks=c("thistle",'tomato',"aquamarine4","peachpuff","lightblue"),transpose=T, legend.position = "topright")
dev.off()

# 2 components
pdf("black_mC_heatmap_comp12_2comp.pdf", width=10, height=10)
cimDiablo(final.diablo.model, comp=1:2, size.legend = 0.8, color.Y=c("goldenrod","dodgerblue"), 
          color.blocks=c("thistle",'tomato',"aquamarine4","peachpuff","lightblue"),transpose=T, legend.position = "topright")
dev.off()


# Save the final model
saveRDS(final.diablo.model, "black_mC_diablo")

