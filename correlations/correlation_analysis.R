setwd("input_dir")
library(corrplot)

# read grouping information
group=read.table("group_data.txt")

# Read the data
intestine_balb=read.table("analysis/intestine_withNA.txt", header=T)
Tcell=read.table("analysis/Tcell_withNA.txt", header=T)
scfa=read.table("analysis/scfa_clr_norm.txt")
behaviour_balb=read.table("analysis/behaviour_NA.txt")
brain=read.table("analysis/brain_withNA.txt")
metabolites=read.table("analysis/Metabolites_imputed.txt")
microbiome=read.table("analysis/MicrobiomeC_phyl_CLR.tx")
alpha_diversity=read.table("analysis/MicrobiomeC_relative_abundance.txt")


# Function for running corrplot for asd vs sibling separately (plotted in upper and lower triangles)
# input: df=corr_df, strain="BALB/C" and n=number of columns to correlate (integer)
corplot=function(df, strain, n){
  asd=df[which(df$status=="ASD+GI" & df$strain==strain),]
  sib=df[which(df$status=="Siblings" & df$strain==strain),]
  print(dim(asd))
  print(dim(sib))
  cor.a = cor(asd[,1:n], method="spearman",use="pairwise.complete.obs")
  test.a = cor.mtest(asd[,1:n], conf.level = 0.95, method="spearman",exact = F,use="pairwise.complete.obs")
  cor.b = cor(sib[,1:n], method="spearman",use="pairwise.complete.obs")
  test.b = cor.mtest(sib[,1:n], conf.level = 0.95, method="spearman",exact = F,use="pairwise.complete.obs")
  cor.a[lower.tri(cor.a)] <- cor.b[lower.tri(cor.b)]
  test.a$p[lower.tri(test.a$p)] <- test.b$p[lower.tri(test.b$p)]
  corrplot(cor.a, p.mat=test.a$p, sig.level=0.05, diag=F)
  title(paste0("\n", "ASD+GI vs siblings ", strain, "\ntop: ASD+GI bottom: sibling"))
}

# Function to calculate P-value for correlations
# Using Fisher's Z-transformation
comparecorrlists_reija <- function(R1, R2, n1, n2){
  
  z1<-0.5*log((1+R1)/(1-R1))
  
  z2<-0.5*log((1+R2)/(1-R2))
  
  sezdiff<-sqrt(1.06/(n1-3) + 1.06/ (n2-3))
  
  ztesti<- (z1-z2)/sezdiff
  
  alfa<-2*(1 - pnorm(abs(ztesti)))
  
  return(alfa)
}

# Function to get the information about the correlation separately for sib and asd
cordat <- function(df, strain, variable1, variable2){
  asd_df=df[which(df$status=="ASD+GI" & df$strain==strain),]
  sib_df=df[which(df$status=="Siblings" & df$strain==strain),]
  cor.1 = cor(asd_df[,c(variable1, variable2)], method="spearman",use="pairwise.complete.obs")
  cor.2 = cor(sib_df[,c(variable1, variable2)], method="spearman",use="pairwise.complete.obs")
  test.asd = cor.mtest(asd_df[,c(variable1, variable2)], conf.level = 0.95, method="spearman",exact = F,use="pairwise.complete.obs")$p[1,2]
  test.sib = cor.mtest(sib_df[,c(variable1, variable2)], conf.level = 0.95, method="spearman",exact = F,use="pairwise.complete.obs")$p[1,2]
  asds=length(intersect(rownames(asd_df[!is.na(asd_df[,variable1]),]), #Number of ASD mice for each variable
                        rownames(asd_df[!is.na(asd_df[,variable2]),])))
  sibs=length(intersect(rownames(sib_df[!is.na(sib_df[,variable1]),]),
                        rownames(sib_df[!is.na(sib_df[,variable2]),]))) #Number of sib mice for each variable
  return(list("cor_asd"=cor.1[1,2], "cor_sib"=cor.2[1,2], "asdN"=asds, "sibN"=sibs, "SigASD"=test.asd, "sigSib"=test.sib))
}


# Generate a dataframe to run correlations for selected variables (edit selected
# columns based on variables of interest)
corr_df=merge(relative_abundance[,c("Shannon", "Chao1")], microbiome, by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=merge(corr_df,intestine[,c("Tph1_ileum","Il6_ileum","Il10_ileum","Tnf_ileum",
                                    "Tph1_colon","Il6_colon","Il10_colon","Tnf_colon")], by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=merge(corr_df,Tcell[,c("Th1/Th2","Th17/Treg","Th1/Treg")], by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=merge(corr_df,brain[,c("Slc6a4_PFC","Htr1a_PFC","Slc6a4_Hippocampus","Htr1a_Hippocampus")], by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=merge(corr_df,behaviour[,c("mZ.Locomotion","mZ.Social.behavior",
                                   "mZ.Anxiety","mZ.Stereotypies")], by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=merge(corr_df,metabolites[,c("X5.HIAA.5.HT",
                                     "Kyn.Trp",
                                     "KA.QA","X3.Indoxylsulfate",
                                     "Indole.3.lactic.acid",
                                     "Indole.3.acetic.acid")], by=0, all=T)
rownames(corr_df)=corr_df$Row.names
corr_df=corr_df[,2:ncol(corr_df)]
corr_df=mutate_all(corr_df, function(x) as.numeric(as.character(x)))

#Add the mouse type (strain) and status info
corr_df$strain=group_data$Strain
corr_df$status=group_data$status

setwd("output_dir")

#run corrplot for for ASD and Sibling, separately for the mouse strains
#balbc
pdf("correlation_balbc_AvsS.pdf",width=20, height=20)
corplot(corr_df, "BALB/c", 25) # 25 indicates the number of variables in the data (- status and strain columns)
dev.off()


# For loop for generating a data frame with the correlation comparisons
# on the corr_df

# initialize vectors to append the data
cor_asd=c()
cor_sib=c()
p_corr_asd=c()
p_corr_sib=c()
variable1=c()
variable2=c()
p_value=c()
n_asd=c()
n_sib=c()
idx=combn((ncol(corr_df)-2),2) # iterate through all the variables in the df
for(i in 1:ncol(idx)){
  var1=colnames(corr_df)[idx[1,i]]
  var2=colnames(corr_df)[idx[2,i]]
  cors=cordat(corr_df,"BALB/c", var1, var2) # change the strain based on mouse type
  p=comparecorrlists_reija(cors[[1]], cors[[2]], cors[[3]], cors[[4]]) #get the correlation comparison p-values
  variable1=append(variable1, var1) # append info to vectors
  variable2=append(variable2, var2) # append info to vectors
  p_value=append(p_value, p) # append info to vectors
  cor_asd=append(cor_asd,cors[[1]]) # append info to vectors
  cor_sib=append(cor_sib,cors[[2]]) # append info to vectors
  n_asd=append(n_asd,cors[[3]]) # append info to vectors
  n_sib=append(n_sib,cors[[4]]) # append info to vectors
  p_corr_asd=append(p_corr_asd, cors[[5]]) # append info to vectors
  p_corr_sib=append(p_corr_sib, cors[[6]]) # append info to vectors
  #print(paste0(var1, ", ", var2, " p: ", p))
}

#Adjusted p-values
p_adjust=p.adjust(p_value, method="fdr")

#Write the data after sorting
balbc_AvsS=data.frame("var1"=variable1, "var2"=variable2, "p_adj_comparison"=p_adjust, "p_val_comparison"=p_value, corAG=cor_asd, corAG_pvalue=p_corr_asd, corSAG=cor_sib, corSAG_pvalue=p_corr_sib, N_AG=n_asd, N_SAG=n_sib)
#Sort by p-values
balbc_AvsS=balbc_AvsS[order(balbc_AvsS$p_val_comparison, decreasing = F),]
write.table(balbc_AvsS, "CompareCorrelations_balbc_AvsS_caecum.tsv", sep="\t", col.names=T, row.names=F)


