### LFMM ###

library(LEA)
library(magrittr)
library(tidyr)
library(plyr)
library(dplyr)
library(RSpectra)
library(lfmm)
#load in environmental data, combined across samples and environmental layers
env.data<-read.csv("envlfmm_combined.csv", header=T)
#subset for each species, based on diem designations
AF_env <- env.data[c(3,4,6:15,17:26,35:44,49:51,53,54,56,57,59:61,63,64,76,88:98,111,112,115:123,125:127,130:136,138,139,155:163),]
AR_env <- env.data[c(1,2,27:34,46,47,58,65:75,77,79:87,99,100,102,103,105:108,114,137,143,145:154),]
H_env <- env.data[c(5,16,45,48,52,55,62,78,101,104,109,110,113,124,128,129,140:142,144),]

# Define LFMM run parameters
CPU = 4
repetitions = 10
#load unlinked SNP file 
lfmm1 = geno2lfmm("RAD_All_NewRef_NoOut_unlinked.geno")
#runsnmf and find best runfor each value of K
project.missing = snmf(lfmm1, K = 1:5, entropy = T, ploidy = 2, repetitions = 100, project= "new")
best1 = which.min(cross.entropy(project.missing, K = 1))
best2 = which.min(cross.entropy(project.missing, K = 2))
best3 = which.min(cross.entropy(project.missing, K = 3))
best4 = which.min(cross.entropy(project.missing, K = 4))
best5 = which.min(cross.entropy(project.missing, K = 5)) 
plot(project.missing)
#impute data for each value of K, and write it out to a file
imputeK1 = impute(project.missing, lfmm1, method = "mode", K = 1, run = best1)
file.rename("RAD_All_NewRef_NoOut_unlinked.lfmm_imputed.lfmm", "LEA_analyses/RAD_All_NewRef_NoOut_unlinkedK1.lfmm_imputed.lfmm")
imputeK2 = impute(project.missing, lfmm1, method = "mode", K = 2, run = best2)
file.rename("RAD_All_NewRef_NoOut_unlinked.lfmm_imputed.lfmm", "LEA_analyses/RAD_All_NewRef_NoOut_unlinkedK2.lfmm_imputed.lfmm ")
imputeK3 = impute(project.missing, lfmm1, method = "mode", K = 3, run = best3)
file.rename("RAD_All_NewRef_NoOut_unlinked.lfmm_imputed.lfmm", "LEA_analyses/RAD_All_NewRef_NoOut_unlinkedK3.lfmm_imputed.lfmm")
imputeK4 = impute(project.missing, lfmm1, method = "mode", K = 4, run = best4)
file.rename("RAD_All_NewRef_NoOut_unlinked.lfmm_imputed.lfmm", "LEA_analyses/RAD_All_NewRef_NoOut_unlinkedK4.lfmm_imputed.lfmm")
imputeK5 = impute(project.missing, lfmm1, method = "mode", K = 5, run = best5)
file.rename("RAD_All_NewRef_NoOut_unlinked.lfmm_imputed.lfmm", "LEA_analyses/RAD_All_NewRef_NoOut_unlinkedK5.lfmm_imputed.lfmm")

##########

# Define LFMM run parameters
# Use LFMM imputed file, only with the relevant K value, (no 9's)
gendata.imp <- read.table("RAD_All_NewRef_NoOut_unlinkedK3.lfmm_imputed.lfmm", sep = " ", header = FALSE) 
dim(gendata.imp)
#separate out the different species, based on diem results
AF <- gendata.imp[c(3,4,6:15,17:26,35:44,49:51,53,54,56,57,59:61,63,64,76,88:98,111,112,115:123,125:127,130:136,138,139,155:163),]
AR <- gendata.imp[c(1,2,27:34,46,47,58,65:75,77,79:87,99,100,102,103,105:108,114,137,143,145:154),]
H <- gendata.imp[c(5,16,45,48,52,55,62,78,101,104,109,110,113,124,128,129,140:142,144),]

# remove all constant SNPs. Goes from 25148 to 19809 SNPs
AF2 <- AF[,sapply(AF, function(v) var(v, na.rm=TRUE)!=0)]
dim(AF2)
AR2 <- AR[,sapply(AR, function(v) var(v, na.rm=TRUE)!=0)]
dim(AR2)
H2 <- H[,sapply(H, function(v) var(v, na.rm=TRUE)!=0)]
All2<-gendata.imp[,sapply(gendata.imp, function(v) var(v, na.rm=TRUE)!=0)]

# write to file
write.table(AF2, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "AFonly_imputedk3_122922.lfmm")
write.table(AR2, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "ARonly_imputedk3_122922.lfmm")
write.table(H2, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "Honly_imputedk3_122922.lfmm")
#load in the previously partitioned data, one at a time- run for loop between so  you don't need to change the variable
gendata.imp <- read.table("AFonly_imputedk3_122922.lfmm", sep = " ", header = FALSE) 
gendata.imp <- read.table("ARonly_imputedk3_122922.lfmm", sep = " ", header = FALSE) 
gendata.imp <- read.table("Honly_imputedk3_122922.lfmm", sep = " ", header = FALSE) 
#or just use the variables
gendata.imp <- H2
env.data <- H_env
#set lambda value
lambda = 1e-05
#define column names for the for-loop below
envdatanames = names(env.data[,2:19])

#define q-value, after inspecting several (0.1, 0.001. 0.0001)
n<-(0.01)
#this for loop will run through each column in the environmental data, and search for candidate genes 
for(q in n){
  for(j in 1:length(envdatanames)){
    envdata = env.data[envdatanames[j]] #function to run all of the environmental datasets, so you don't have to start a separate run for each of them
    
    project2 = lfmm_ridge(Y = gendata.imp, X = envdata, lambda = lambda, K = 3)
    pv <- lfmm_test(Y = gendata.imp, X = envdata, lfmm = project2, calibrate = "gif")
    pvalues <- pv$calibrated.pvalue 
    pvalues[is.na(pvalues)] <- 1
    
    name = paste("QQplot",envdatanames[j], "H", "q",q, sep="_")
    qqplot(rexp(length(pvalues), rate = log(10)),
           -log10(pvalues), xlab = "Expected quantile",
           pch = 19, cex = .4, main = name)
    abline(0,1)
    dev.copy(pdf, paste0(name, ".pdf"))
    dev.off()
    
    name = paste("Histogram", envdatanames[j], "H","q",q, sep="_")
    hist(pvalues, col = "red", main = name, breaks = 10)
    dev.copy(pdf, paste0(name, ".pdf"))
    dev.off()
    
    for (i in 1:ncol(envdata)){
      L = dim(gendata.imp)[2]
      w = which(sort(pvalues) < q*(1:L)/L) # Benjamini-Hochberg algorithm
      v = sort(pvalues)
      
      if(v[1] < (q*(1:L)/L)[1]){
        v.cutoff = length(v[v < q*(1:L)/L])
        pvalues.candidates = v[1:v.cutoff]
        candidates.list = order(pvalues)[w] # candidate SNPs sorted by adjusted p-value
        env.variable = rep(paste(i), v.cutoff)
        name = paste("Hcandidates.list.k3", envdatanames[j], "env", i, "q",q,sep="_")
        df = assign(name, cbind(candidates.list,pvalues.candidates,env.variable))
        write.csv(df, paste0(name, ".csv"), row.names = FALSE)
      }
      if(v[1] >= (q*(1:L)/L)[1]){
        v.cutoff = length(v[v < q*(1:L)/L])
        pvalues.candidates = NA
        candidates.list = NA
        env.variable = NA
        name = paste("Hcandidates.list.K3", envdatanames[j], "env", i,"q",q, sep="_")
        df = assign(name, cbind(candidates.list,pvalues.candidates,env.variable))
        write.csv(df, paste0(name, ".csv"), row.names = FALSE)
      }
    }
    
    lf <- list.files(pattern=glob2rx(paste0("*_", envdatanames[j], "_env*", ".csv")))
    df_list <- lapply(lf, read.csv, header = TRUE)
    total.candidate.SNPs <- bind_rows(df_list)
    write.csv(total.candidate.SNPs, file = paste0("H_", envdatanames[j], "q",q,"_K3_candidate_SNPs_perenv.csv"), row.names=FALSE)
    total.candidate.SNPs.unique <- unique(total.candidate.SNPs[,1])
    write.csv(total.candidate.SNPs.unique, file = paste0("UNIQUE_AR_", envdatanames[j],"q",q, "_K3_candidate_SNPs_perenv.csv"), row.names=FALSE)
    
    name = paste("Manhattan_plot1",envdatanames[j] , "H","q",q, sep="_")
    plot(-log10(pvalues), 
         pch = 19, 
         cex = .2, 
         xlab = "SNP", ylab = "-Log P", main = name,
         col = "grey")
    dev.copy(pdf, paste0(name, ".pdf"))
    dev.off()
    
    name = paste("Manhattan_plot2", envdatanames[j], "H","q",q, sep="_")
    points(as.numeric(df[,1]), 
           -log10(pvalues)[as.numeric(df[,1])], 
           type = "h", 
           col = "blue")
    
    dev.copy(pdf, paste0(name, ".pdf"))
    dev.off()
  }
  
}

####use predict function
q<-0.1
x<-env.data$BO22_phosphatemean_ss

mean(env.data$BO22_phosphatemean_ss)
y<-gendata.imp
mod = lfmm_ridge(Y = y, X = x, lambda = lambda, K = 2)
pv <- lfmm_test(Y = gendata.imp, X = envdata, lfmm = project2, calibrate = "gif")
pvalues <- pv$calibrated.pvalue 
pvalues[is.na(pvalues)] <- 1

L = dim(gendata.imp)[2]
w = which(sort(pvalues) < q*(1:L)/L) # Benjamini-Hochberg algorithm
v = sort(pvalues)

if(v[1] < (q*(1:L)/L)[1]){
  v.cutoff = length(v[v < q*(1:L)/L])
  pvalues.candidates = v[1:v.cutoff]
  candidates.list = order(pvalues)[w] # candidate SNPs sorted by adjusted p-value
  env.variable = rep(paste(i), v.cutoff)
  name = paste("candidates.list.K3", "dissox", "env", i, "q",q,sep="_")
  df = assign(name, cbind(candidates.list,pvalues.candidates,env.variable))
  write.csv(df, paste0(name, ".csv"), row.names = FALSE)
}

candidates <- as.numeric(df[,1]) #causal loci
b.values <- effect_size(Y = y, X = data.frame(x), lfmm.object = mod) 
x.pred <- scale(y[,candidates], scale = F)%*% matrix(b.values[candidates])

#population differentation tests
project.missing <- load.snmfProject("RAD_All_NewRef_NoOut_unlinked.snmfProject")
#try different K values and lambda values, then choose one. 
p3.2.5 = snmf.pvalues(project.missing, entropy = TRUE, ploidy = 2, K = 3, genomic.control = F, lambda = 2.5) 

?snmf.pvalues
pvalues3.2.5 = p3.2.5$pvalues
write.csv()
length(pvalues3.2.5)

par(mfrow = c(1,1)) 
hist(pvalues3.2.5, col = "orange") 
plot(-log10(pvalues3.2.5), pch = 19, col = "blue", cex = .5)

## FDR control: Benjamini-Hochberg at level q
L = length(pvalues3.2.5)
L
q = 1e-5
w = which(sort(pvalues3.2.5) < q * (1:L)/L)
candidates = order(pvalues3.2.5)[w] #list of candidates
length(candidates)
candidates
plot(-log10(pvalues3.2.5), main="Fst Outlier SNPs", xlab = "Locus", cex = .7, col = "grey")
points(candidates, -log10(pvalues2.2.5)[candidates], pch = 19, cex = .7, col = "red")
write.csv(candidates, file = "fstcandidates3.2.5_122922.csv")

#plot only LFMM Outliers
#Arubens
#Need to get from index SNP number to ACTUAL SNP number from full dataset 
Allcandidates <- read.csv("All_fstoverlap_010323.csv", header= FALSE)[,1]
Allcandidates
cand.names<-names(All2)[Allcandidates]

library(stringr)
cand.names<-as.numeric(str_remove(cand.names,"V"))
length(cand.names)
cand.names
cn1 <- unique(cand.names)
length(cn1)

L
L = length(pvalues3.2.5)
q = 1e-5
w = which(sort(pvalues3.2.5) < q * (1:L)/L)
candidates = order(pvalues3.2.5)[w] #list of candidates
candidates
length(candidates)
plot(-log10(pvalues3.2.5), main="LFMM and Fst outliers", xlab = "Locus", cex = .7, col = "grey")
points(cand.names, -log10(pvalues3.2.5)[cand.names], pch = 19, cex = .7, col = "red")
write.csv(candidates, file = "fstcandidates3.2.5.csv")
