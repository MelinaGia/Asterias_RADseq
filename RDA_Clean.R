library(vegan)
library(psych)
#load in data, imputed in LEA (see LFMM script)
gen.imp <- read.table("RAD_All_NewRef_NoOut_unlinkedK3.lfmm_imputed.lfmm", sep = " ", header = FALSE) 
sum(is.na(gen.imp))
#load in environmental layers for each sample with population information
env<-read.csv("envlfmm_combined_shortname_withpops.csv", header=T)
env$samples <- as.character(env$samples) # Make individual names characters
env$site <- as.factor(env$site) 
env$pop <- as.factor(env$pop)
env$HI <- as.factor(env$HI)
#investigate correlations
#this is derived from Brenda R. Forrester's tutorial: https://popgen.nescent.org/2018-03-27_RDA_GEA.html
identical(rownames(gen.imp), env[,1]) 
rownames(gen.imp)
rownames(env[,1])
pairs.panels(env, pch = ".", gap=0, scale = TRUE)
pred1 <- subset(env, select=-c(samples, site, pop, HI))
pairs.panels(pred1, pch = ".", gap=0, scale = TRUE)

?pairs.panels
#remove correlated: ppmean (PRO), phyto (PHY), Av temp (TAV)
predA <- subset(env, select=-c(PRO, PHY, OXY, TMAX, TMIN, TAV, CHL, samples, site, pop, HI))
str(predA)
pairs.panels(predA [,4:11], pch = ".", gap=0, scale = TRUE)

pred <- subset(env, select=-c(TMIN, TMAX, OXY, CLU, CHL, PHY, PRO, PHO, samples, site, pop, HI))
str(pred)
pairs.panels(pred, pch = ".", gap=0, scale = TRUE)
ast.rda <- rda(gen.imp ~ ., data=pred, scale=T)
ast.rda
RsquareAdj(ast.rda)
summary(eigenvals(ast.rda, model = "constrained"))
screeplot(ast.rda)
signif.full <- anova.cca(ast.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
signif.axis <- anova.cca(ast.rda, by="axis", parallel=getOption("mc.cores"), first = TRUE)
signif.axis
vif.cca(ast.rda)
plot(ast.rda, scaling=3)
plot(ast.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3
load.rda <- scores(ast.rda, choices=c(1:3), display="species") 
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 
levels(env$pop) <- c("R","H","F")

pops <- (env$pop)
pops
bg <- c("#1f78b4","purple","#ff7f00") # 6 nice colors for our pops
# axes 1 & 2
plot(ast.rda, type="n", scaling=3)
points(ast.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ast.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[pops]) # the wolves
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(pops), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ast.rda, type="n", scaling=3, choices=c(1,3))
points(ast.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(ast.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[pops], choices=c(1,3))
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(pops), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

levels(env$site) <- c("NL","ME","MA", "SC", "EU", "NS", "NH", "NY", "PEI", "NC")
site <- env$site
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c", "purple", "pink", "darkolivegreen", "gold") # 6 nice colors for our pops
# axes 1 & 2
plot(ast.rda, type="n", scaling=3)
points(ast.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ast.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[site]) # the wolves
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(site), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ast.rda, type="n", scaling=3, choices=c(1,3))
points(ast.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(ast.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[site], choices=c(1,3))
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(site), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

load.rda <- scores(ast.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda[,1],3) #148
cand2 <- outliers(load.rda[,2],3) #825
cand3 <- outliers(load.rda[,3],3) #137
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand #1128

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=10)  # 11 columns for 11 predictors
colnames(foo) <- c("CAL","CRV","IRO","NIT","RAD","PH","SAL","SIL","TAV","TRA")
foo
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
head(cand)
str(cand) #to see structure of dataframe- in below loop, start at first env variable (4) and end at last one.
length(cand$snp)
length(cand$snp[duplicated(cand$snp)])  # no duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,14] <- names(which.max(abs(bar[4:13]))) # gives the variable
  cand[i,15] <- max(abs(bar[4:13]))              # gives the correlation
}

colnames(cand)[14] <- "predictor"
colnames(cand)[15] <- "correlation"

table(cand$predictor)

sel <- cand$snp
env <- cand$predictor
env[env=="CAL"] <- '#1f78b4'
env[env=="TAV"] <- '#a6cee3'
env[env=="CRV"] <- '#6a3d9a'
env[env=="IRO"] <- '#e31a1c'
env[env=="NIT"] <- '#33a02c'
env[env=="RAD"] <- '#ffff33'
env[env=="PH"] <- '#fb9a99'
#env[env=="PHO"] <- '#b2df8a'
env[env=="SAL"] <- 'purple'
env[env=="SIL"] <- 'pink'
env[env=="TRA"] <- 'gold'

# color by predictor:
col.pred <- rownames(ast.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("V",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')

# axes 1 & 2
plot(ast.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(ast.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(ast.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("CAL","CRV","IRO","NIT","RAD","PH","SAL","SIL","TAV","TRA"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ast.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(ast.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(ast.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(ast.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("CAL","CRV","IRO","NIT","RAD","PH","SAL","SIL","TAV","TRA"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#Try to do this with each population separately 

pred.R<-pred[env.data$pop=="R",]
pred.F<-pred[env.data$pop=="F",]
pred.H<-pred[env.data$pop=="H",]


# call the genetic data file
head(gen.imp)
dim(gen.imp)
#separate out the different species, designated based on diem data
gen.AF <- gen.imp[c(3,4,6:15,17:26,35:44,49:51,53,54,56,57,59:61,63,64,76,88:98,111,112,115:123,125:127,130:136,138,139,155:163),]
gen.AR <- gen.imp[c(1,2,27:34,46,47,58,65:75,77,79:87,99,100,102,103,105:108,114,137,143,145:154),]
gen.H <- gen.imp[c(5,16,45,48,52,55,62,78,101,104,109,110,113,124,128,129,140:142,144),]
# write to file
write.table(gen.AF, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "AFonly_imputedk3_122322.lfmm")
write.table(gen.AR, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "ARonly_imputedk3_122322.lfmm")
write.table(gen.H, sep = " ", 
            col.names = FALSE, 
            row.names = FALSE,
            file = "Honly_imputedk3_122322.lfmm")
#Repeat everything from each data subset



