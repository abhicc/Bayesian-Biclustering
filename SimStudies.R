#==========================================================================================
#==========================================================================================
# This is additional code for simulation studies 6.1.3, 6.2, and 6.3
#==========================================================================================
#==========================================================================================

library(reshape2);library(tidyverse);library(ggpubr)

#==========================================================================================
# dataset generation
#==========================================================================================
I=100;J=200   # total number of rows and columns in dataset
R=5;C=5     # "true" number of row and column blocks

set.seed(1)
truthCs = rep(1:R, each=(I/R))
truthDs = rep(1:C, each=(J/C))
mus = runif(R*C,-3,3)
mus = matrix(c(mus),nrow=R,ncol=C,byrow=FALSE)
X = matrix(rnorm(I*J,mean=0,sd=5),nrow=I,ncol=J)

musmatrix = matrix(NA,nrow=I,ncol=J)
for(i in 1:max(truthCs)){
  for(j in 1:max(truthDs)){ 
    X[truthCs==i,truthDs==j] = X[truthCs==i,truthDs==j]+mus[i,j]
    musmatrix[truthCs==i,truthDs==j]<-mus[i,j]
  } 
}	

ids=matrix(c(1:(R*C)),nrow=R,ncol=C,byrow=TRUE);idsmatrix=matrix(nrow=I,ncol=J)
for(i in 1:max(truthCs)){for(j in 1:max(truthDs)){idsmatrix[truthCs==i,truthDs==j]<-ids[i,j]}}	
Es=matrix(nrow = I*J,ncol=2);Es[,1]=rep(1:I,J);Es[,2]=rep(1:J,each=I)
truthEs=rep(NA,I*J);for(i in 1:(I*J)){truthEs[i]=idsmatrix[Es[i,1],Es[i,2]]}

saveRDS(X,"SimStudies_X.rds")

#==========================================================================================
# Performing Bayesian biclustering (with real dataset SimStudies_X.rds) using either the first or 
# second version (see master codes) of the MCMC sampler. 
# Use 'tCs = truthCs', 'tDs = truthDs', and 'tEs=truthEs' in functions 'row_ri', 'col_ri', and 'element_ri'.
# Obtain (save) the results and plots.
# The following lines of code can be used to replicate the results and plots from the simulation studies.
#==========================================================================================
setwd(".........")   # set working directory
load("....RData")   # load saved results

#==========================================================================================
# Simulation study 6.1.3
#==========================================================================================

library(coda)
# Geweke's convergence diagnostic
gdsse = as.vector(geweke.diag(as.mcmc(sse_dat$SSE)))   # for total SSE values
gdsse

gdlp = as.vector(geweke.diag(as.mcmc(lp_dat$log_posterior)))   # for log posterior values
gdlp

# Use 'gelman.diag' function for obtaining Gelman and Rubin's convergence diagnostic for multiple chains.

#==========================================================================================
# Simulation study 6.2
#==========================================================================================

# perform k-means clustering
km.Cs = kmeans(train,5,nstart=20)$cluster
km.Ds = kmeans(t(train),5,nstart=20)$cluster
km.mus = matrix(NA,nrow=I,ncol=J)

for(i in 1:I){
  for(j in 1:J){
    km.mus[i,j] = mean(train[km.Cs==km.Cs[i],km.Ds==km.Ds[j]])		
  }
}

# plot true underlying means matrix
musmatrix = data.frame(musmatrix)
plotmus = expand.grid(Rows = rownames(musmatrix), Columns = colnames(musmatrix))
vec = as.vector(as.matrix(musmatrix))
plotmus$Z = vec

color = c(rainbow(14,start=0.45,end=0.6)[14:1],rainbow(16,start=0.05,end=0.2)[16:1])
g.mus = ggplot(plotmus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotmus$Rows)))) + 
  ggtitle("True underlying means") + 
  theme(plot.title = element_text(hjust = 0.5))
g.mus + rremove("xy.text") 

# plot means matrix from k-means implementation
km.mus = data.frame(km.mus)
plotkm.mus = expand.grid(Rows = rownames(km.mus), Columns = colnames(km.mus))
vec = as.vector(as.matrix(km.mus))
plotkm.mus$Z = vec

color = c(rainbow(14,start=0.45,end=0.6)[14:1],rainbow(16,start=0.05,end=0.2)[16:1])
g.km.mus = ggplot(plotkm.mus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotkm.mus$Rows)))) + 
  ggtitle("5-means clustering") + 
  theme(plot.title = element_text(hjust = 0.5))
g.km.mus + rremove("xy.text") 

# plot means matrix from Bayesian biclustering implementation
bbc.mus = matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  for(j in 1:J){
    # change criterion (see paper and master codes)
    bbc.mus[i,j] = mean(train[id_row_mat[,which.max(avgRowRI)]==id_row_mat[i,which.max(avgRowRI)],
                         id_col_mat[which.max(avgColRI),]==id_col_mat[which.max(avgColRI),j]])		
  }
}

bbc.mus = data.frame(bbc.mus)
plotbbc.mus = expand.grid(Rows = rownames(bbc.mus), Columns = colnames(bbc.mus))
vec = as.vector(as.matrix(bbc.mus))
plotbbc.mus$Z = vec

color = c(rainbow(14,start=0.45,end=0.6)[14:1],rainbow(16,start=0.05,end=0.2)[16:1])
g.bbc.mus = ggplot(plotbbc.mus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotbbc.mus$Rows)))) + 
  ggtitle("Bayesian biclustering") + 
  theme(plot.title = element_text(hjust = 0.5))
g.bbc.mus + rremove("xy.text") 

# plot means matrices on the same plot
ggarrange(g.mus+rremove("xy.text"), g.km.mus+rremove("xy.text"), g.bbc.mus+rremove("xy.text"),
          labels=c("(a)","(b)","(c)"),ncol=3,common.legend = TRUE,font.label = list(size=5))

#==========================================================================================
# Simulation study 6.3
#==========================================================================================

#sparse biclustering
library(sparseBC)
bicluster = sparseBC(train,5,5,0)
bc.mus = bicluster$mus

# plot means matrix from sparse biclustering implementation
bc.mus=data.frame(bc.mus)
plotbc.mus <- expand.grid(Rows = rownames(bc.mus), Columns = colnames(bc.mus))
vec = as.vector(as.matrix(bc.mus))
plotbc.mus$Z <- vec

color = c(rainbow(14,start=0.45,end=0.6)[14:1],rainbow(16,start=0.05,end=0.2)[16:1])
g.bc.mus=ggplot(plotbc.mus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotbc.mus$Rows)))) + 
  ggtitle("Sparse biclustering") + 
  theme(plot.title = element_text(hjust = 0.5))
g.bc.mus + rremove("xy.text") 

# plot means matrices on the same plot
ggarrange(g.bc.mus+rremove("xy.text"),g.bbc.mus+rremove("xy.text"),ncol=2,
          labels=c("(a)","(b)"),common.legend = TRUE,font.label = list(size=8))
#==========================================================================================
