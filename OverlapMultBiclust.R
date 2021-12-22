#==========================================================================================
#==========================================================================================
# This is additional code for the simulation study on overlapping multiplicative biclusters setting (section 6.5)
#==========================================================================================
#==========================================================================================

library(reshape2);library(tidyverse);library(ggpubr)

#==========================================================================================
# dataset generation
#==========================================================================================
set.seed(5)
u1<-c(10,9,8,7,6,5,4,3,rep(2,17),rep(0,75))
v1<-c(10,-10,8,-8,5,-5,rep(3,5),rep(-3,5),rep(0,34))
u1<-u1/sqrt(sum(u1^2))
v1<-v1/sqrt(sum(v1^2))
u2<-c(rep(0,13),10,9,8,7,6,5,4,3,rep(2,17),rep(0,62))
v2<-c(rep(0,9),10,-9,8,-7,6,-5,rep(4,5),rep(-3,5),rep(0,25))
u2<-u2/sqrt(sum(u2^2))
v2<-v2/sqrt(sum(v2^2))
d<-50
mus<-(d*u1%*%t(v1))+(d*u2%*%t(v2))
binaryX<-(mus!=0)*1
X<-mus+matrix(rnorm(100*50),100,50)
X<-X-mean(X)

saveRDS(X,"OverlapMultBiclust_X.rds")

#==========================================================================================
# Perform Bayesian biclustering (with real dataset OverlapMultBiclust_X.rds) using either the first 
# or second version (see master codes) of the MCMC sampler. Save the results and plots obtained.
# The following lines of code can be used to plot the means matrices.
#==========================================================================================
setwd(".........")   # set working directory
load("....RData")   # load saved results

# plot true underlying means matrix
mus = data.frame(mus)
plotmus = expand.grid(Rows = rownames(mus), Columns = colnames(mus))
vec = as.vector(as.matrix(mus))
plotmus$Z = vec

color = c(rainbow(15,start=0.5,end=0.6)[15:1],"white",rainbow(15,start=0.1,end=0.2)[15:1])
g.mus = ggplot(plotmus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotmus$Rows)))) + 
  ggtitle("True underlying means matrix") + 
  theme(plot.title = element_text(hjust = 0.5))
g.mus + rremove("xy.text")

# plot means matrix estimated from Bayesian biclustering
bbc.mus = matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  for(j in 1:J){
    # change criterion (see paper and master codes)
    bbc.mus[i,j] = mean(train[id_row_mat[,which.max(avgRowRI)]==id_row_mat[i,which.max(avgRowRI)],
                              id_col_mat[which.max(avgColRI),]==id_col_mat[which.max(avgColRI),j]])		
  }
}

bbc.mus = round(data.frame(bbc.mus),2)
plotbbc.mus = expand.grid(Rows = rownames(bbc.mus), Columns = colnames(bbc.mus))
vec = as.vector(as.matrix(bbc.mus))
plotbbc.mus$Z = vec

color = c(rainbow(15,start=0.5,end=0.6)[15:1],"white",rainbow(15,start=0.1,end=0.2)[15:1])
g.bbc.mus = ggplot(plotbbc.mus, aes(x = Columns, y = Rows)) + geom_tile(aes(fill = Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotbbc.mus$Rows)))) + 
  ggtitle("Estimated means matrix") + 
  theme(plot.title = element_text(hjust = 0.5))
g.bbc.mus + rremove("xy.text") 

# plot both means matrices on the same plot
ggarrange(g.mus+rremove("xy.text"), g.bbc.mus+rremove("xy.text"),
          labels=c("(a)","(b)"),ncol=2,common.legend = TRUE,font.label = list(size=8))
#==========================================================================================
