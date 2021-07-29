#==========================================================================================
#==========================================================================================
# This is the missing-at-random (MAR) implementation of Bayesian biclustering with 'AgYield' dataset.
#==========================================================================================
#==========================================================================================

library(phyclust);library(reshape2);library(plyr);library(tidyverse);library(Rcpp);library(ggpubr)
setwd(".........")   # set working directory

#==========================================================================================
train_original=as.matrix(readRDS("AgYield.rds"))
I=nrow(train_original);J=ncol(train_original)   # number of rows and columns in dataset
train=train_original

#==========================================================================================
# plot of the dataset
#==========================================================================================
plotdat <- expand.grid(Rows=rownames(data.frame(train_original)), Columns=colnames(data.frame(train_original)))
vec=as.vector(train_original)
plotdat$Z <- vec
color=c(rainbow(10.5,start=0.5,end=0.6)[10.5:1],rainbow(11.5,start=0.05,end=0.20)[11.5:1])
colnames(plotdat)=c("Varieties", "Locations","Z")

g.raw = ggplot(plotdat, aes(x=Locations, y=Varieties)) + geom_tile(aes(fill= Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotdat$Varieties)))) + 
  ggtitle("Raw data matrix") +  font("title", size = 10, face = "bold") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(g.raw + rremove("xy.text"),nrow=1,ncol=1,common.legend = TRUE)

#==========================================================================================
# define required functions
#==========================================================================================
n_func=function(d,drow,dcol,vec){val=length(unlist(d[drow[,J+1]==vec[1],dcol[,I+1]==vec[2]]));return(val)}

mean_func=function(d,drow,dcol,vec)
{
  val=ifelse(n_func(d,drow,dcol,vec)==0,0,
             ifelse(is.nan(mean(unlist(d[drow[,J+1]==vec[1],dcol[,I+1]==vec[2]]),na.rm = TRUE))==TRUE,0,
                    mean(unlist(d[drow[,J+1]==vec[1],dcol[,I+1]==vec[2]]),na.rm = TRUE)))
  return(val)
}

sd_func=function(d,drow,dcol,vec)
{
  val=ifelse(n_func(d,drow,dcol,vec)==0,0,
             ifelse(is.na(sd(unlist(d[drow[,J+1]==vec[1],dcol[,I+1]==vec[2]]),na.rm = TRUE))==TRUE,0,
                    sd(unlist(d[drow[,J+1]==vec[1],dcol[,I+1]==vec[2]]),na.rm = TRUE)))
  return(val)
}

post_var=function(rho,sigma,n)
{
  prior_prec=1/(rho^2)
  data_prec=n/(sigma^2)
  return(1/(prior_prec+data_prec))
}

post_mean=function(rho,mu_prior,sigma,xbar,n)
{
  num=((n*xbar)/(sigma^2))+(mu_prior/(rho^2))
  den=1/post_var(rho,sigma,n)
  return(list(val=rnorm(1,num/den,1/sqrt(den)),m=num/den,s=1/sqrt(den)))
}

stat=function(d,drow,dcol)
{
  m1=matrix(0,nrow=R,ncol=C);m2=matrix(0,nrow=R,ncol=C);m3=matrix(0,nrow=R,ncol=C)
  for(i in 1:R)
  {
    for(j in 1:C)
    {
      m1[i,j]=n_func(d,drow,dcol,c(i,j))
      m2[i,j]=mean_func(d,drow,dcol,c(i,j))
      m3[i,j]=sd_func(d,drow,dcol,c(i,j))
    }
  }
  return(list(n=c(t(m1)),means=c(t(m2)),sds=c(t(m3))))
}

poststat=function(d,drow,dcol)
{
  m1=rep(0,R*C);m2=rep(0,R*C);stats=stat(d,drow,dcol)
  for(i in 1:(R*C))
  {
    m1[i]=post_mean(rho,mu_prior,sigma,stats$means[i],stats$n[i])$val
    m2[i]=post_var(rho,sigma,stats$n[i])
  }
  return(list(postmean=m1,postvar=m2))
}

row_id_update=function(dcol,rownum,reps)
{
  d=matrix(c(dcol[,rownum],dcol[,I+1]),nrow=J,ncol=2)
  mat=matrix(0,nrow=R,ncol=J+1);a=rep(0,R)
  matmean=matrix(post_mean_mat[reps-1,],nrow=R,ncol=C,byrow=T)
  for(i in 1:R)
  {
    for(j in 1:J)
    {
      mat[i,j]=((d[j,1]-matmean[i,d[j,2]])/sigma)^2
    }
    a[i]=(-0.5*sum(mat[i,(1:J)],na.rm = TRUE))
  }
  for(i in 1:R){mat[i,J+1]=(sum(train_row[,J+1][-rownum]==i)+(alpha/R))*exp(a[i]-max(a))}
  vec=rep(0,R);for(i in 1:R){vec[i]=mat[i,J+1]/sum(mat[1:R,J+1])}
  return(vec)
}

col_id_update=function(drow,colnum,reps)
{
  d=matrix(c(drow[,colnum],drow[,J+1]),nrow=I,ncol=2)
  mat=matrix(0,nrow=I+1,ncol=C);a=rep(0,C)
  matmean=matrix(post_mean_mat[reps-1,],nrow=R,ncol=C,byrow=T)
  for(i in 1:C)
  {
    for(j in 1:I)
    {
      mat[j,i]=((d[j,1]-matmean[d[j,2],i])/sigma)^2
    }
    a[i]=(-0.5*sum(mat[(1:I),i],na.rm = TRUE))
  }
  for(i in 1:C){mat[I+1,i]=(sum(train_col[,I+1][-colnum]==i)+(beta/C))*exp(a[i]-max(a))}
  vec=rep(0,C);for(i in 1:C){vec[i]=mat[(I+1),i]/sum(mat[(I+1),1:C])}
  return(vec)
}

row_clust_freq=function()
{
  mat=matrix(nrow=I,ncol=I)
  for(i in 1:I)
  {
    for(j in i:I)
    {
      mat[i,j]=round(sum(id_row_mat[i,][-1]==id_row_mat[j,][-1])/nreps,3)
    }
  }
  return(data.frame(mat))
}

col_clust_freq=function()
{
  mat=matrix(nrow=J,ncol=J)
  for(i in 1:J)
  {
    for(j in i:J)
    {
      mat[i,j]=round(sum(id_col_mat[,i][-1]==id_col_mat[,j][-1])/nreps,3)
    }
  }
  return(data.frame(mat))
}

sse=function(d,drow,dcol)
{
  datmean=matrix(stat(d,drow,dcol)$means,nrow=R,ncol=C,byrow = TRUE)
  mat=matrix(0,nrow=R,ncol=C)
  for(i in 1:R)
  {
    for(j in 1:C)
    {
      mat[i,j]=sum((unlist(d[drow[,J+1]==i,dcol[,I+1]==j])-datmean[i,j])^2,na.rm = TRUE)
    }
  }
  return(sum(as.vector(mat)))
}

row_ri=function(idrowmat)
{
  rowrisucs=rep(0,nreps)
  for(i in 2:(nreps+1)){rowrisucs[i-1]=RRand(idrowmat[,i-1],idrowmat[,i])[[1]]}
  return(RowRIsucs=rowrisucs)
}

col_ri=function(idcolmat)
{
  colrisucs=rep(0,nreps)
  for(i in 2:(nreps+1)){colrisucs[i-1]=RRand(idcolmat[i-1,],idcolmat[i,])[[1]]}
  return(ColRIsucs=colrisucs)
}

element_ri=function()
{
  id_element_mat=matrix(nrow=I*J,ncol=nreps+1)
  Es=matrix(nrow = I*J,ncol=2);Es[,1]=rep(1:I,J);Es[,2]=rep(1:J,each=I)
  for(j in 1:(nreps+1))
  {
    l=length(unique(id_col_mat[j,]))
    for(i in 1:(I*J)){id_element_mat[i,j]=(l*(id_row_mat[Es[i,1],j]-1))+id_col_mat[j,Es[i,2]]}
  }  
  elementrisucs=rep(0,nreps)
  for(i in 2:(nreps+1)){elementrisucs[i-1]=RRand(id_element_mat[,i-1],id_element_mat[,i])[[1]]}
  return(list(id_element_mat=id_element_mat,ElementRIsucs=elementrisucs))
}

sourceCpp("cppFunc.cpp")

data_update=function(drow,dcol,reps)
{
  matmean=matrix(post_mean_mat[reps,],nrow=R,ncol=C,byrow=T);vec=rep(0,nrow(na_index))
  for(i in 1:nrow(na_index))
  {
    vec[i]=rnorm(1,matmean[drow[,J+1][na_index[i,1]],dcol[,I+1][na_index[i,2]]],sigma)
  }
  return(vec)
}

update_train=function(reps)
{
  for(i in 1:nrow(na_index))
  {
    train[unobs_mat[i,1],unobs_mat[i,2]]=unobs_mat[i,reps+2]
    train_row[unobs_mat[i,1],unobs_mat[i,2]]=unobs_mat[i,reps+2]
    train_col[unobs_mat[i,2],unobs_mat[i,1]]=unobs_mat[i,reps+2]
  }
  return(list(t=train,tr=train_row,tc=train_col))
}

prediction=function()
{
  mat=matrix(nrow=I*J,ncol=nreps+3)
  mat[,1]=rep(1:I,J);mat[,2]=rep(1:J,each=I)
  for(reps in 1:(nreps+1))
  {
    for(i in 1:nrow(mat))
    {
      idr=id_row_mat[mat[i,1],reps];idc=id_col_mat[reps,mat[i,2]]
      mat[i,reps+2]=post_mean_mat[reps,((C*(idr-1))+idc)]
    }
  }
  return(mat)
}

log_posterior=function(d,drow,dcol,reps)
{
  vecr=rep(NA,R)
  for(i in 1:R){vecr[i]=ifelse(sum(drow[,J+1]==i)==0,1,prod(seq(alpha/R,(alpha/R)+(sum(drow[,J+1]==i)-1))))}
  gr=log(prod(vecr))
  
  vecc=rep(NA,C)
  for(i in 1:C){vecc[i]=ifelse(sum(dcol[,I+1]==i)==0,1,prod(seq(beta/C,(beta/C)+(sum(dcol[,I+1]==i)-1))))}
  hc=log(prod(vecc))
  
  vec=rep(NA,R*C)
  stats=stat(d,drow,dcol);ms=post_mean(rho,mu_prior,sigma,stats$means,stats$n)
  for(i in 1:(R*C))
  {
    vec[i]=dnorm(post_mean_mat[reps,i],ms$m[i],ms$s[i])
  }
  return(gr+hc+sum(log(vec)))
}

#==========================================================================================
# set parameter values for the sampler
#==========================================================================================
nreps=1e4   # number of iterations T (to be set by user)

sigma=3;mu_prior=0;rho=10   # defined in the paper (to be set by user)

R=15;C=15;alpha=20;beta=20   # defined in the paper (to be set by user)

# initial row and cluster assignment vectors, identical to the ones used with informative missingness assumption
# id_row_0=sample(c(1:R),I,replace = TRUE,prob=rep(1/R,R))
# id_col_0=sample(c(1:C),J,replace = TRUE,prob=rep(1/C,C))
id_row_0=readRDS("idr0.rds");id_col_0=readRDS("idc0.rds")
train_row=cbind(train,idr=id_row_0);train_col=cbind(t(train),idc=id_col_0)

# initial bicluster means
stats=stat(train,train_row,train_col);poststats=poststat(train,train_row,train_col)

#==========================================================================================
# define matrices to store required values
#==========================================================================================
na_index=which(is.na(train),TRUE);unobs_mat=cbind(na_index,matrix(nrow=nrow(na_index),ncol=nreps+1))
n_mat=matrix(nrow=nreps+1,ncol=R*C);mean_mat=matrix(nrow=nreps+1,ncol=R*C);sd_mat=matrix(nrow=nreps+1,ncol=R*C)
post_mean_mat=matrix(nrow=nreps+1,ncol=R*C);post_var_mat=matrix(nrow=nreps+1,ncol=R*C)
id_row_mat=matrix(nrow=I,ncol=nreps+1);id_col_mat=matrix(nrow=nreps+1,ncol=J)
sse_mat=rep(NA,nreps+1);lp=rep(NA,nreps+1)

unobs_mat[,3]=mean(as.vector(train_original),na.rm=TRUE)   # initialize unobserved y's
n_mat[1,]=stats$n;mean_mat[1,]=stats$means;sd_mat[1,]=stats$sds
post_mean_mat[1,]=poststats$postmean;post_var_mat[1,]=poststats$postvar
id_row_mat[,1]=id_row_0;id_col_mat[1,]=id_col_0
sse_mat[1]=sse(train,train_row,train_col);lp[1]=log_posterior(train,train_row,train_col,1)

trainlist=update_train(1)   
train=trainlist$t;train_row=trainlist$tr;train_col=trainlist$tc

#==========================================================================================
# the MCMC sampler
#==========================================================================================
for(reps in 2:(nreps+1))
{
  # row id updates
  row_id_update_prob_mat=matrix(0,nrow=I,ncol=R+1)
  for(i in 1:I)
  {
    row_id_update_prob_mat[i,1:R]=row_id_update(train_col,i,reps)
    row_id_update_prob_mat[i,R+1]=sample(c(1:R),1,prob=row_id_update_prob_mat[i,1:R])
    train_row[i,J+1]=row_id_update_prob_mat[i,R+1]
  }
  id_row_mat[,reps]=train_row[,J+1]
  
  # column id updates
  col_id_update_prob_mat=matrix(0,nrow=J,ncol=C+1)
  for(j in 1:J)
  {
    col_id_update_prob_mat[j,1:C]=col_id_update(train_row,j,reps)
    col_id_update_prob_mat[j,C+1]=sample(c(1:C),1,prob=col_id_update_prob_mat[j,1:C])
    train_col[j,I+1]=col_id_update_prob_mat[j,C+1]
  }
  id_col_mat[reps,]=train_col[,I+1]
  
  # update bicluster means
  stats=stat(train,train_row,train_col)
  n_mat[reps,]=stats$n
  mean_mat[reps,]=stats$means
  sd_mat[reps,]=stats$sds

  poststats=poststat(train,train_row,train_col)
  post_mean_mat[reps,]=poststats$postmean
  post_var_mat[reps,]=poststats$postvar
  
  # update unobserved y's
  unobs_mat[,reps+2]=data_update(train_row,train_col,reps)

  trainlist=update_train(reps)
  train=trainlist$t
  train_row=trainlist$tr
  train_col=trainlist$tc
  
  lp[reps]=log_posterior(train,train_row,train_col,reps)
  
  sse_mat[reps]=sse(train_original,train_row,train_col)
}

#==========================================================================================
# obtaining results
#==========================================================================================
# compute cluster frequencies
rcf=row_clust_freq();ccf=col_clust_freq()

# compute Rand indices
rowRI=row_ri(id_row_mat);colRI=col_ri(id_col_mat);elementRI=element_ri()

# C++ implementation to compute average Rand indices
mr=arr(id_row_mat,nreps,RRand);mc=acr(id_col_mat,nreps,RRand);me=aer(elementRI$id_element_mat,nreps,RRand)
mr[lower.tri(mr)]=t(mr)[lower.tri(mr)];mc[lower.tri(mc)]=t(mc)[lower.tri(mc)];me[lower.tri(me)]=t(me)[lower.tri(me)]
diag(mr)=NA;diag(mc)=NA;diag(me)=NA
avgRowRI=apply(mr,2,function(x) mean(x,na.rm = TRUE));avgColRI=apply(mc,2,function(x) mean(x,na.rm = TRUE))
avgElementRI=apply(me,2,function(x) mean(x,na.rm = TRUE))

sse_dat=data.frame(Iteration=c(1:(nreps+1)),SSE=sse_mat)

lp_dat=data.frame(Iteration=c(1:(nreps+1)),log_posterior=lp)

RI=data.frame(Iteration=c(1:nreps), RowRIsuccessive=rowRI, ColumnRIsuccessive=colRI, 
              ElementRIsuccessive=elementRI$ElementRIsucs)
RI=melt(RI, id.vars = "Iteration")

#==========================================================================================
# obtaining plots
#==========================================================================================
# total SSE plot
ggplot(data=sse_dat, aes(x=Iteration, y=SSE)) + geom_line(color="blue", size=1) +
  ggtitle("SSE Plot") + theme(plot.title = element_text(hjust = 0.5))

# log of joint posterior plot
ggplot(data=lp_dat, aes(x=Iteration, y=log_posterior)) + geom_line(color="blue", size=1) +
  ggtitle("Log Posterior Plot") + theme(plot.title = element_text(hjust = 0.5))

# plots of Rand indices
ggplot(data=RI, aes(x=Iteration, y=value, group=variable)) +
  geom_line(aes(color=variable), size=0.7) +
  geom_point(aes(color=variable), size=1) + facet_grid(variable~.) +
  theme(legend.position = "none")

ggplot(data=RI, aes(x=value)) + geom_histogram() + facet_grid(variable~.) + 
  theme_bw() + ggtitle("Histograms of Rand Indices") + theme(plot.title = element_text(hjust = 0.5))

#==========================================================================================
# plot biclusters (reported only for max average RI criterion)
# to replicate for min total SSE and max log posterior criteria, please consult master code 'BayesianBiclustering_WithMissing.R'

# maximum average Rand index
train_clust_2=train_original[order(id_row_mat[,which.max(avgRowRI)]),order(id_col_mat[which.max(avgColRI),])]
train_clust_plot_2=expand.grid(Rows=rownames(data.frame(train_clust_2)), Columns=colnames(data.frame(train_clust_2)))
val=as.vector(train_clust_2)
train_clust_plot_2$Z=val
v=cumsum(as.vector(table(id_col_mat[which.max(avgColRI),])))+0.5;vline_coords_2=data.frame(v=v[-length(v)])
h=cumsum(as.vector(table(id_row_mat[,which.max(avgRowRI)])))+0.5;hline_coords_2=data.frame(h=h[-length(h)])
clust_plot_2=list(data=train_clust_plot_2,vlines=vline_coords_2,hlines=hline_coords_2)

colnames(clust_plot_2$data)=c("Varieties", "Locations","Z")
color=c(rainbow(10.5,start=0.5,end=0.6)[10.5:1],rainbow(11.5,start=0.05,end=0.20)[11.5:1])
g2.MAR = ggplot(data=clust_plot_2$data, aes(x=Locations, y=Varieties)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=clust_plot_2$vlines,aes(xintercept = v)) +
  geom_hline(data=clust_plot_2$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Biclusters maximizing average Row and Column RI (MAR)") + rremove("xy.text") +
  theme(plot.title = element_text(hjust = 0.5)) +  font("title", size = 10, face = "bold") 
g2.MAR

#==========================================================================================
# generate predicted values for unobserved/missing y's
# predictions are generated from average of the corresponding posterior means ("prediction" function) 
# across all iterations
#==========================================================================================
vec=as.vector(train_original)
full_mat=matrix(nrow=length(vec),ncol=3)
full_mat[,1]=rep(1:I,J);full_mat[,2]=rep(1:J,each=I);full_mat[,3]=vec
full_mat=data.frame(full_mat)
colnames(full_mat)=c("row","col","observed")
pred_mat=prediction()
B=0   # burn-in iterations (change as needed)
pred=data.frame(cbind(pred_mat[,1:2],apply(pred_mat[,-(1:2)],1,mean)))
colnames(pred)=c("row","col","predicted_from_mu")
full_dat=join(full_mat,pred,by=c("row","col"),type="left")

observed_pred=ifelse(is.na(full_dat$observed),full_dat$predicted_from_mu,full_dat$observed)
mat_observed_pred_MAR=data.frame(matrix(observed_pred,nrow = I,ncol = J))

#==========================================================================================
# plot biclusters with predictions (reported only for max average RI criterion)
# to replicate for min total SSE and max log posterior criteria, please consult master code 'BayesianBiclustering_WithMissing.R'

train_pred_clust_5=mat_observed_pred_MAR[order(id_row_mat[,which.max(avgRowRI)]),order(id_col_mat[which.max(avgColRI),])]
train_pred_clust_plot_5 <- expand.grid(Rows=rownames(train_pred_clust_5), Columns=colnames(train_pred_clust_5))
vec=as.vector(as.matrix(train_pred_clust_5))
train_pred_clust_plot_5$Z <- vec;colnames(train_pred_clust_plot_5)=c("Varieties","Locations","Z")
pred_clust_plot_5=list(data=train_pred_clust_plot_5,vlines=vline_coords_2,hlines=hline_coords_2)

g5.predMAR=ggplot(data=pred_clust_plot_5$data, aes(x=Locations, y=Varieties)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=pred_clust_plot_5$vlines,aes(xintercept = v)) +
  geom_hline(data=pred_clust_plot_5$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Biclusters with predictions for max average RI (MAR)") + 
  font("title", size = 10, face = "bold") +
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text")

ggarrange(g2.MAR, g5.predMAR,nrow=2,ncol=1,labels=c("(a)","(b)"),common.legend = TRUE)

#==========================================================================================
# obtain summary and plot of numbers of non-empty row-column cluster combinations
#==========================================================================================
pq=data.frame(Iteration=c(1:nreps),p=apply(id_row_mat[,-1],2,function(x) length(unique(x))),
              q=apply(id_col_mat[-1,],1,function(x) length(unique(x))),
              SSE=sse_mat[-1],log_posterior=lp[-1])
pq_summary=ddply(pq,.(p,q),summarize,freq=length(SSE),meanSSE=mean(SSE),minSSE=min(SSE),maxSSE=max(SSE),
                 mean_log_post=mean(log_posterior,na.rm = TRUE),min_log_post=min(log_posterior,na.rm = TRUE),
                 max_log_post=max(log_posterior,na.rm = TRUE))

pq_plot=melt(select(pq,c("Iteration","p","q")), id.vars = "Iteration")

ggplot(data=pq_plot, aes(x=Iteration, y=value, group=variable)) +
  geom_line(aes(color=variable), size=0.7) +
  geom_point(aes(color=variable), size=1) + facet_grid(variable~.) +
  ggtitle("p-q Plot") + theme(plot.title = element_text(hjust = 0.5))

#==========================================================================================
save.image(file="AgYield_MAR.RData")
#==========================================================================================