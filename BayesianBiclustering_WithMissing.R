#==========================================================================================
#==========================================================================================
# This is the master code for the MCMC sampler for datasets with missing/unobserved entries
#==========================================================================================
#==========================================================================================

library(phyclust);library(reshape2);library(plyr);library(tidyverse);library(ggpubr);library(Rcpp)
setwd(".........")   # set working directory

#==========================================================================================
# simulate a rectangular dataset with missing entries from a given means matrix, variance, "true" row and column blocks
# (to be used only for simulation)
#==========================================================================================
# # set.seed(333)
# 
# I=20;J=20   # total number of rows and columns in dataset
# R=2;C=2     # "true" number of row and column blocks (to be used only to simulate a dataset)
# sigma=1     # data variance
# 
# truthCs<-rep(1:R, each=(I/R));truthDs<-rep(1:C, each=(J/C))
# Es=matrix(nrow = I*J,ncol=2);Es[,1]=rep(1:I,J);Es[,2]=rep(1:J,each=I)
# 
# # generate means matrix
# # lb=-3;ub=3;mus<-runif(R*C,lb,ub)   # randomly choose means between a lower and upper bound
# mus=c(-3,3,-10,10)   # specify means
# mus<-matrix(c(mus),nrow=R,ncol=C,byrow=TRUE);ids=matrix(c(1:(R*C)),nrow=R,ncol=C,byrow=TRUE)
# 
# # generate observations
# x<-matrix(rnorm(I*J,mean=0,sd=sigma),nrow=I,ncol=J)
# musmatrix<-matrix(NA,nrow=I,ncol=J);idsmatrix=matrix(nrow=I,ncol=J)
# for(i in 1:max(truthCs)){
#   for(j in 1:max(truthDs)){ 
#     x[truthCs==i,truthDs==j]<-x[truthCs==i,truthDs==j]+mus[i,j]
#     musmatrix[truthCs==i,truthDs==j]<-mus[i,j]
#     idsmatrix[truthCs==i,truthDs==j]<-ids[i,j]
#   } 
# }	
# truthEs=rep(NA,I*J);for(i in 1:(I*J)){truthEs[i]=idsmatrix[Es[i,1],Es[i,2]]}
# 
# # generate missing/unobserved y's
# p=function(z,l1,l2){l=exp((z-l1)/l2)/(1+exp((z-l1)/l2));return(l)}
# val=as.vector(x)
# p_obs=p(val,summary(as.vector(mus))[[4]],5);p_missing=1-p_obs
# # p_obs=p(val,lb-sd(as.vector(mus)),1);p_missing=1-p_obs
# # p_obs=p(val,ub+sd(as.vector(mus)),1);p_missing=1-p_obs
# 
# dat=matrix(nrow=length(val),ncol=6)
# dat[,1]=rep(1:I,J);dat[,2]=rep(1:J,each=I);dat[,3]=val;dat[,4]=p_obs;dat[,5]=p_missing
# for(i in 1:length(val)){dat[i,6]=sample(c(dat[i,3],NA),size = 1,prob = c(dat[i,4],dat[i,5]))}
# x=matrix(dat[,6],nrow=I,ncol=J)
# 
# # randomly rearrange rows and columns
# samp_row=sample(1:I,I);samp_col=sample(1:J,J)
# train_original=x[samp_row,samp_col]   # original dataset to be worked with
# train=train_original   
# 
# # "true" row, column, element-wise clusterings
# tCs=truthCs[samp_row];tDs=truthDs[samp_col]
# tEs=rep(NA,I*J);for(i in 1:(I*J)){tEs[i]=idsmatrix[samp_row,samp_col][Es[i,1],Es[i,2]]}
# 
# 100*(sum(is.na(train))/(I*J))   # % of missing entries

#==========================================================================================
# when working with a real dataset (load the data array)
#==========================================================================================
train_original=as.matrix(readRDS("....rds"))   # dataset is in .rds format and must be loaded as a matrix 
I=nrow(train_original);J=ncol(train_original)   # number of rows and columns in dataset
train=train_original

100*(sum(is.na(train))/(I*J))   # % of missing entries

#==========================================================================================
# plot of the data array
#==========================================================================================
plotdat <- expand.grid(Rows=rownames(data.frame(train_original)), Columns=colnames(data.frame(train_original)))
vec=as.vector(train_original)
plotdat$Z <- vec

color=c(rainbow(8,start=0.5,end=0.6)[8:1],rainbow(17,start=0.05,end=0.20)[17:1])   # can be specified by user
g.data = ggplot(plotdat, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(as.factor(plotdat$Rows)))) + ggtitle("Raw data Matrix") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.data

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
  return(RowRIsucs=rowrisucs)   # comment out when using following lines of code
  # following lines of code to be used only when "true" row clustering 'tCs' is known
  # rowritrue=rep(0,nreps)
  # for(i in 2:(nreps+1)){rowritrue[i-1]=RRand(tCs,idrowmat[,i])[[1]]}   
  # return(list(RowRIsucs=rowrisucs,RowRItrue=rowritrue))
}

col_ri=function(idcolmat)
{
  colrisucs=rep(0,nreps)
  for(i in 2:(nreps+1)){colrisucs[i-1]=RRand(idcolmat[i-1,],idcolmat[i,])[[1]]}
  return(ColRIsucs=colrisucs)   # comment out when using following lines of code
  # following lines of code to be used only when "true" column clustering 'tDs' is known
  # colritrue=rep(0,nreps)
  # for(i in 2:(nreps+1)){colritrue[i-1]=RRand(tDs,idcolmat[i,])[[1]]}
  # return(list(ColRIsucs=colrisucs,ColRItrue=colritrue))
}

element_ri=function()
{
  id_element_mat=matrix(nrow=I*J,ncol=nreps+1)
  # comment out following line of code if working with a simulated dataset
  Es=matrix(nrow = I*J,ncol=2);Es[,1]=rep(1:I,J);Es[,2]=rep(1:J,each=I)
  for(j in 1:(nreps+1))
  {
    l=length(unique(id_col_mat[j,]))
    for(i in 1:(I*J)){id_element_mat[i,j]=(l*(id_row_mat[Es[i,1],j]-1))+id_col_mat[j,Es[i,2]]}
  }  
  elementrisucs=rep(0,nreps)
  for(i in 2:(nreps+1)){elementrisucs[i-1]=RRand(id_element_mat[,i-1],id_element_mat[,i])[[1]]}
  return(list(id_element_mat=id_element_mat,ElementRIsucs=elementrisucs))   # comment out when using following lines of code
  # following lines of code to be used only when "true" element-wise clustering 'tEs' is known
  # elementritrue=rep(0,nreps)
  # for(i in 2:(nreps+1)){elementritrue[i-1]=RRand(tEs,id_element_mat[,i])[[1]]}
  # return(list(id_element_mat=id_element_mat,ElementRIsucs=elementrisucs,ElementRItrue=elementritrue))
}

sourceCpp("cppFunc.cpp")   # C++ implementation to compute average Rand indices

# following three functions are R implementations to compute average Rand indices
# average_row_ri=function()
# {
#   mat=matrix(NA,nrow=nreps+1,ncol=nreps+1)
#   for(i in 1:(nrow(mat)-1))
#   {
#     for(j in (i+1):ncol(mat))
#     {
#       mat[i,j]=RRand(id_row_mat[,i],id_row_mat[,j])[[1]]
#     }
#   }
#   mat[lower.tri(mat)]=t(mat)[lower.tri(mat)]
#   vec=apply(mat,2,function(x) mean(x,na.rm = TRUE))
#   return(vec)
# }
# average_col_ri=function()
# {
#   mat=matrix(NA,nrow=nreps+1,ncol=nreps+1)
#   for(i in 1:(nrow(mat)-1))
#   {
#     for(j in (i+1):ncol(mat))
#     {
#       mat[i,j]=RRand(id_col_mat[i,],id_col_mat[j,])[[1]]
#     }
#   }
#   mat[lower.tri(mat)]=t(mat)[lower.tri(mat)]
#   vec=apply(mat,2,function(x) mean(x,na.rm = TRUE))
#   return(vec)
# }
# average_element_ri=function()
# {
#   iem=elementRI$id_element_mat
#   mat=matrix(NA,nrow=nreps+1,ncol=nreps+1)
#   for(i in 1:(nrow(mat)-1))
#   {
#     for(j in (i+1):ncol(mat))
#     {
#       mat[i,j]=RRand(iem[,i],iem[,j])[[1]]
#     }
#   }
#   mat[lower.tri(mat)]=t(mat)[lower.tri(mat)]
#   vec=apply(mat,2,function(x) mean(x,na.rm = TRUE))
#   return(vec)
# }

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

log_density_lambda1_lambda2=function(l1,k)
{
  y_all=as.vector(train)
  y_seen=as.vector(train_original);y_seen=y_seen[!is.na(y_seen)]
  num=sum((y_seen-l1)/exp(k))+(-(l1^2)/(2*(rho_lambda1^2)))+(-(k^2)/(2*(rho_lambda2^2)))
  den=sum(log(1+exp((y_all-l1)/exp(k))))
  return(num-den)
}

update_lambda1=function(l1_current,k)
{
  l1_proposed=rnorm(1,l1_current,gamma_lambda1)
  log_accept_ratio=min(0,log_density_lambda1_lambda2(l1_proposed,k)-log_density_lambda1_lambda2(l1_current,k))
  u=runif(1,0,1)
  val=ifelse(log(u) < log_accept_ratio,l1_proposed,l1_current)
  return(val)
}

update_k=function(l1,k_current)
{
  k_proposed=rnorm(1,k_current,gamma_lambda2)
  log_accept_ratio=min(0,log_density_lambda1_lambda2(l1,k_proposed)-log_density_lambda1_lambda2(l1,k_current))
  u=runif(1,0,1)
  val=ifelse(log(u) < log_accept_ratio,k_proposed,k_current)
  return(val)
}

density_y=function(y,m,l1,l2){d=exp(-((y-m)^2)/(2*(sigma^2)))/(1+exp((y-l1)/l2));return(d)}

update_y=function(reps)
{
  mat=matrix(post_mean_mat[reps,],nrow=R,ncol=C,byrow=TRUE);vec=rep(0,nrow(na_index))
  for(i in 1:nrow(na_index))
  {
    y_current=unobs_mat[i,reps+2-1]
    r_i=id_row_mat[na_index[i,1],reps];c_j=id_col_mat[reps,na_index[i,2]]
    y_proposed=rnorm(1,y_current,gamma_y)
    accept_ratio=min(1,density_y(y_proposed,mat[r_i,c_j],lambda1[reps],lambda2[reps])/density_y(y_current,mat[r_i,c_j],lambda1[reps],lambda2[reps]))
    u=runif(1,0,1)
    vec[i]=ifelse(u < accept_ratio,y_proposed,y_current)
  } 
  return(vec)
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
  
  l1_l2=log_density_lambda1_lambda2(lambda1[reps],k[reps])
  
  return(gr+hc+sum(log(vec))+l1_l2)
}

#==========================================================================================
# set parameter values for the sampler
#==========================================================================================
nreps=1e4   # number of iterations T (to be set by user)

# defined in the paper (to be set by user)
sigma=1;mu_prior=0;rho=10
rho_lambda1=10;rho_lambda2=10
gamma_y=4;gamma_lambda1=2;gamma_lambda2=1                     

R=4;C=4;alpha=10;beta=10   # defined in the paper (to be set by user)

# initial row and column cluster assignment vectors, lambda1 and lambda2
id_row_0=sample(c(1:R),I,replace = TRUE,prob=rep(1/R,R))
id_col_0=sample(c(1:C),J,replace = TRUE,prob=rep(1/C,C))
train_row=cbind(train,idr=id_row_0);train_col=cbind(t(train),idc=id_col_0)
lambda1_0=0;lambda2_0=1   

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
lambda1=rep(NA,nreps+1);lambda2=rep(NA,nreps+1);k=log(lambda2)

unobs_mat[,3]=mean(as.vector(train_original),na.rm=TRUE)   # initialize unobserved y's
n_mat[1,]=stats$n;mean_mat[1,]=stats$means;sd_mat[1,]=stats$sds
post_mean_mat[1,]=poststats$postmean;post_var_mat[1,]=poststats$postvar
id_row_mat[,1]=id_row_0;id_col_mat[1,]=id_col_0
sse_mat[1]=sse(train,train_row,train_col);lp[1]=log_posterior(train,train_row,train_col,1)
lambda1[1]=lambda1_0;lambda2[1]=lambda2_0;k[1]=log(lambda2[1])

trainlist=update_train(1) 
train=trainlist$t;train_row=trainlist$tr;train_col=trainlist$tc

#==========================================================================================
# the MCMC sampler
#==========================================================================================
for(reps in 2:(nreps+1))
{
  # update lambda1 and log(lambda2) as in step 5(i)
  lambda1[reps]=update_lambda1(lambda1[reps-1],log(lambda2[reps-1]))
  k[reps]=update_k(lambda1[reps],k[reps-1])
  lambda2[reps]=exp(k[reps])
  
  # row id updates as in step 5(ii)
  row_id_update_prob_mat=matrix(0,nrow=I,ncol=R+1)
  for(i in 1:I)
  {
    row_id_update_prob_mat[i,1:R]=row_id_update(train_col,i,reps)
    row_id_update_prob_mat[i,R+1]=sample(c(1:R),1,prob=row_id_update_prob_mat[i,1:R])
    train_row[i,J+1]=row_id_update_prob_mat[i,R+1]
  }
  id_row_mat[,reps]=train_row[,J+1]

  # column id updates as in step 5(iii)
  col_id_update_prob_mat=matrix(0,nrow=J,ncol=C+1)
  for(j in 1:J)
  {
    col_id_update_prob_mat[j,1:C]=col_id_update(train_row,j,reps)
    col_id_update_prob_mat[j,C+1]=sample(c(1:C),1,prob=col_id_update_prob_mat[j,1:C])
    train_col[j,I+1]=col_id_update_prob_mat[j,C+1]
  }
  id_col_mat[reps,]=train_col[,I+1]

  # update bicluster means as in step 5(iv)
  stats=stat(train,train_row,train_col)
  n_mat[reps,]=stats$n
  mean_mat[reps,]=stats$means
  sd_mat[reps,]=stats$sds
  
  poststats=poststat(train,train_row,train_col)
  post_mean_mat[reps,]=poststats$postmean
  post_var_mat[reps,]=poststats$postvar
  
  # update unobserved y's as in step 5(v)
  unobs_mat[,reps+2]=update_y(reps)
  
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

# following three operations compute average Rand indices (slower)
# avgRowRI=average_row_ri();avgColRI=average_col_ri();avgElementRI=average_element_ri()

# below is a C++ implementation to compute average Rand indices (faster)
mr=arr(id_row_mat,nreps,RRand);mc=acr(id_col_mat,nreps,RRand);me=aer(elementRI$id_element_mat,nreps,RRand)
mr[lower.tri(mr)]=t(mr)[lower.tri(mr)];mc[lower.tri(mc)]=t(mc)[lower.tri(mc)];me[lower.tri(me)]=t(me)[lower.tri(me)]
diag(mr)=NA;diag(mc)=NA;diag(me)=NA
avgRowRI=apply(mr,2,function(x) mean(x,na.rm = TRUE));avgColRI=apply(mc,2,function(x) mean(x,na.rm = TRUE))
avgElementRI=apply(me,2,function(x) mean(x,na.rm = TRUE))

# following operation to be modified based on whether the dataset is simulated or real
# see 'row_ri', 'column_ri', 'element_ri' functions
RI=data.frame(Iteration=c(1:nreps), 
              # RowRItrue=rowRI$RowRItrue, RowRIsuccessive=rowRI$RowRIsucs, 
              RowRIsuccessive=rowRI, 
              # ColumnRItrue=colRI$ColRItrue, ColumnRIsuccessive=colRI$ColRIsucs,
              ColumnRIsuccessive=colRI,
              # ElementRItrue=elementRI$ElementRItrue,
              ElementRIsuccessive=elementRI$ElementRIsucs)

# summaries of numbers of non-empty row-column cluster combinations
pq=data.frame(Iteration=c(1:nreps),
              p=apply(id_row_mat[,-1],2,function(x) length(unique(x))),
              q=apply(id_col_mat[-1,],1,function(x) length(unique(x))),
              SSE=sse_mat[-1],
              log_posterior=lp[-1]) %>%
  mutate(pSSE = log10(SSE+1)+p*q)   # calculate penalized SSE

pq_summary=ddply(pq,.(p,q),summarize,freq=length(SSE),
                 meanSSE=mean(SSE),minSSE=min(SSE),maxSSE=max(SSE),
                 mean_log_post=mean(log_posterior,na.rm = TRUE),
                 min_log_post=min(log_posterior,na.rm = TRUE),
                 max_log_post=max(log_posterior,na.rm = TRUE),
                 mean_pSSE=mean(pSSE),min_pSSE=min(pSSE),max_pSSE=max(pSSE))

#==========================================================================================
# obtaining plots
#==========================================================================================
# total SSE plot
sse_dat=data.frame(Iteration=c(1:(nreps+1)),SSE=sse_mat)

g.sse = ggplot(data=sse_dat, aes(x=Iteration, y=SSE)) + geom_line(color="blue", size=1) +
  ggtitle("SSE Plot") + theme(plot.title = element_text(hjust = 0.5))
g.sse

# log of joint posterior plot
lp_dat=data.frame(Iteration=c(1:(nreps+1)),log_posterior=lp)

g.lp = ggplot(data=lp_dat, aes(x=Iteration, y=log_posterior)) + geom_line(color="blue", size=1) +
  ggtitle("Log Posterior Plot") + theme(plot.title = element_text(hjust = 0.5))
g.lp

# plots of Rand indices
RI=melt(RI, id.vars = "Iteration")

g.ri = ggplot(data=RI, aes(x=Iteration, y=value, group=variable)) +
  geom_line(aes(color=variable), size=0.7) +
  geom_point(aes(color=variable), size=1) + facet_grid(variable~.) +
  theme(legend.position = "none")
g.ri

g.ri.hist = ggplot(data=RI, aes(x=value)) + geom_histogram() + facet_grid(variable~.) + 
  theme_bw() + ggtitle("Histograms of Rand Indices") + theme(plot.title = element_text(hjust = 0.5))
g.ri.hist

# plot of numbers of non-empty row-column cluster combinations
pq_plot=melt(select(pq,c("Iteration","p","q")), id.vars = "Iteration")

g.pq = ggplot(data=pq_plot, aes(x=Iteration, y=value, group=variable)) +
  geom_line(aes(color=variable), size=0.7) +
  geom_point(aes(color=variable), size=1) + facet_grid(variable~.) +
  ggtitle("p-q Plot") + theme(plot.title = element_text(hjust = 0.5))
g.pq

#==========================================================================================
# plot biclusters, biclusters determined by 4 criteria

# minimum total SSE 
train_clust_1=train_original[order(id_row_mat[,which.min(sse_mat[-1])]),order(id_col_mat[which.min(sse_mat[-1]),])]
train_clust_plot_1=expand.grid(Rows=rownames(data.frame(train_clust_1)), Columns=colnames(data.frame(train_clust_1)))
val=as.vector(train_clust_1)
train_clust_plot_1$Z=val
v=cumsum(as.vector(table(id_col_mat[which.min(sse_mat[-1]),])))+0.5;vline_coords_1=data.frame(v=v[-length(v)])
h=cumsum(as.vector(table(id_row_mat[,which.min(sse_mat[-1])])))+0.5;hline_coords_1=data.frame(h=h[-length(h)])
clust_plot_1=list(data=train_clust_plot_1,vlines=vline_coords_1,hlines=hline_coords_1)

g.1 = ggplot(data=clust_plot_1$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=clust_plot_1$vlines,aes(xintercept = v)) +
  geom_hline(data=clust_plot_1$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_clust_plot_1$Rows)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Biclusters identified for min total SSE") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.1

# maximum average Rand index
train_clust_2=train_original[order(id_row_mat[,which.max(avgRowRI)]),order(id_col_mat[which.max(avgColRI),])]
train_clust_plot_2=expand.grid(Rows=rownames(data.frame(train_clust_2)), Columns=colnames(data.frame(train_clust_2)))
val=as.vector(train_clust_2)
train_clust_plot_2$Z=val
v=cumsum(as.vector(table(id_col_mat[which.max(avgColRI),])))+0.5;vline_coords_2=data.frame(v=v[-length(v)])
h=cumsum(as.vector(table(id_row_mat[,which.max(avgRowRI)])))+0.5;hline_coords_2=data.frame(h=h[-length(h)])
clust_plot_2=list(data=train_clust_plot_2,vlines=vline_coords_2,hlines=hline_coords_2)

g.2 = ggplot(data=clust_plot_2$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=clust_plot_2$vlines,aes(xintercept = v)) +
  geom_hline(data=clust_plot_2$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_clust_plot_2$Rows)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Biclusters identified for max average RI") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.2

# maximum log of joint posterior
train_clust_3=train_original[order(id_row_mat[,which.max(lp[-1])+1]),order(id_col_mat[which.max(lp[-1])+1,])]
train_clust_plot_3=expand.grid(Rows=rownames(data.frame(train_clust_3)), Columns=colnames(data.frame(train_clust_3)))
val=as.vector(train_clust_3)
train_clust_plot_3$Z=val
v=cumsum(as.vector(table(id_col_mat[which.max(lp[-1])+1,])))+0.5;vline_coords_3=data.frame(v=v[-length(v)])
h=cumsum(as.vector(table(id_row_mat[,which.max(lp[-1])+1])))+0.5;hline_coords_3=data.frame(h=h[-length(h)])
clust_plot_3=list(data=train_clust_plot_3,vlines=vline_coords_3,hlines=hline_coords_3)

g.3 = ggplot(data=clust_plot_3$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=clust_plot_3$vlines,aes(xintercept = v)) +
  geom_hline(data=clust_plot_3$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_clust_plot_3$Rows)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Biclusters identified for max log posterior") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.3

# minimum penalized SSE 
train_clust_4=train_original[order(id_row_mat[,which.min(pq$pSSE)+1]),order(id_col_mat[which.min(pq$pSSE)+1,])]
train_clust_plot_4=expand.grid(Rows=rownames(data.frame(train_clust_4)), Columns=colnames(data.frame(train_clust_4)))
val=as.vector(train_clust_4)
train_clust_plot_4$Z=val
v=cumsum(as.vector(table(id_col_mat[which.min(pq$pSSE)+1,])))+0.5;vline_coords_4=data.frame(v=v[-length(v)])
h=cumsum(as.vector(table(id_row_mat[,which.min(pq$pSSE)+1])))+0.5;hline_coords_4=data.frame(h=h[-length(h)])
clust_plot_4=list(data=train_clust_plot_4,vlines=vline_coords_4,hlines=hline_coords_4)

g.4 = ggplot(data=clust_plot_4$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=clust_plot_4$vlines,aes(xintercept = v)) +
  geom_hline(data=clust_plot_4$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color, na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_clust_plot_4$Rows)))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Biclusters identified for min penalized SSE") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.4

#==========================================================================================
# row cluster frequency plot
rcf[lower.tri(rcf)]=t(rcf)[lower.tri(rcf)]
rcfdat <- expand.grid(Rows = rownames(rcf), Columns = colnames(rcf))
vec = as.vector(as.matrix(rcf))
rcfdat$Z <- vec

g.rcf = ggplot(data = rcfdat, aes(Rows, Columns, fill = Z))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", limit = c(0,1)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed() +
  scale_y_discrete(limits = rev(levels(as.factor(rcfdat$Columns)))) +
  xlab("Rows") + ylab("Rows") + ggtitle("Row Cluster Frequency") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")
g.rcf

# column cluster frequency plot
ccf[lower.tri(ccf)]=t(ccf)[lower.tri(ccf)]
ccfdat <- expand.grid(Rows = rownames(ccf), Columns = colnames(ccf))
vec = as.vector(as.matrix(ccf))
ccfdat$Z <- vec

g.ccf = ggplot(data = ccfdat, aes(Rows, Columns, fill = Z))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", limit = c(0,1)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1))+
  coord_fixed() +
  scale_y_discrete(limits = rev(levels(as.factor(ccfdat$Columns)))) +
  xlab("Columns") + ylab("Columns") + ggtitle("Column Cluster Frequency") + 
  theme(plot.title = element_text(hjust = 0.5)) + rremove("xylab")
g.ccf


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
mat_observed_pred=matrix(observed_pred,nrow = I,ncol = J)   # matrix with observed/actual entries and predictions for missing entries

#==========================================================================================
# plot biclusters with predictions, biclusters determined by 4 criteria
#==========================================================================================
# minimum total SSE 
train_pred_clust_1=mat_observed_pred[order(id_row_mat[,which.min(sse_mat[-1])]),order(id_col_mat[which.min(sse_mat[-1]),])]
train_pred_clust_plot_1 <- expand.grid(Rows=rownames(data.frame(train_pred_clust_1)), Columns=colnames(data.frame(train_pred_clust_1)))
vec=as.vector(train_pred_clust_1)
train_pred_clust_plot_1$Z <- vec
pred_clust_plot_1=list(data=train_pred_clust_plot_1,vlines=vline_coords_1,hlines=hline_coords_1)

g.1.pred = ggplot(data=pred_clust_plot_1$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=pred_clust_plot_1$vlines,aes(xintercept = v)) +
  geom_hline(data=pred_clust_plot_1$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_discrete(limits = rev(levels(as.factor(train_pred_clust_plot_1$Rows)))) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Biclusters with predictions for min total SSE") + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.1.pred

# maximum average Rand index
train_pred_clust_2=mat_observed_pred[order(id_row_mat[,which.max(avgRowRI)]),order(id_col_mat[which.max(avgColRI),])]
train_pred_clust_plot_2 <- expand.grid(Rows=rownames(data.frame(train_pred_clust_2)), Columns=colnames(data.frame(train_pred_clust_2)))
vec=as.vector(train_pred_clust_2)
train_pred_clust_plot_2$Z <- vec
pred_clust_plot_2=list(data=train_pred_clust_plot_2,vlines=vline_coords_2,hlines=hline_coords_2)

g.2.pred = ggplot(data=pred_clust_plot_2$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=pred_clust_plot_2$vlines,aes(xintercept = v)) +
  geom_hline(data=pred_clust_plot_2$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_pred_clust_plot_2$Rows)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Biclusters with predictions for max average RI") + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.2.pred

# maximum log of joint posterior
train_pred_clust_3=mat_observed_pred[order(id_row_mat[,which.max(lp)]),order(id_col_mat[which.max(lp),])]
train_pred_clust_plot_3 <- expand.grid(Rows=rownames(data.frame(train_pred_clust_3)), Columns=colnames(data.frame(train_pred_clust_3)))
vec=as.vector(train_pred_clust_3)
train_pred_clust_plot_3$Z <- vec
pred_clust_plot_3=list(data=train_pred_clust_plot_3,vlines=vline_coords_3,hlines=hline_coords_3)

g.3.pred = ggplot(data=pred_clust_plot_3$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=pred_clust_plot_3$vlines,aes(xintercept = v)) +
  geom_hline(data=pred_clust_plot_3$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  scale_y_discrete(limits = rev(levels(as.factor(train_pred_clust_plot_3$Rows)))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Biclusters with predictions for max log posterior") + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.3.pred

# minimum penalized SSE 
train_pred_clust_4=mat_observed_pred[order(id_row_mat[,which.min(pq$pSSE)+1]),order(id_col_mat[which.min(pq$pSSE)+1,])]
train_pred_clust_plot_4 <- expand.grid(Rows=rownames(data.frame(train_pred_clust_4)), Columns=colnames(data.frame(train_pred_clust_4)))
vec=as.vector(train_pred_clust_4)
train_pred_clust_plot_4$Z <- vec
pred_clust_plot_4=list(data=train_pred_clust_plot_4,vlines=vline_coords_4,hlines=hline_coords_4)

g.4.pred = ggplot(data=pred_clust_plot_4$data, aes(x=Columns, y=Rows)) + geom_tile(aes(fill= Z)) +
  geom_vline(data=pred_clust_plot_4$vlines,aes(xintercept = v)) +
  geom_hline(data=pred_clust_plot_4$hlines,aes(yintercept = h)) +
  theme_bw() +
  scale_fill_gradientn(colours = color,na.value = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_discrete(limits = rev(levels(as.factor(train_pred_clust_plot_4$Rows)))) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Biclusters with predictions for min penalized SSE") + rremove("xy.text") + rremove("xylab")
# + xlab("...") + ylab("...")
g.4.pred

#==========================================================================================
# compare Rand indices for partitions based on the 3 criteria
# (to used only when "true" clusterings 'tCs','tDs',and 'tEs' are known)
#==========================================================================================
# # for rows
# RRand(id_row_mat[,which.min(sse_mat[-1])],tCs);RRand(id_row_mat[,which.max(avgRowRI)],tCs)
# RRand(id_row_mat[,which.max(lp[-1])+1],tCs); RRand(id_row_mat[,which.min(pq$pSSE)+1],tCs)
# 
# # for columns
# RRand(id_col_mat[which.min(sse_mat[-1]),],tDs);RRand(id_col_mat[which.max(avgColRI),],tDs)
# RRand(id_col_mat[which.max(lp[-1])+1,],tDs); RRand(id_col_mat[which.min(pq$pSSE)+1,],tDs)
# 
# # for elements
# RRand(tEs,elementRI$id_element_mat[,which.min(sse_mat[-1])]);RRand(tEs,elementRI$id_element_mat[,which.max(avgElementRI)])
# RRand(tEs,elementRI$id_element_mat[,which.max(lp[-1])+1]); RRand(tEs,elementRI$id_element_mat[,which.min(pq$pSSE)+1]) 

#==========================================================================================
save.image(file="BayesianBiclustering_WithMissing.RData")
#==========================================================================================