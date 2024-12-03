# example 1 (T=20,sigma=0.7)
rm(list=ls())
library(fpc)
library(GenSA)
library(MASS)
library(glmnet)
library(fda)
library(splines)
library(fpc)
library(mclust)
library(igraph)
#library(tidyverse)
library(broom)
library(flexmix)
#require(ggplot2)
generate_data <- function(sam){
  set.seed(sam)
  t <- seq(0,2,length.out=m)
  t1 <- c(t[c(t<=1)])
  t2 <- c(t[c(t>1)])
  Bf <- create.fourier.basis(c(0,2),5,2)
  f_bas<-eval.basis(t, Bf, returnMatrix=T)
  Bsp=create.bspline.basis(rangeval=c(0,2),norder = 4,nbasis =d)
  bb=eval.basis(t, Bsp, returnMatrix=T)
  #bb <- bs(t,df = d)
  #bb <- as.matrix(bb)
  #beta_1 <- 1 - 1.5*t
  #beta_2 <- 1 + 1.5*t
  beta_3<-(ifelse(t%in%t1,1,0))*t+(ifelse(t%in% t2,1,0))
  #beta_4<-2*((ifelse(t%in%t1,1,0))*(abs(1-t)))
  #beta_5<-3-2*f_bas[,1]-0.8*f_bas[,2]-f_bas[,3]-1.1*f_bas[,4]-0.5*f_bas[,5]
  beta_6<-2*f_bas[,1]+0.8*f_bas[,2]+f_bas[,3]+1.1*f_bas[,4]+0.5*f_bas[,5]
  #beta_7<-2*(1-abs(t-1))
  beta_8<-2*(abs(1-t))
  BETA <- rbind(t(beta_3),t(beta_6),t(beta_8))
  beta <- matrix(c(rep(beta_3,40),
                   rep(beta_6,30),
                   rep(beta_8,20)),
                 nrow=n,ncol=m,byrow=TRUE)
  noise<-matrix(rnorm(N,mu,sigma),n,m)
  betarel <- beta + noise
  return(list( betarel= betarel,BB=bb,g1=beta,t0=t,BETA=BETA))
}

#objective function####
f_optim <- function(beta_vec){ 
  khm <- 0
  beta_optim<-matrix(beta_vec,n,d)
  rss<-sum((betareln-(beta_optim%*%t(B)))^2)/m  
  for (i in 1:K){
    khm <- khm + (apply((beta_optim-rep(1,n)%o%gamma0[i,])^2,1,sum))^s 
  }
  pen<-lambda1*sum((khm/K)^(1/s))
  fun_value<-rss+pen 
  fun_value
} 

#update center####
gamhat<-function(betaest,gammaest){ 
  dis<-matrix(0,nrow=1,ncol=K)
  wei<-matrix(0,n,K)
  for (l in 1:n){
    for(i in 1:K){
      dis[1,i]<-sum((betaest[l,]-gammaest[i,])^2)
    }
    W<-((1/K)^(1/s))*((sum(dis^s)))^(1/s -1)
    wei[l,]<-t(W*(dis^(s-1)))
  }
  wei[is.nan (wei)] <- 0
  gamma_up<-t(wei)%*%betaest/(apply(wei,2,sum)%o%rep(1,(d))) 
  return(gamma_up)
}


#The first derivative of the objective function####
f_d<-function(beta_vec){
  beta<-matrix(beta_vec,n,(d))
  dis<-matrix(0,nrow=1,ncol=K)
  f_d<-f_d1<-f_d2<-matrix(0,n,(d))
  WW <- matrix(0,nrow=d,ncol=K)
  for (j in 1:n){
    betarel1j<-betareln[j,]
    xj<-B 
    f_d1[j,]<--2/m*(t(xj)%*%(betarel1j-xj%*%beta[j,]))
    for(i in 1:K){
      dis[1,i]<-sum((beta[j,]-gamma0[i,])^2)
    }
    W<-(1/K*(sum(dis^s)))^(1/s-1)
    for (l in 1:K)
    {
      WW[,l] <- W*dis[1,l]^(s-1)*2*(beta[j,]-gamma0[l,])/K
    }
    f_d2[j,]<-apply(WW,1,sum)
    f_d[j,]<-f_d1[j,]+lambda1*f_d2[j,]
  }
  return(as.vector(f_d))
}

#RSS####
function_rss <- function(a2_1){
  M_new <- (dim(betareln))[2]
  data1 <- cbind(a2_1,betareln,br)
  RSS1<-0
  RSS1mse<-0
  l0<-1
  l2<-0
  y_n<-rep(NULL,n)
  data2<-data1[order(data1[,1],decreasing=F),]
  k_v <- sort(unique(a2_1))
  for (i in k_v){
    l1 <- length(which((data2[,1]==i)))
    l11 <- length(which((data2[,1]==i)))
    l1<-l1+l2
    y1<-data2[c(l0:l1),c(2:(M_new+1))]
    br_pre <- data2[c(l0:l1),c((M_new+2):(2*M_new+1))]
    v_v<-is.vector(y1)
    y1_v<-as.vector(t(y1))
    br_pre1 <- as.vector(t(br_pre))
    X_RE<-(matrix(rep(t(B),l11),l11*M_new,(d),byrow=TRUE))
    RSS1 <- RSS1+sum(((summary(lm(y1_v~X_RE-1)))$residual)^2)
    beta_pre <- (lm(y1_v~X_RE-1))$coef
    RSS1mse <- RSS1mse+sum((X_RE%*%beta_pre-br_pre1)^2)
    if(v_v==FALSE){
      l0<-l0+dim(y1)[1]
    }else{
      l0<-l0+1
    }
    l2<-l1
  }
  return(list(RSS1=RSS1,RSS1mse=RSS1mse))
}

#main####
n <- 90;m <- 20;n0 <- 100;N <- n*m
K_vector <- c(2:5)
L1 <- length(K_vector)
s=-1;d_vec <- c(4:7)
lambda_vector <- c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,1)
lambda_vector_vv <- c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,1,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,1,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,1,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,1)
epis<-10^(-15);emax=200
sigma <- 0.7;mu=0
a<-c(rep(1,40),rep(2,30),rep(3,20))
L2 <- length(lambda_vector)
RI<-matrix(0,n0,2)
MISE00<-matrix(0,n0,1)
MISE<-matrix(0,n0,1)
MISE0<-matrix(0,n0,1)
lambda_number <- matrix(0,nrow=n0,ncol=1)#Store the results of lambda selection
t1 <- proc.time()
K_F1 <- matrix(0,n0,1)
d_number <- matrix(0,nrow=n0,ncol=1)

for (sam in 1:n0){
  BIC_FINAL <- rep(NULL,L1)
  RI_FINAL <- rep(NULL,L1)
  MISE1_FINAL <- rep(NULL,L1)
  MISE2_FINAL <- rep(NULL,L1)
  MISE_oracle_FINAL <- rep(NULL,L1)#oracle
  for(d in d_vec){
    ss1 <- 0
    data<-generate_data(sam)
    B<-data$BB
    br <- data$g1
    betareln <- data$betarel
    beta0<-matrix(Inf,n,(d))
    for(j in 1:n)
    {
      yj<-betareln[j,]
      xj <- B
      beta0[j,]<-lm(yj~xj-1)$coef
    }
    beta00 <- beta0
    RRI1 <- matrix(Inf,nrow=L1,ncol=L2 + 1)
    MISE1 <- matrix(Inf,nrow=L1,ncol=L2,dimnames=list(c("2","3","4","5")
                                                      ,c("0.0001","0.0005","0.001","0.005","0.01","0.05","0.1","1")))#MSEs
    MISE2 <- matrix(Inf,nrow=L1,ncol=L2,dimnames=list(c("2","3","4","5")
                                                      ,c("0.0001","0.0005","0.001","0.005","0.01","0.05","0.1","1")))#MSEs
    MISE_oracle<- matrix(Inf,nrow=L1,ncol=L2,dimnames=list(c("2","3","4","5")
                                                           ,c("0.0001","0.0005","0.001","0.005","0.01","0.05","0.1","1")))#MSEs
    BIC_value1 <- matrix(Inf,nrow=L1,ncol=L2,dimnames=list(c("2","3","4","5")
                                                           ,c("0.0001","0.0005","0.001","0.005","0.01","0.05","0.1","1")))#MSEs
    
    for(K in K_vector){
      beta0<-beta00
      ss1 <-ss1+1 
      ss2 <- 0
      clu<-kmeans(beta0,K,nstart = 50)
      gamma00 <-clu$center
      a1<-clu$cluster
      RRI1[ss1,1]<-compare(a,a1,method="rand")
      #Lambda
      for (lambda1 in lambda_vector){
        gamma0 <- gamma00
        beta0<-beta00
        ss2 <- ss2 + 1 #lambda+1
        e1<-1;q1=1
        while(e1>epis&q1<=emax){ 
          beta_vec0<-as.vector(beta0)
          outcome <- optim(par=beta_vec0,fn=f_optim,gr=f_d,method = "BFGS",control = list(maxit=100))
          beta_vec1<-outcome$par####beta
          beta1<-matrix(beta_vec1,n,(d))
          q1<-q1+1
          #gamma
          E1 <- try({
            gamma1<-gamhat(beta1,gamma0)
            e1<-max((beta1-beta0)^2)+max((gamma1-gamma0)^2)
            #update center
            gamma0<-gamma1
            beta0<-beta1},silent = FALSE)
          if ("try-error"%in%class(E1)==TRUE){
            break}
        }
        if ("try-error"%in%class(E1)==TRUE){
          next}
        a2<-kmeans(beta1,K,nstart = 50)$cluster
        RRI1[ss1,(ss2 + 1)]<-compare(a2,a,method="rand") 
        RSS_MSE <- function_rss(a2)
        RSS <- RSS_MSE$RSS1
        #BIC
        BIC_value1[ss1,ss2]<-log(RSS/N)+1.3*(d)*((log(m*n))/(m*n))*(K) 
        ####MSEs
        if (K==3){
          MISE1[ss1,ss2] <- ((function_rss(a1))$RSS1mse)/(n*m)
          MISE2[ss1,ss2] <- RSS_MSE$RSS1mse/(n*m)
          MISE_oracle[ss1,ss2]   <-  ((function_rss(a))$RSS1mse)/(n*m)
        }
      } #lambda end of loop
    } #K end of loop
    BIC_FINAL <- cbind(BIC_FINAL,BIC_value1)
    RI_FINAL <- cbind(RI_FINAL,RRI1)
    MISE1_FINAL <- cbind(MISE1_FINAL, MISE1)
    MISE2_FINAL <- cbind(MISE2_FINAL, MISE2)
    MISE_oracle_FINAL <- cbind(MISE_oracle_FINAL, MISE_oracle)
  }
  
  #Selecting K and LAMBDA with BIC
  RI_N1 <- ((which(BIC_FINAL== min(BIC_FINAL), arr.ind = TRUE))[1,1])####Find the location of the minimum value in the BIC
  RI_N2 <- ((which(BIC_FINAL== min(BIC_FINAL), arr.ind = TRUE))[1,2])####Find the location of the minimum value in the BIC
  d_number[sam,1] <- d_vec[ceiling(RI_N2/L2)] 
  K_F1[sam,1]<- RI_N1 + 1####Different seeds were used to select the optimal K
  rinumber1 <- ((ceiling((RI_N2)/(L2)))-1)*(L2+1)+ 1
  rinumber2 <- ((ceiling((RI_N2)/(L2))))+ RI_N2
  RI[sam,1] <- RI_FINAL[RI_N1,rinumber1]#####RI corresponding to minimum BIC
  RI[sam,2] <- RI_FINAL[RI_N1,rinumber2 ]#####RI corresponding to minimum BIC
  lambda_number[sam,1] <- lambda_vector_vv[RI_N2]####LAMBDA selected through BIC
  if( RI_N1[1]==2){#####MSEs
    MISE00[sam,] <- MISE_oracle_FINAL[RI_N1,RI_N2]
    MISE0[sam,1] <-  MISE1_FINAL[RI_N1,RI_N2]
    MISE[sam,1] <-  MISE2_FINAL[RI_N1,RI_N2]
  }
  print(sam)  
}

t2 <- proc.time()
t <- t2 - t1
table(K_F1)
table(d_number)
apply(lambda_number,2,sum)
mean(MISE00) #oracle RMSE
mean(MISE0) #kmeans RMSE
mean(MISE) #khm-va RMSE
(apply(RI,2,sum))/(n0) #RI for kmeans and khm-va
print(paste0('时间',t[3][[1]],'???'))


