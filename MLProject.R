#################################################################
#                     Math 9880 Final Project                   #
#                           R  code                             #
#      Group: Jiajing Niu, Xiyan Tan, Feng Gao, Jingjing Wang   #
#################################################################

############################################################
# Some functions used in the project
############################################################
# 1. Lasso Function
#############################################################
Lasso_fun = function(X, Y, lambda, maxIt = 1e4, tol = 1e-6){
  ##########################################################
  # Input: X = predictor matrix
  #        Y = Response vector
  #        lambda = tuning parameter
  #        maxIt = maximum iteration
  #        tol = convergence tolerance level
  # Output: b = penalized coefficient estimate
  ##########################################################
  # Initialization
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(p > n){
    b = rep(0,p)
  }else{
    b = solve(t(X)%*%X)%*%t(X)%*%Y #From the OLS estimate, if p <= n
  }
  
  b_old = b;
  nIter = 0;
  
  # Precompute X'X and X'Y
  XTX = t(X)%*%X;
  XTY = t(X)%*%Y;
  
  # Shooting loop
  while (nIter < maxIt){
    nIter = nIter+1
    
    for (j in 1:p){
      S0 = XTX[j,-j]%*%b[-j] - XTY[j];  
      if (S0 > lambda){
        b[j] = (lambda-S0) / norm(X[,j],"2")^2
      }else if(S0 < -lambda){
        b[j] = -(lambda+S0) / norm(X[,j],"2")^2
      }else{
        b[j] = 0
      }
    }
    delta = norm(b-b_old,"1");    # Norm change during successive iterations
    
    if (delta < tol){
      break
    }
    b_old = b;
  }
  
  if (nIter == maxIt){
    print("Maximum number of iteration reached, shooting may not converge.")
  }
  return (b)
}

#########################################################################
# 2. GLasso function: Lasso Graphical Model Algorithm
#########################################################################
glasso_fun = function(S, lambda, maxIt = 1e2, tol = 1e-6){
  #############################################################
  ### Graphical lasso algoritm
  ### Ref: Friedman et al. (2007) Sparse inverse covariance estimation with the
  ### graphical lasso. Biostatistics.
  ### Note: This function needs to call an algorithm that solves the Lasso
  # Input:
  #  S -- sample covariance matrix (precision matrix)
  #  lambda --  regularization parameter
  #  maxIt -- maximum number of iterations
  #  tol -- convergence tolerance level
  #
  # Output:
  # Theta -- inverse covariance matrix estimate
  # W -- regularized covariance matrix estimate, W = Theta^-1
  ###############################################################
  
  p = dim(S)[1]
  
  ### 1. Initialization W = S + lambda*I
  W = S + lambda * diag(p) ## diagonal of W remains unchanges
  
  W_old = W
  niter = 0
  
  ## Graphical lasso 
  while(niter < maxIt){
    niter = niter + 1
    prog = sprintf("-------%s th itheration------", niter)
    print(prog)
    ## 2: Repeat for j=1,2,...,p... until convergence
    for (j in 1:p){
      ### a. partition matrix W into part 1
      W11 = W[-j, -j]
      w22 = W[j, j]
      w12 = W[-j, j]
      
      S11 = S[-j, -j]
      s22 = S[j, j]
      s12 = S[-j, j]
      
      ### b. Solve the estimating equations
      d = eigen(W11, symmetric = T)$values
      V = eigen(W11, symmetric = T)$vectors
      
      X = V %*% diag(sqrt(d)) %*% t(V) # W11^(1/2)
      Y = V %*% solve(diag(sqrt(d))) %*% t(V) %*% s12  # W11^(-1/2) * s_12
      
      library(glmnet)
      library(mvtnorm)
      # Way 1. using glmnet lasso for the modifed lasso
      # fit.lasso = glmnet(X, Y, family="gaussian", alpha=1, lambda = lambda)
      # b = as.vector(fit.lasso$beta)
      # Way 2. using Lasso_fun
      b = Lasso_fun(X, Y, lambda = lambda)
      w12 = W11%*%b
      
      ### c. Update W
      W[-j, j] = w12
      W[j, -j] = t(w12)
    }
    if(norm(W-W_old, type = "1")<tol){
      break
    }
    print(norm(W-W_old, type = "1"))
    W_old = W
  }
  if (niter == maxIt){
    print("Maximum number of iteration reached, glasso may not converge.")
  } 
  
  Theta = solve(W)
  
  results  = list()
  results$Theta = Theta
  results$W = W
  return (results)
}


#####################################################################
#  Two Simulations
#####################################################################
# Install the following packages
library(matrixcalc) #pd
library(ggcorrplot)  #plot correlation
library(DensParcorr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
library(lattice)
library(glmnet)
library(mvtnorm)
library(pracma)
library(glasso) 
library(tictoc) # measure running time

########################################################
## Simulation 1:graph estimation for samll p: p=6, N=200
########################################################
# model 1: tridiagonal matrix as follows,
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]  1.0  0.4  0.4  0.0  0.0  0.0
# [2,]  0.4  1.0  0.4  0.4  0.0  0.0
# [3,]  0.4  0.4  1.0  0.4  0.4  0.0
# [4,]  0.0  0.4  0.4  1.0  0.4  0.4
# [5,]  0.0  0.0  0.4  0.4  1.0  0.4
# [6,]  0.0  0.0  0.0  0.4  0.4  1.0
p = 6 #number of nodes
Net.1 = diag(p)
Net.1[abs(row(Net.1) - col(Net.1)) == 1] = 0.4
Net.1[abs(row(Net.1) - col(Net.1)) == 2] = 0.4
Net.1

###########################################
# model 2: precision matrix as follows,
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    2    1    1    0    0    0
# [2,]    1    2    1    1    1    0
# [3,]    1    1    2    1    1    0
# [4,]    0    1    1    2    1    0
# [5,]    0    1    1    1    2    0
# [6,]    0    0    0    0    0    2
p = 6
delta_tilda = 2
mat = matrix(1, nrow = p, ncol = p)
x = c(1,1,1,0,1,1,0,1,1,1,0,0,0,0,0)
mat[upper.tri(mat, diag = F)] = x
mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
Net.2 = mat
diag(Net.2) = delta_tilda
Net.2

##True precision matrix
Net.true = Net.1
#Net.true = Net.2

##check positive definite
is.positive.definite(Net.true)

##generate data from this precision matrix
Sigma.true = solve(Net.true)
N = 200 ##number of replications
Y = mvrnorm(n = N, mu = rep(0, p),  Sigma = Sigma.true)
##add noise
Y = Y + mvrnorm(N, mu = rep(0, p), 0.1*diag(p))
# sample correlation matrix
S = cov(Y)

###################################################
# estimted matrice for different lambda 0.1,0.3,0.6
tic("p=6 model 1, 0.1")
res1 = glasso_fun(S, lambda = 0.1)
Theta_hat_1 = res1$Theta
toc()
# p=6 model 1, 0.1: 0.033 sec elapsed
# p=6 model 2, 0.1: 0.048 sec elapsed

tic("p=6 model 1, 0.3")
res2 = glasso_fun(S, lambda = 0.3)
Theta_hat_2 = res2$Theta
toc()
# p=6 model 1, 0.2: 0.024 sec elapsed
# p=6 model 2, 0.2: 0.027 sec elapsed

tic("p=6 model 1, 0.6")
res3 = glasso_fun(S, lambda = 0.6)
Theta_hat_3 = res3$Theta
toc()
# p=6 model 1, 0.3: 0.007 sec elapsed
# p=6 model 2, 0.3: 0.014 sec elapsed

#######################################################################
## we can use glasso package check our results, For example: lambda=0.1
tic("p=6 model 1, 0.1 glasso")
res = glasso(S, nobs = N, rho = 0.1) # rho=lambda
Theta_hat = res$wi
toc()
# p=6 model 1, 0.1 glasso: 0.011 sec elapsed
#######################################################################

#######################################
# Plot the correlation matrices
#######################################
# plot the true precision matrix
mat = abs(DensParcorr:::prec2part(Net.true))
mat.true = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE, title = "True precision matrix")+
  ggplot2:::scale_x_discrete(limits = nrow(mat):1)+
  ggplot2:::scale_y_continuous(trans = "reverse", breaks = unique(1:nrow(mat)))

# plot the estimated matrices
mat = abs(DensParcorr:::prec2part(Theta_hat_1))
mat.hat.1 = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+labs(title = expression(paste(lambda,'=0.1')))+
  ggplot2:::scale_x_discrete(limits = nrow(mat):1)+
  ggplot2:::scale_y_continuous(trans = "reverse", breaks = unique(1:nrow(mat)))


mat = abs(DensParcorr:::prec2part(Theta_hat_2))
mat.hat.2 = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+labs(title = expression(paste(lambda,'=0.3')))+
  ggplot2:::scale_x_discrete(limits = nrow(mat):1)+
  ggplot2:::scale_y_continuous(trans = "reverse", breaks = unique(1:nrow(mat)))


mat = abs(DensParcorr:::prec2part(Theta_hat_3))
mat.hat.3 = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+labs(title = expression(paste(lambda,'=0.6')))+
  ggplot2:::scale_x_discrete(limits = nrow(mat):1)+
  ggplot2:::scale_y_continuous(trans = "reverse", breaks = unique(1:nrow(mat)))

#############################
# plot the undirected graphs
#############################
# plot the true network
net = network(Net.true, directed = FALSE)
g.true = ggnet2(net, node.size = p,label = TRUE, edge.size = 1, node.color = "tomato",
                edge.color = "tomato",mode = "circle")+labs(title="True graph")

# plot the estimated networks
Net.hat.1 = ifelse(Theta_hat_1<1e-6, 0,1)
net = network(Net.hat.1, directed = FALSE)
g.hat.1 = ggnet2(net, node.size = p,label = TRUE, edge.size = 1, node.color = "grey",
                 edge.color = "grey",mode = "circle")+labs(title=expression(paste(lambda,'=0.1')))

Net.hat.2 = ifelse(Theta_hat_2<1e-6, 0,1)
net = network(Net.hat.2, directed = FALSE)
g.hat.2 = ggnet2(net, node.size = p,label = TRUE, edge.size = 1, node.color = "grey",
                 edge.color = "grey",mode = "circle")+labs(title=expression(paste(lambda,'=0.3')))

Net.hat.3 = ifelse(Theta_hat_3<1e-6, 0,1)
net = network(Net.hat.3, directed = FALSE)
g.hat.3 = ggnet2(net, node.size = p,label = TRUE, edge.size = 1, node.color = "grey",
                 edge.color = "grey",mode = "circle")+labs(title=expression(paste(lambda,'=0.6')))

# Report the plots of true and estimated correlation matrices and undirected graphs
grid.arrange(mat.true,mat.hat.1, mat.hat.2, mat.hat.3,g.true,g.hat.1,g.hat.2,g.hat.3,nrow = 2, ncol = 4)


#####################################################################
# Simulation 2: graph estimation for large p: p=50,100,150, N=100
#####################################################################

##################################
## Three Models for Simulation 2
##################################
ModelGen <- function(p, Net.num){
  if(Net.num == 1){
    #model 1: band matrix
    Net.1 = diag(p)
    Net.1[abs(row(Net.1) - col(Net.1)) == 1] = 0.4
    Net.1[abs(row(Net.1) - col(Net.1)) == 2] = 0.4
    
    Net = Net.1
  }else if(Net.num == 2){
    #model 2: band matrix with some isolated nodes
    #p = 14 #p needs to be >10
    Net.2 = diag(p)
    
    subnet = diag(10)
    subnet[abs(row(subnet) - col(subnet)) == 1] = 0.4
    subnet[abs(row(subnet) - col(subnet)) == 2] = 0.4
    
    for (i in seq(1, floor(p/10), by = 2)){
      block = seq((i-1)*10+1, (i-1)*10+10)
      Net.2[block, block] = subnet
    }
    
    Net = Net.2
  }else if (Net.num == 3){
    #model 3: random sparse matrix 
    Net.3 = rbinom(p*p, 1, 0.1)*0.5
    dim(Net.3) = c(p,p)
    Net.3[lower.tri(Net.3)] = t(Net.3)[lower.tri(Net.3)] #symmetric matrix
    delta_tilda = 4 ##choose to garantee theta.true is positive definite
    diag(Net.3) = delta_tilda
    
    is.positive.definite(Net.3)
    Net = Net.3
  }
  return(Net)
}
  
# Change p=50,100,150 and model=1,2,3 to generate different structures
p_seq = c(50,100,150)
Model_seq = c(1,2,3)


for (Net.num in Model_seq){
  for (p in p_seq){
    Net = ModelGen(p=p, Net.num=Net.num)
    
    # lambda sequence
    lambda_sq = seq(0.05, 1 ,by = 0.02)
    Net.true = Net
    
    Theta_hat = rep(NA, length(lambda_sq)*p*p)
    dim(Theta_hat) = c(length(lambda_sq),p,p)
    
    ##generate data from this precision matrix
    for (i in 1:length(lambda_sq)){
      Sigma.true = solve(Net.true)
      N = 100 ##number of replications
      Y = mvrnorm(n = N, mu = rep(0, p),  Sigma = Sigma.true)
      
      ##add noise
      Y = Y + mvrnorm(N, mu = rep(0, p), 0.5*diag(p))
      S = cov(Y)
      
       res = glasso_fun(S, lambda = lambda_sq[i])
       Theta_hat[i,,] = res$Theta
      # check using glasso package
      # res = glasso(S, nobs = N, rho = lambda_sq[i])
      # Theta_hat[i,,] = res$wi
    }
    Results = list()
    Results$Theta_hat = Theta_hat
    Results$Theta.true = Net.true
    save(Results,file = sprintf("Results_Model%s_p_%s.RData",Net.num,p))
  }
}
# p=50, model 1, 0.1: 29.439 sec elapsed
# p=100, model 1, 0.1: 36.005 sec elapsed
# p=150, model 1, 0.1: 124.544 sec elapsed
# package total: 10.171 sec elapsed

#########################################################################
##  ROC plot and AUC comparison
#########################################################################
for (Net.num in c(1,2,3)){
  lambda_sq = seq(0.05, 1,by = 0.02)
  p_seq = c(50,100,150)
  res.data = c()
  AUC = rep(-99,length(p_seq))
  for (i in 1:length(p_seq)){
    p = p_seq[i]
    load(file = sprintf("Results_Model%s_p_%s.RData",Net.num,p))
    res.Theta = Results$Theta_hat
    Net.true = Results$Theta.true
    
    ##plot the ROC 
    pred.mat = c()
    for (k in 1:length(lambda_sq)){
      Network = res.Theta[k,,]
      pred.mat = cbind(pred.mat, ifelse(as.vector(Network)<1e-6,0,1))
    }
    
    TPR = rep(-99, length(lambda_sq))
    FPR = rep(-99, length(lambda_sq))
    
    Labels = as.vector(Net.true)
    for (s in 1:length(lambda_sq)){
      Predictions = pred.mat[,s]
      TP = sum(Labels > 0 & Predictions == 1)
      FP = sum(Labels == 0 & Predictions == 1)
      FN = sum(Labels > 0 & Predictions == 0)
      TN = sum(Labels == 0 & Predictions == 0)
      TPR[s] = TP/(TP+FN)
      FPR[s] = FP/(FP+TN)
    }
    
    TPR[length(lambda_sq)+1] = 0
    FPR[length(lambda_sq)+1] = 0
    TPR = append(1, TPR)
    FPR = append(1, FPR)
    p_num = rep(p,(length(lambda_sq)+2))
    temp = cbind(FPR,TPR,p_num)
    res.data = data.frame(rbind(res.data,temp))
    # calculate area under curve (AUC)
    AUC[i] = abs(trapz(FPR, TPR))
  }
  print(AUC)
  res.data$p_num = factor(res.data$p_num)
  p = ggplot(res.data, aes(x = FPR, y = TPR, colour = p_num)) + 
    geom_line(size = 1) + 
    ylab(label="TPR") + 
    xlab("FPR") + 
    scale_colour_manual(values=c("red","blue","yellow")) +
    labs(title=sprintf("Model %s: p=50,100,150", Net.num))
  nam = paste("p", Net.num,sep = "")
  assign(nam, p)
}

grid.arrange(p1, p2, p3, nrow = 1)
# AUC output
# [1] 0.9715476 0.9672992 0.9663809
# [1] 0.9756400 0.9710829 0.9691902
# [1] 0.6741142 0.6346550 0.6588740


#####################################################################
#  One Application: EEG Data
#####################################################################
# Full dataset ------ Palmetto 9.69G
library(rgl)
library(eegkit)
library(eegkitdata)
# Note: since the data are very big, we extracted them on Palmetto 
# #(1)# download the eeg_full data from UCI Machine Learning Repository
#  https://archive.ics.uci.edu/ml/datasets/eeg+database
# #(2)# extract condition "S1" and save as .rda
eegS1 = geteegdata(indir="/home/xtan/9880ML/project/eeg_full/", cond="S1",filename="eegfullS1")
# #(3)# extract condition "S2m" and save as .rda
eegS2m = geteegdata(indir="/home/xtan/9880ML/project/eeg_full/",cond="S2m",filename="eegfullS2m")
# #(4)# extract condition "S2n" and save as .rda
eegS2n = geteegdata(indir="/home/xtan/9880ML/project/eeg_full/",cond="S2n",filename="eegfullS2n")
# #(5)# combine conditions
eegdata = rbind(eegS1,eegS2m,eegS2n)
# types of each variable
levels(eegdata$group)  #"a" and "c"
levels(eegdata$channel) #64
summary(eegdata$time)  #256: 0-255
levels(eegdata$subject) #122

########################################################
#  We consider the average of all trials 
# for each subject under the single stimulus conditions
########################################################
# Data :eegS1
# Two groups of dataset: alcoholic and control
alcoholic = eegS1[eegS1$group=="a",]
control = eegS1[eegS1$group=="c",]

# eegS1
levels(eegS1$group)  #2
levels(eegS1$channel) #64
summary(eegS1$time)  #256: 0-255
levels(eegS1$subject) #122

# alcoholic 
summary(alcoholic$group)  #a
levels(alcoholic$channel) #64
summary(alcoholic$time)  #256: 0-255
summary(alcoholic$subject) #77

# control
summary(control$group)  #c
levels(control$channel) #64
summary(control$time)  #256: 0-255
summary(control$subject) #45

# plot EEG signals: using 6 examples/positions
par(mfrow=c(3,2))
for (i in c("AF1","C5","FP1","F8","PZ","O1")){
  data("eegdata")
  # get "PZ" electrode from "eegdata" data
  idx = which(eegdata$channel==i)
  eegdata = eegdata[idx,]
  # get average and standard error (note se=sd/sqrt(n))
  eegmean = tapply(eegdata$voltage,list(eegdata$time,eegdata$group),mean)
  eegse = tapply(eegdata$voltage,list(eegdata$time,eegdata$group),sd)/sqrt(50)
  # plot results with legend
  tseq = seq(0,1000,length.out=256)
  eegtime(tseq,eegmean[,2],voltageSE=eegse[,2],ylim=c(-10,10),main=i)
  eegtime(tseq,eegmean[,1],vlty=1,vcol="red",voltageSE=eegse[,1],scol="pink",add=TRUE)
  legend("bottomright",c("controls","alcoholics"),lty=c(1,1),
         lwd=c(2,2),col=c("blue","red"),bty="n")
}

########################################################################################
# average of voltage
alcoholmean = tapply(alcoholic$voltage,list(alcoholic$subject,alcoholic$channel),mean)
alcoholmean = na.omit(alcoholmean)

controlmean = tapply(control$voltage,list(control$subject,control$channel),mean)
controlmean = na.omit(controlmean)

# median of voltage (similar as mean, not in the report)
alcoholmedian = tapply(alcoholic$voltage,list(alcoholic$subject,alcoholic$channel),median)
alcoholmedian = na.omit(alcoholmedian)

controlmedian = tapply(control$voltage,list(control$subject,control$channel),median)
controlmedian = na.omit(controlmedian)

# trapz integration of voltage
temp.a = rep(-99,nrow(alcoholic)/256)
for(i in (0:(nrow(alcoholic)/256-1))){
  temp.a[i+1] = trapz(alcoholic$time[(i*256+1):(i*256+256)],alcoholic$voltage[(i*256+1):(i*256+256)])
}

data_a = alcoholic[alcoholic$time==0,]

data_a$int = temp.a

alcoholint = tapply(data_a$int,list(data_a$subject,data_a$channel),mean)
alcoholint = na.omit(alcoholint)

temp.c = rep(-99,nrow(control)/256)
for(i in (0:(nrow(control)/256-1))){
  temp.c[i+1] = trapz(control$time[(i*256+1):(i*256+256)],control$voltage[(i*256+1):(i*256+256)])
}

data_c = control[control$time==0,]

data_c$int = temp.c
controlint = tapply(data_c$int,list(data_c$subject,data_c$channel),mean)
controlint = na.omit(controlint)

# Standardize the data
Standard = function(X){
  xb = apply(X,2,mean)
  sx = apply(X,2,sd)
  Xs = t((t(X)-xb)/sx)
  return(Xs)
}
alcoholmean = Standard(alcoholmean)
controlmean = Standard(controlmean)
alcoholint = Standard(alcoholint)
controlint = Standard(controlint)

#######################
# Estimated matrix
#######################
# mean case
S1 = cov(alcoholmean)
S2 = cov(controlmean)
####################
# integration case
S1 = cov(alcoholint)
S2 = cov(controlint)

# plot correlation matrix
# example: alcoholic lambda = 0.01
tic("EEG: p=64")
res = glasso_fun(S1, lambda = 0.01)
Theta_hat = res$theta
toc()
# EEG: p=64: 760.237 sec elapsed
mat = abs(DensParcorr:::prec2part(Theta_hat))
mat.hat = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+
          ggplot2:::scale_y_continuous(trans = "reverse", breaks = c(0,20,40,60))


#########################################
# glasso in EEG
#########################################
lambda_seq = seq(0.05,1,by=0.02)

# total number of edges
choose(64,2)
# we would like to keep approximate 2.5% connected edges
choose(64,2)*0.025 #50.4

##########################
# Alcoholic 
count = c()
for (lambda in lambda_seq){
   res = glasso_fun(S1, lambda = lambda)
   Theta_hat = res$Theta
  # check using package
  # res = glasso(S1, nobs = 77, rho = lambda)
  # Theta_hat = res$wi
  Net.hat = ifelse(Theta_hat<1e-6, 0,1)
  count = cbind(count, sum(as.vector(Net.hat)))
}
# Estimate Theta for alcoholic
lambda_opt_a = lambda_seq[min(which(count < 164))]
lambda_opt_a
res_a = glasso_fun(S1, lambda = lambda_opt_a)
Theta_hat_a = res_a$Theta

###########################
# Control
count = c()
for (lambda in lambda_seq){
   res = glasso_fun(S2, lambda = lambda)
   Theta_hat = res$Theta
  # check using package
  # res = glasso(S2, nobs = 45, rho = lambda)
  # Theta_hat = res$wi
  Net.hat = ifelse(Theta_hat<1e-6, 0,1)
  count = cbind(count, sum(as.vector(Net.hat)))
}
# Estimate Theta for control
lambda_opt_c = lambda_seq[min(which(count < 164))]
lambda_opt_c
res_c = glasso_fun(S2, rho = lambda_opt_c)
Theta_hat_c = res_c$Theta

# plot estimated matrice
mat = abs(DensParcorr:::prec2part(Theta_hat_a))
mat.hat.a = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+
            labs(title = expression(paste('Alcoholic, ',lambda,'=0.19')))+
            ggplot2:::scale_y_continuous(trans = "reverse", breaks = c(0,20,40,60))
   

mat = abs(DensParcorr:::prec2part(Theta_hat_c))
mat.hat.c = ggcorrplot(mat[1:nrow(mat),], show.legend = FALSE)+
            labs(title = expression(paste('Control, ',lambda,'=0.27')))+
            ggplot2:::scale_y_continuous(trans = "reverse", breaks = c(0,20,40,60)) 

grid.arrange(mat.hat.a, mat.hat.c, nrow = 1)

# plot estimated graphs
Net.hat.a = ifelse(Theta_hat_a<1e-6, 0,1)
Net.hat.c = ifelse(Theta_hat_c<1e-6, 0,1)
rownames(Net.hat.a) = rownames(S1)
colnames(Net.hat.a) = colnames(S1)
rownames(Net.hat.c) = rownames(S2)
colnames(Net.hat.c) = colnames(S2)


# igraph plot : Fix coordinate
library(igraph)
# plot EEG coordinates
par(mfrow=c(1,1))
data(eegcoord)
# plot custom subset of electrodes
myelectrodes <- c(rownames(Net.hat.a),"IZ","A1","A2")
eegcap(myelectrodes)
eegcoord$names = rownames(eegcoord)
eegcoord$names[1] = "X"
eegcoord$names[2] = "Y"
eegcoord$xproj[1] = -11
eegcoord$xproj[2] = 11
eegcoord$names[eegcoord$names=="IZ"] = "nd"
eegcoord = subset(eegcoord, names %in% rownames(Net.hat.a))
rownames(eegcoord) = eegcoord$names
enames = rownames(eegcoord)
plot(eegcoord[,4],eegcoord[,5],cex=2,col="green",pch=19)
text(eegcoord[,4],eegcoord[,5],labels=enames,col="blue",cex=0.6)

# define the corrdinate
eegcoord <- eegcoord[rownames(Net.hat.a),]
coords <- eegcoord[,4:5]
coords <- as.matrix(coords)
colnames(coords) <- NULL
rownames(coords) <- NULL
l = coords

#par(mfrow=c(1,2))
# 1. alcoholic
Net.a = Net.hat.a
diag(Net.a) = 0
net.a = graph_from_adjacency_matrix(Net.a,mode="Undirected")
plot(net.a, layout=l,vertex.size=10, vertex.label.cex=0.7,vertex.color="orange",edge.curved=.1,
     vertex.frame.color="#ffffff",edge.color="orange",vertex.label.color="black")


# 2. control
Net.c = Net.hat.c
diag(Net.c) = 0
net.c = graph_from_adjacency_matrix(Net.c,mode="Undirected")
plot(net.c, layout=l,vertex.size=10, vertex.label.cex=0.7,vertex.color="orange",edge.curved=.1,
     vertex.frame.color="#ffffff",edge.color="orange",vertex.label.color="black")

# edges comparison
# > E(net.a)
# + 41/41 edges from 6786922 (vertex names):
#   [1] AF1--P2  AF1--P4  AF1--P5  AF2--P2  AF2--P3  AF2--P4  AF2--POZ AF7--O2  AF7--P2  AF7--PO1 AF7--PO8
# [12] AF8--CZ  AFZ--P2  AFZ--P3  AFZ--P5  AFZ--PO1 CZ --FPZ CZ --X   F1 --PO8 F2 --POZ F3 --O2  F3 --OZ 
# [23] FP1--P4  FP1--PO7 FP1--TP8 FP2--P1  FP2--P2  FP2--P3  FP2--P4  FP2--PO7 FZ --O2  FZ --P3  FZ --P4 
# [34] FZ --P5  FZ --PO7 FZ --PO8 O2 --X   P4 --X   P5 --X   P7 --X   PO7--X  

# > E(net.c)
# + 45/45 edges from 98a8fd1 (vertex names):
#   [1] AF1--CP4 AF1--P6  AF1--PO1 AF2--CPZ AF2--P2  AF8--P1  AF8--PO2 AF8--POZ AFZ--CP2 AFZ--CP4 AFZ--P6 
# [12] CP1--F2  CP2--F1  CP2--F5  CP2--FCZ CP2--FZ  CP4--F3  CP4--FC1 CP4--FCZ CP4--FP1 CPZ--F2  CPZ--F4 
# [23] CPZ--FC6 F2 --P2  F4 --P2  F8 --P2  FC2--P2  FC3--POZ FC4--P2  FCZ--P2  FCZ--P4  FCZ--P6  FCZ--PZ 
# [34] FP1--PO1 FP1--PO2 FP1--POZ FP2--P7  FP2--PO2 FPZ--O1  FPZ--PO1 FPZ--PO2 FZ --P2  FZ --P6  PO2--X  
# [45] PZ --Y  
