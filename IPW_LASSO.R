##############################
###--- Loading Libraries ---##
##############################

library(Matrix)
library(MASS)
library(glmnet)
library(lars)

#######################################################################################################################

#(B)

#######################################################
##--- Simulation for Full Observed Data 500 times ---##
#######################################################


Overall.Dataset=vector(mode = "list",length = 500)
for(k in 1:500){
  
  ###--- Correlation Matrix ---###
  x=matrix(0,10,10)
  for(i in 1:10)
  {
    for(j in 1:10)
    {
      if(i<j)
      {
        x[i,j]=0.5^(abs(i-j))
      }
      else if (i==j)
      {
        x[i,j]=1
      }
      else
      {
        x[i,j]=0
      }
    }
  }
  
  corr.matrix=as.matrix(forceSymmetric(x)) # Correlation matrix
  
  ###--- 1st 10  samples form normal dist. ---### 
  Data=MASS::mvrnorm(200,rep(0,10),corr.matrix)
  
  ###--- Simulating 11th and 12th variable from Exp(1/2) and U(0,1) respectively ---##
  AA=matrix(c(rexp(200,1/2),runif(200,0,1)),ncol = 2,byrow = F)
  
  ###--- Overall Covariates ---###
  Data=cbind(Data,AA)
  
  ###--- Random Error ---###
  error=rnorm(200,0,2)
  
  ###--- Model Paraameter ---###
  beta=matrix(c(1,1,0,0,2.5,1,0,1,0,0,0,0),ncol = 1)
  
  ###--- Responce Vector ---###
  y=Data%*%beta+error
  
  ###--- Overall data set ---###
  FULL.Data=cbind(y,Data,error)
  
  ###--- Converting into a data frame ---###
  FULL.Data=as.data.frame(FULL.Data)
  names(FULL.Data)=c("Y",paste(rep("X",12),1:12,sep="_"),"ERROR") ## Setting Names of the variables
  
  ###--- Storing in to a List variable ---##
  Overall.Dataset[[k]]=FULL.Data
  
  ###--- Also Exporting the 500 simulated data sets into a text file ---###
  C=paste("Full Data",k,sep = "_")
  s=paste(C,"text",sep = ".")
  b=paste("D:","Users","User","OneDrive","OneDrive","Documents","M_SC_Project Data",s,sep="/")
  write.table(FULL.Data,b,sep="\t")
}

####################################################################################################################################

#matrix(c(rep(0,3),0.5,0.5,0,rep(0,6)),nrow=12,ncol = 1,byrow = T)

#matrix(c(rep(0,3),0,0.7,0,rep(0,6)),nrow=12,ncol = 1,byrow = T)

#matrix(c(rep(0,4),0.96,rep(0,7)),nrow=12,ncol = 1,byrow = T)

#(C)

####################################################
##--- Writing A Function Relivant For Our Work ---##
####################################################


Our.Method=function(Our.Original.Data, Logistic.Coeficient){
  
  m.row= nrow(Our.Original.Data)
  n.col= ncol(Our.Original.Data)
  
  ##=============================##
  ## Simulation for Missing Data ##
  ##=============================##
  
  ###--- Calling the Full observed dataset into R ---###
  Sim.Data= Our.Original.Data
  
  ###--- Logistic Model Parameter for MNAR ---###
  #psi=matrix(c(rep(0,3),0,0.7,0,rep(0,6)),nrow=12,ncol = 1,byrow = T)
  
  #psi=matrix(c(rep(0,4),0.96,rep(0,7)),nrow=12,ncol = 1,byrow = T)
  psi = Logistic.Coeficient
  
  ###--- Missing Probability ---###
  missing.Probability=function(data){
    prob_V=NULL
    for(i in 1:nrow(data))
    {
      store=(exp(t(psi)%*%t(data)[,i]))/(1+(exp(t(psi)%*%t(data)[,i])))  
      prob_V=c(prob_V,store)
    }
    
    return(prob_V)
  }
  X = as.matrix(Sim.Data[,2:13])  # All Covariates
  p = missing.Probability(X)       # Missing Probability
  
  ###--- Missing Indicator ---###
  Missing.Indecator = rbinom(m.row,1,missing.Probability(X))
  
  ###--- Missing Percentage ---###
  percentage = (1-sum(Missing.Indecator)/length(Missing.Indecator))*100
  
  ###--- Creating New Missing Data set And Exporting ---###
  
  New.X_5 = NULL                    #Deleting values of X_5 from the Full observed Simulated Data 
  for(i in 1:m.row)
  {
    if(Missing.Indecator[i]==1)
    {
      New.X_5[i]=Sim.Data$X_5[i]
    }
    else
    {
      New.X_5[i]=NA
    }
  }
  #New.X_5 is the new X_5 with missing observations
  
  ###--- Converting into a data frame and exporting ---###
  Missing.Data = as.data.frame(cbind(Sim.Data[,1:5],New.X_5,Sim.Data[,7:14]))
  names(Missing.Data) = c("Y",paste(rep("X",12),1:12,sep="_"),"ERROR")
  
  ##====================##
  ## Complete Case Data ##
  ##====================##
  
  ###--- Calling the Missing dataset into R ---###
  Miss.Data=Missing.Data
  
  ##---- Complete Case Data ----##
  Complete.Data=na.omit(Miss.Data)     # Deleting missing Observations
  
  
  ##====================##
  ## Lasso on Full Data ##
  ##====================##
  
  ##--- Preparing Data For LASSO ---##
  X.Full = as.matrix(Sim.Data[,2:13])
  Y.Full = Sim.Data$Y
  
  ###-- LASSO Model --###
  LASSO.FULL=lars(X.Full,Y.Full,type ="lasso")
  
  ###--- Choice of Lambda 10-fold CV ---###
  CV.FULL = lars::cv.lars(X.Full,Y.Full,type = "lasso",K=10,mode = "fraction")
  
  ###--- lambda 1-se ---###
  limit=min(CV.FULL$cv) + CV.FULL$cv.error[which.min(CV.FULL$cv)]   # Lambda Limit
  s.cv=CV.FULL$index[min(which(CV.FULL$cv < limit))]                # Choise of Lambda 
  
  ###--- Lasso Coefficient ---###
  #coef(LASSO.FULL, s = s.cv, mode = "fraction")                     # Lasso coefficients
  
  
  ##========================##
  ## Lasso on Complete Case ##
  ##========================##
  
  ##--- Preparing Data For LASSO ---##
  X.Com = as.matrix(Complete.Data[,2:13])
  Y.Com = Complete.Data$Y
  
  ###-- LASSO Model --###
  LASSO.Com=lars(X.Com,Y.Com,type ="lasso")
  
  ###--- Choice of Lambda 10-fold CV ---###
  CV.Com = lars::cv.lars(X.Com,Y.Com,type = "lasso",K=10,mode = "fraction")
  
  ###--- lambda 1-se ---###
  limit.Com = min(CV.Com$cv) + CV.Com$cv.error[which.min(CV.Com$cv)]   # Lambda Limit
  s.cv.Com = CV.Com$index[min(which(CV.Com$cv < limit.Com))]                # Choise of Lambda 
  
  ###--- Lasso Coefficient ---###
  #coef(LASSO.Com, s = s.cv.Com, mode = "fraction")                     # Lasso coefficients
  
  
  ##====================================##
  ## IPW-Lasso with known Probabilities ##
  ##====================================##
  
  ##======== Transforming Data =========##
  
  ##---- Finding Missing Observations -----##
  SS=seq(1,m.row,1)
  
  MDU=Missing.Indecator*SS
  
  MDU.o=MDU[(MDU!=0)]
  
  MDU.m=setdiff(1:m.row,MDU.o)
  
  ##---- Finding Weight(Missing Prob.) ----##
  pi=p[MDU.o]
  
  ##---- Transforming Design Matrix ----##
  x.New=as.matrix(Complete.Data[,2:13])
  X_trans=matrix(NA,nrow=length(pi),ncol=12)
  Y_trans = NULL
  for(i in 1:length(pi))
  {
    X_trans[i,]=x.New[i,]/sqrt(pi[i])
    Y_trans[i]=Complete.Data$Y[i]/sqrt(pi[i])
  }
  
  ##--- Exporting the Transformed data set in to a text file ---##
  transform.Data=as.data.frame(cbind(Y_trans,X_trans))
  names(transform.Data)=c("Y_trans",paste(rep("X",12),1:12,rep("trans",12),sep="_"))
  
  
  
  ##========== Applying LASSO ===========##
  
  ##--- Preparing Data For LASSO ---##
  X.Trans.Known.Prob = as.matrix(transform.Data[,2:13])
  Y.Trans.Known.Prob = transform.Data$Y
  
  ###-- LASSO Model --###
  IPW_LASSO.Known.Prob = lars(X.Trans.Known.Prob,Y.Trans.Known.Prob,type ="lasso")
  
  ###--- Choice of Lambda 10-fold CV ---###
  CV.IPW_LASSO.Known.Prob = lars::cv.lars(X.Trans.Known.Prob,Y.Trans.Known.Prob,
                                          type = "lasso",K=10,mode = "fraction")
  
  ###--- lambda 1-se ---###
  limit.IPW_LASSO.Known.Prob = min(CV.IPW_LASSO.Known.Prob$cv) + 
    CV.IPW_LASSO.Known.Prob$cv.error[which.min(CV.IPW_LASSO.Known.Prob$cv)]   # Lambda Limit
  s.cv.IPW_LASSO.Known.Prob = CV.IPW_LASSO.Known.Prob$index[min(which(CV.IPW_LASSO.Known.Prob$cv <
                                                                        limit.IPW_LASSO.Known.Prob))] # Choise of Lambda 
  
  ###--- Lasso Coefficient ---###
  #coef(IPW_LASSO.Known.Prob, s = s.cv.IPW_LASSO.Known.Prob, mode = "fraction")                     # Lasso coefficients
  
  
  
  ##========================================##
  ## IPW-Lasso with Estimated Probabilities ##
  ##========================================##
  
  ##========= Estimating Logistic Regression Coefficients =========##
  
  ##--- Preparing Data ---##
  Design    = Missing.Data[,c(2:5,7:13)]
  Dependent = Missing.Indecator
  Data_for_Logistic_fit = as.data.frame(cbind(Dependent,Design))
  
  ##--- Fitting Logistic Regression on Complete case Covariates ---## 
  Logistic.Fit = glm(Dependent~.-1,data=Data_for_Logistic_fit,family = "binomial")
  
  ##--- Calculating Estimated Missing Probability ---##
  psi_hat = matrix(coefficients(Logistic.Fit),nrow=11,ncol = 1,byrow = T)
  
  p_hat = NULL
  for(i in 1:nrow(Design))
  {
    store = (exp(t(psi_hat)%*% t(Design)[,i]))/(1+(exp(t(psi_hat)%*% t(Design)[,i])))  
    p_hat = c(p_hat,store)
  }
  
  
  ##======== Transforming Data =========##
  
  ##---- Finding Missing Observations -----##
  SS=seq(1,m.row,1)
  
  MDU=Missing.Indecator*SS
  
  MDU.o=MDU[(MDU!=0)]
  
  MDU.m=setdiff(1:m.row,MDU.o)
  
  
  ##---- Finding Weight(Missing Prob.) ----##
  pi_hat = p_hat[MDU.o]
  
  ##---- Transforming Design Matrix ----##
  x.New_1 = as.matrix(Complete.Data[,2:13])
  X_trans_1 = matrix(NA,nrow=length(pi_hat),ncol=12)
  Y_trans_1 = NULL
  for(i in 1:length(pi_hat))
  {
    X_trans_1[i,]=x.New[i,]/sqrt(pi_hat[i])
    Y_trans_1[i] = Complete.Data$Y[i]/sqrt(pi_hat[i])
  }
  
  
  ##--- Exporting the Transformed data set in to a text file ---##
  transform.Data_for_estimated_Pi = as.data.frame(cbind(Y_trans_1,X_trans_1))
  dim(transform.Data_for_estimated_Pi)
  names(transform.Data_for_estimated_Pi)=c("Y_trans",paste(rep("X",12),1:12,rep("trans_Estimted_Pi",12),sep="_"))
  
  
  ##========== Applying LASSO ===========##
  
  ##--- Preparing Data For LASSO ---##
  X.Trans.Estimated.Prob = as.matrix(transform.Data_for_estimated_Pi[,2:13])
  Y.Trans.Estimated.Prob = transform.Data_for_estimated_Pi$Y_trans
  
  ###-- LASSO Model --###
  IPW_LASSO.Estimated.Prob = lars(X.Trans.Estimated.Prob,Y.Trans.Estimated.Prob,type ="lasso")
  
  ###--- Choice of Lambda 10-fold CV ---###
  CV.IPW_LASSO.Estimated.Prob = lars::cv.lars(X.Trans.Estimated.Prob,Y.Trans.Estimated.Prob,type = "lasso",K=10,mode = "fraction")
  
  ###--- lambda 1-se ---###
  limit.IPW_LASSO.Estimated.Prob = min(CV.IPW_LASSO.Estimated.Prob$cv) + 
    CV.IPW_LASSO.Estimated.Prob$cv.error[which.min(CV.IPW_LASSO.Estimated.Prob$cv)]   # Lambda Limit
  s.cv.IPW_LASSO.Estimated.Prob = CV.IPW_LASSO.Estimated.Prob$index[min(which(CV.IPW_LASSO.Estimated.Prob$cv <
                                                                                limit.IPW_LASSO.Estimated.Prob))] # Choise of Lambda 
  
  
  ###--- Lasso Coefficient ---###
  #coef(IPW_LASSO.Estimated.Prob, s = s.cv.IPW_LASSO.Estimated.Prob, mode = "fraction")  # Lasso coefficients
  
  
  
  ##===============================================================##
  ## IPW-Lasso with Estimated Probabilities (Using Logistic Lasso) ##
  ##===============================================================##
  
  ##========= Estimating Logistic-Lasso Regression Coefficients =========##
  
  ##--- Preparing Data ---##
  Design    = Missing.Data[,c(2:5,7:13)]
  Dependent = Missing.Indecator
  
  ##--- Fitting Logistic Lasso Regression on Complete case Covariates ---## 
  Logistic_Lasso.Fit       = lars(as.matrix(Design),as.vector(Dependent),type = "lasso")
  CV.Logistic_Lasso.Fit    = lars::cv.lars(as.matrix(Design),as.vector(Dependent),type = "lasso",K=10,mode = "fraction")
  limit.Logistic_Lasso.Fit = min(CV.Logistic_Lasso.Fit$cv) + CV.Logistic_Lasso.Fit$cv.error[which.min(CV.Logistic_Lasso.Fit$cv)]   
  s.cv.Logistic_Lasso.Fit  = CV.Logistic_Lasso.Fit$index[min(which(CV.Logistic_Lasso.Fit$cv < limit.Logistic_Lasso.Fit))]                
  Cof= coef(Logistic_Lasso.Fit, s = s.cv.Logistic_Lasso.Fit, mode = "fraction")                     # Lasso coefficients
  
  ##--- Calculating Estimated Missing Probability ---##
  psi_hat_Lasso = matrix(Cof,nrow=11,ncol = 1,byrow = T)
  
  p_hat_Lasso = NULL
  for(i in 1:nrow(Design))
  {
    store = (exp(t(psi_hat_Lasso)%*% t(Design)[,i]))/(1+(exp(t(psi_hat_Lasso)%*% t(Design)[,i])))  
    p_hat_Lasso = c(p_hat_Lasso,store)
  }
  
  
  ##======== Transforming Data =========##
  
  ##---- Finding Missing Observations -----##
  
  SS=seq(1,m.row,1)
  
  MDU=Missing.Indecator*SS
  
  MDU.o=MDU[(MDU!=0)]
  
  MDU.m=setdiff(1:m.row,MDU.o)
  
  ##---- Finding Weight(Missing Prob.) ----##
  pi_hat_Lasso = p_hat_Lasso[MDU.o]
  
  ##---- Transforming Design Matrix ----##
  x.New_2 = as.matrix(Complete.Data[,2:13])
  X_trans_2 = matrix(NA,nrow=length(pi_hat_Lasso),ncol=12)
  Y_trans_2 = NULL
  for(i in 1:length(pi_hat))
  {
    X_trans_2[i,]=x.New[i,]/sqrt(pi_hat_Lasso[i])
    Y_trans_2[i] = Complete.Data$Y[i]/sqrt(pi_hat_Lasso[i])
  }
  
  ##--- Exporting the Transformed data set ---##
  transform.Data_for_estimated_Pi_LASSO = as.data.frame(cbind(Y_trans_2,X_trans_2))
  names(transform.Data_for_estimated_Pi_LASSO)=c("Y_trans_lasso",paste(rep("X",12),
                                                                       1:12,
                                                                       rep("trans_Estimted_Pi_lasso",12),
                                                                       sep="_"))
  
  
  
  ##========== Applying LASSO ===========##
  
  ##--- Preparing Data For LASSO ---##
  X.Trans.Estimated.Prob_lasso = as.matrix(transform.Data_for_estimated_Pi_LASSO[,2:13])
  Y.Trans.Estimated.Prob_lasso = transform.Data_for_estimated_Pi_LASSO$Y_trans_lasso
  
  ###-- LASSO Model --###
  IPW_LASSO.Estimated.Prob_lasso = lars(X.Trans.Estimated.Prob_lasso,Y.Trans.Estimated.Prob_lasso,type ="lasso")
  
  ###--- Choice of Lambda 10-fold CV ---###
  CV.IPW_LASSO.Estimated.Prob_lasso = lars::cv.lars(X.Trans.Estimated.Prob_lasso,Y.Trans.Estimated.Prob_lasso,type = "lasso",K=10,mode = "fraction")
  
  ###--- lambda 1-se ---###
  limit.IPW_LASSO.Estimated.Prob_lasso = min(CV.IPW_LASSO.Estimated.Prob_lasso$cv) + 
    CV.IPW_LASSO.Estimated.Prob_lasso$cv.error[which.min(CV.IPW_LASSO.Estimated.Prob_lasso$cv)]   # Lambda Limit
  s.cv.IPW_LASSO.Estimated.Prob_lasso = CV.IPW_LASSO.Estimated.Prob_lasso$index[min(which(CV.IPW_LASSO.Estimated.Prob_lasso$cv < 
                                                                                            limit.IPW_LASSO.Estimated.Prob_lasso))] # Choise of Lambda 
  
  
  ###--- Lasso Coefficient ---###
  #coef(IPW_LASSO.Estimated.Prob_lasso, s = s.cv.IPW_LASSO.Estimated.Prob_lasso, mode = "fraction")    # Lasso coefficients
  
  
  ###-- Storing the LASSO Outputs into a table ---###
  Ulti.R.Table = as.data.frame(matrix(c(coef(LASSO.FULL, s = s.cv, mode = "fraction"),
                                        coef(LASSO.Com, s = s.cv.Com, mode = "fraction"),
                                        coef(IPW_LASSO.Known.Prob, s = s.cv.IPW_LASSO.Known.Prob, mode = "fraction"),
                                        coef(IPW_LASSO.Estimated.Prob, s = s.cv.IPW_LASSO.Estimated.Prob, mode = "fraction"),
                                        coef(IPW_LASSO.Estimated.Prob_lasso, s = s.cv.IPW_LASSO.Estimated.Prob_lasso, mode = "fraction")),
                                      ncol = 5,byrow = F))
  
  names(Ulti.R.Table) = c("Lasso On Full Dataset",                ## Giving names of the columns of this table
                          "Lasso on Complete Case",
                          "IPW-LASSO With Known Prob",
                          "IPW-LASSO With Estimated Prob by MLE",
                          "IPW-LASSO With Estimated Prob by Logistic LASSO")
  
  
  ##--- OVERALL OUTPUT RESULT---##
  Output.Result=list("Missing Probability"= p,
                     "Missing Percentage"= percentage,
                     "Missing Indicator" = Missing.Indecator,
                     "Missing Data set" = Missing.Data,
                     "MLE of Psi" = coefficients(Logistic.Fit),
                     "Logistic Lasso Estimate of Psi" = Cof,
                     "Estimated Missing Probabilities by MLE" = pi_hat,
                     "Estimated Missing Probabilities by Logistic LASSO" = pi_hat_Lasso,
                     "Ultimate Output" = Ulti.R.Table)
  
  return(Output.Result)
}

#########################################################################################################################################################
##WE.R = Our.Method(Our.Original.Data = Overall.Dataset[[400]],Logistic.Coeficient = matrix(c(rep(0,3),0.5,0.5,0,rep(0,6)),nrow=12,ncol = 1,byrow = T))
##WE.R$`Missing Percentage`
##View(WE.R$`Ultimate Output`)

#(D)

##################################################################
##-- Applying this function on the 500 simulated data set and --## 
##--             calculating necessary results.               --##
##################################################################


All.Coeficient.Result = vector(mode = "list",length = 500)
All.Missing.Percentages = NULL
for(l in 1:500){
  
  Result = Our.Method(Our.Original.Data = Overall.Dataset[[l]],
                      Logistic.Coeficient = matrix(c(rep(0,3),0.5,0.5,0,rep(0,6)),
                                                   nrow=12,ncol = 1,byrow = T))
  
  
  All.Coeficient.Result[[l]] = Result$`Ultimate Output`       ## all LASSO coeficients for 500 simulated dataset
  
  All.Missing.Percentages[l] = Result$`Missing Percentage`    ## Missing Percentage for 500 simulated datasets
}

##################################################################################################################################################

#(E)

##########################################
##-- Constructing The necessary Table --##
##########################################

##(1) Table of Miss-Match for 500 simulation

## TRUE  => SELECTED
## FALSE => REJECTED

Original.beta.indicator = c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE)

Miss.Match.Table.For.all.LASSO = list(Miss.Match.Matrix.for_LASSO.FULL = NULL,
                                      Miss.Match.Matrix.for_LASSO.Com = NULL,
                                      Miss.Match.Matrix.for_IPW_LASSO.for.Known.Prob = NULL,
                                      Miss.Match.Matrix.for_IPW_LASSO.for.Estimated.MLE.Prob = NULL,
                                      Miss.Match.Matrix.for_IPW_LASSO.for.Estimated.Logistic.LASSO.Prob = NULL)

for(k in 1:5){
  
  Miss.Match.Matrix = matrix(NA,500,12)
  
  for(i in 1:500){
    
    for(j in 1:12){
      
      if(All.Coeficient.Result[[i]][j,k]!=0){
        ANS = TRUE
      }
      else{
        ANS = FALSE 
      }
      
      if(ANS==Original.beta.indicator[j]){
        
        Miss.Match.Matrix[i,j] = 0                 ## 0 => Correctly Selected
      }
      else{
        Miss.Match.Matrix[i,j] = 1                 ## 1 => Falsely selected
      }
    }
  }
  
  Miss.Match.Table = as.data.frame(Miss.Match.Matrix,row.names = c(paste("Sim",1:500,sep = "_")))
  names(Miss.Match.Table)=c(paste("beta",1:12,sep = "_"))
  
  Miss.Match.Table.For.all.LASSO[[k]] = Miss.Match.Table  ## Miss-match indicator table for all Lasso coef. 
  ## for 500 simulated data sets for different methods.
}



##(2) Creating Beta hat Estimated values 

Beta_VAL_for_all_500_simulation_for_all_5_type_LASSO = list("LASSO.FULL" = NULL,
                                                            "LASSO.COMPLETE" = NULL,
                                                            "IPW-LASSO for Known Prob." =NULL,
                                                            "IPW-LASSO for MLE Prob." = NULL,
                                                            "IPW-LASSO for Logistic-LASSO Prob." =NULL)

for(j in 1:5){
  
  beta_1_Val = NULL
  beta_2_Val = NULL
  beta_3_Val = NULL
  beta_4_Val = NULL
  beta_5_Val = NULL
  beta_6_Val = NULL
  beta_7_Val = NULL
  beta_8_Val = NULL
  beta_9_Val = NULL
  beta_10_Val = NULL
  beta_11_Val = NULL
  beta_12_Val = NULL
  
  for(i in 1:500){
    
    beta_1_Val[i] = All.Coeficient.Result[[i]][1,j]
    beta_2_Val[i] = All.Coeficient.Result[[i]][2,j]
    beta_3_Val[i] = All.Coeficient.Result[[i]][3,j]
    beta_4_Val[i] = All.Coeficient.Result[[i]][4,j]
    beta_5_Val[i] = All.Coeficient.Result[[i]][5,j]
    beta_6_Val[i] = All.Coeficient.Result[[i]][6,j]
    beta_7_Val[i] = All.Coeficient.Result[[i]][7,j]
    beta_8_Val[i] = All.Coeficient.Result[[i]][8,j]
    beta_9_Val[i] = All.Coeficient.Result[[i]][9,j]
    beta_10_Val[i] = All.Coeficient.Result[[i]][10,j]
    beta_11_Val[i] = All.Coeficient.Result[[i]][11,j]
    beta_12_Val[i] = All.Coeficient.Result[[i]][12,j]
    
  }
  
  
  Beta_VAL_for_all_500_simulation =as.data.frame(cbind(beta_1_Val,
                                                       beta_2_Val,
                                                       beta_3_Val,
                                                       beta_4_Val,
                                                       beta_5_Val,
                                                       beta_6_Val,
                                                       beta_7_Val,
                                                       beta_8_Val,
                                                       beta_9_Val,
                                                       beta_10_Val,
                                                       beta_11_Val,
                                                       beta_12_Val))
  
  names(Beta_VAL_for_all_500_simulation) = c(paste("beta",1:12,sep = "_"))
  
  Beta_VAL_for_all_500_simulation_for_all_5_type_LASSO[[j]] = Beta_VAL_for_all_500_simulation
  
}
View(Beta_VAL_for_all_500_simulation_for_all_5_type_LASSO$LASSO.FULL)

Beta_hat_matrix = matrix(NA,12,6)
for(i in 1:5){
  beta_mean = NULL
  beta_mean = colMeans(Beta_VAL_for_all_500_simulation_for_all_5_type_LASSO[[i]])
  Beta_hat_matrix[,i+1] = beta_mean
}
Beta_hat_matrix[,1] = c(1,1,0,0,2.5,1,0,1,0,0,0,0)
Beta_hat_Estimate = as.data.frame(Beta_hat_matrix,
                                  row.names = paste("beta",1:12,sep = "_"))  ##---- Beta hat average table
names(Beta_hat_Estimate) = c(c("true Beta",
                               "Lasso On Full Dataset",
                               "Lasso on Complete Case",
                               "IPW-LASSO With Known Prob",
                               "IPW-LASSO With Estimated Prob by MLE",
                               "IPW-LASSO With Estimated Prob by Logistic LASSO"))


#View(Beta_hat_Estimate)    
View(Beta_hat_Estimate)
##(3) false selection Table

False_Selection_Matrix = matrix(NA,12,5)
for(i in 1:5){
  False_retio = NULL
  False_retio = colMeans(Miss.Match.Table.For.all.LASSO[[i]])
  False_Selection_Matrix[,i] = round(False_retio,4) 
}
False_Selection_Table = as.data.frame(False_Selection_Matrix,row.names = paste("beta",1:12,sep = "_"))
names(False_Selection_Table) = c(c("Lasso On Full Dataset",
                                   "Lasso on Complete Case",
                                   "IPW-LASSO With Known Prob",
                                   "IPW-LASSO With Estimated Prob by MLE",
                                   "IPW-LASSO With Estimated Prob by Logistic LASSO"))

View(False_Selection_Table)
class(False_Selection_Table)

View(Miss.Match.Table.For.all.LASSO$Miss.Match.Matrix.for_LASSO.FULL)

##(4) Average False Selection for different methods.
round(colMeans(False_Selection_Table),4) 
