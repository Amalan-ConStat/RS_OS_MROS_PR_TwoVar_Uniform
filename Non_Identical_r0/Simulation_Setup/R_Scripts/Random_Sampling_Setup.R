library(here)
library(compiler)
library(LaplacesDemon)

# Cordeiro Bias Estimation ----
Cordeiro<-function(XData,With_bias)
{
  p <- as.vector(exp(XData%*%as.vector(With_bias)))
  W <- diag(p)
  inverse_term <- solve(t(XData)%*%W%*%XData)
  
  Term1 <- inverse_term%*%t(XData)%*%W
  Term2 <- diag(diag(XData%*%(inverse_term)%*%t(XData))) %*% rep(-0.5,nrow(XData))
  
  bias <- as.vector(Term1%*%Term2)
  return(bias)
}

# Using Random Sampling Sub-sample from Big Data ----
Run_RandomSample <- function(Replicates,FullData,Subsample_Size,N,Choices,No_of_Models,Theta)
{
  #RandomSample<-list()
  #Sample_RS<-list()
  
  All_Parameter<-list()
  Est_Parameter <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))
  
  All_Bias<-list()
  Bias <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))
  
  All_optimality<-list()
  optimality <- matrix(NA,nrow = length(Choices),ncol = 3)
  
  Final_Parameter<-Final_optimality<-Final_Bias<-list()
  
  for (k in 1:No_of_Models) 
  {
    cat("Model :",k,"\n")
    for(i in 1:Replicates)
    {
      Sampled<-sample(1:N,size = Subsample_Size)
      Temp_Data<-as.data.frame(FullData[Sampled,])
      #RandomSample[[i]]<-Temp_Data
      
      for (j in 1:length(Choices))
      {
        Temp_Temp_Data<-Temp_Data[1:Choices[j],]
        Results<-glm(Y~.-1,data=Temp_Temp_Data,family = "poisson")
        Est_Parameter[j,]<-c(Choices[j],Results$coefficients)
        pi<-c(exp(as.matrix(Temp_Temp_Data[,-1]) %*% Theta))
        Mx<-t(as.matrix(Temp_Temp_Data[,-1])) %*% (as.matrix(Temp_Temp_Data[,-1]) * pi)
        optimality[j,]<-c(Choices[j], psych::tr(vcov(Results)), det(Mx) )
        Bias[j,]<-c(Choices[j],Cordeiro(XData = as.matrix(Temp_Temp_Data[,-1]),
                                        With_bias = as.vector(Results$coefficients)))
      }
      
      All_Parameter[[i]]<-Est_Parameter
      All_optimality[[i]]<-optimality
      All_Bias[[i]]<-Bias
      #Sample_RS[[i]]<-cbind(i,RandomSample[[i]])
      print(i)
    }
    
    Final_Parameter[[k]]<-do.call(rbind,All_Parameter)
    Final_optimality[[k]]<-do.call(rbind,All_optimality)
    Final_Bias[[k]]<-do.call(rbind,All_Bias)
    #Final_Sample_RS<-do.call(rbind,Sample_RS)
  }
  
  return(list("EstPar"=Final_Parameter,
              "Utility"=Final_optimality,
              "Bias"=Final_Bias#,"Subsampled_Data"=Final_Sample_RS
  )
  )
}

Run_RandomSample<-cmpfun(Run_RandomSample)

#Save the Random Sample function for No Correlated Covariate Data ----
save(Run_RandomSample,Cordeiro,
     file=here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling","Run_RandomSample.RData"))

rm(list = ls())

# Run the Random sampling method ----
source(here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling",
            "Simulation_Results_Random_Sampling.R"))
