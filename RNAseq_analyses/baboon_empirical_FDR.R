library(biomaRt)

#permutation function
perm.fdr=function(input_df,perm_df,Pvals_col_name,plot=T,details=T){
  library(ggplot2)
  library(reshape)
  library(qvalue)
  library(minpack.lm)
  options(width=10000,max.print=1000000000)
  
  min=0.0000001
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df
  rownames(ro)=rownames(input_df)
  ro<-ro[order(ro[,pvals_index]),]
  names=rownames(ro)
  #names=rownames(input_df)
  
  
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  n_perm=ncol(perm_df)
  n_genes=nrow(input_df[,1])
  p_vector<-perm_df[,1]
  if(n_perm>1){
    for(i in 2:n_perm){
      p_vector<-rbind(p_vector,perm_df[,i])
    }}
  p_vector<-data.frame(p_vector[order(p_vector)])
  
  if(plot==T)
  {
    comando=("mkdir -p outputs/figures/")
    system(comando)
  }
  if(details==T)
  {
    comando=("mkdir -p outputs/files/")
    system(comando)
  }
  
  
  
  p_rand=data.frame(p_vector)
  p_rand$p_common='p_random'
  colnames(p_rand)=c("pvalue","p_common")
  p_obs$p_common='p_observed'
  
  p_rand=data.frame(p_vector)
  p_rand$p_common='p_random'
  colnames(p_rand)=c("pvalue","p_common")
  p_obs$p_common='p_observed'
  
  p_hists <- rbind(p_rand,p_obs)
  if(plot==T)
  {
    densidades<-ggplot(p_hists, aes(pvalue, fill = p_common)) + geom_density(alpha = 0.4)
    pdf(paste0(getwd(), "/outputs/figures/P_value_densities.pdf"))
    print(densidades)
    dev.off()
    histogramas<-ggplot(p_hists, aes(pvalue, fill = p_common)) + geom_histogram(alpha = 0.4, aes(y = ..density..), binwidth=.01, position = 'identity')
    pdf(paste0(getwd(), "/outputs/figures/P_value_histograms.pdf"))
    print(histogramas)
    dev.off()
  }
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  print("Obtaining p-value distributions and estimators for the fraction of true nulls...")
  for(i in 1:observed)
  {
    #  print(c(i,j))
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
  }
  pi_hat=(1-F)/(1-F_o)
  Fdr_BH_perm=F_o/F
  
  
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  general_scaled_logistic <- function(pars, x) {(pars$q^pars$n-pars$k)/(pars$q^pars$n-1)+((pars$k-1)*pars$q^pars$n)/((pars$q^pars$n-1)*((1-(1-pars$q)*exp(-tan(pi*x/2)))^pars$n))}
  residFun <- function(p, observed, xx) observed - general_scaled_logistic(p,xx)
  parStart <- list(q=2,n=1,k=0.9)
  print("Fitting true null fraction...")
  
  f_hat <- nls.lm(par=parStart, fn = residFun, observed = tabla$pi_hat,xx = tabla$pval,control = nls.lm.control(nprint=-1,maxiter=1000))
  f_hat_serie<-general_scaled_logistic(f_hat$par,p_obs[,1])
  
  
  
  pi_o=min(f_hat$par[[3]],1)
  pi_o=max(f_hat$par[[3]],0)
  
  
  
  Fndr_ST_perm=1-pi_o/pi_hat
  print("Calculating Fdrs and Fndrs...")
  
  for(i in 1:length(p_obs[,1]))
  {
    #	print(i)
    if(Fdr_BH_perm[i]>1)
      Fdr_BH_perm[i]=1
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_BH_perm[i-j]>Fdr_BH_perm[i]){Fdr_BH_perm[i-j]=Fdr_BH_perm[i]}else{break}
      }
    }
    #Fndr_unif[i]=(1-pi_o_un-F[i]+pi_o_un*p_obs[i,1])/(1-F[i])
    if(Fndr_ST_perm[i]<0){Fndr_ST_perm[i]=0}
    
  }
  Fdr_ST_perm=pi_o*Fdr_BH_perm
  
  for(i in length(p_obs[,1]):1)
  {
    #	print(i)
    if(i<length(p_obs[,1]))
    {
      for(j in 1:(length(p_obs[,1])-i))
      {
        if(Fndr_ST_perm[i+j]>Fndr_ST_perm[i]){Fndr_ST_perm[i+j]=Fndr_ST_perm[i]}else{break}
      }
    }
  }
  
  if(details==T)
  {
    sink(paste(getwd(),"/outputs/files/true_null_fractions_fit_info.txt",sep=""))
    print(paste("PI_o from permuted null: ",pi_o))
    print("PI_o fit to generalized, re-scaled logistic function summary: ")
    print(summary(f_hat))
    sink()
  }
  
  
  fdrs_df <-data.frame(ro,q_BH_perm=Fdr_BH_perm,q_ST_perm=Fdr_ST_perm,Fndr_ST_perm=Fndr_ST_perm)
  rownames(fdrs_df)=names
  #sink(paste(getwd(),"/outputs/files/fdrs.txt",sep=""))
  #fdrs_df
  #sink()
  #rownames(fdrs_df)=names
  fdrs_df=fdrs_df[order(rownames(fdrs_df)),]
  
  if(plot==T)
  {
    denom=as.integer(length(fdrs_df[,1])/2000)
    indices=c(1:length(fdrs_df[,1]))
    bins<-which(indices%%denom==0)
    fdrs_df_subsetted<-rbind(fdrs_df[bins,],fdrs_df[c((length(fdrs_df[,1])-200):length(fdrs_df[,1])),])
    pval_column=which(colnames(fdrs_df_subsetted)==Pvals_col_name)
    fdrs_df_subsetted <- melt(fdrs_df_subsetted[,c(pval_column,(length(ro)+1):length(fdrs_df_subsetted))] ,  id = Pvals_col_name, variable_name = 'series')
    
    # plot on same grid, each series colored differently -- 
    # good if the series have same scale
    
    colnames(fdrs_df_subsetted)[1]="P.Values"
    fdrs_plot<-ggplot(fdrs_df_subsetted, aes(P.Values,value)) + geom_line(aes(colour = series))
    pdf(paste0(getwd(), "/outputs/figures/fdrs.pdf"))
    print(fdrs_plot)
    dev.off()
  }
  
  if(details==T)
  {
    fdrs<-data.frame(F,F_o,Fdr_BH_perm,Fdr_ST_perm,Fndr_ST_perm,pi_hat_fit=f_hat_serie,pi_hat,pval=p_obs[,1])
    rownames(fdrs)=names
    #sink(paste(getwd(),"/outputs/files/fdrs_building.txt",sep=""))
    #fdrs
    #sink()
    #rownames(fdrs)=rownames(ro)
    write.table(fdrs,paste(getwd(),"/outputs/files/fdrs_building.txt",sep=""))
  }
  
  if(plot==T & details==T)
  {
    denom=as.integer(length(fdrs[,1])/2000)
    indices=c(1:length(fdrs[,1]))
    bins<-which(indices%%denom==0)
    fdrs<-rbind(fdrs[bins,],fdrs[c((length(fdrs[,1])-200):length(fdrs[,1])),])
    
    fdrs <- melt(fdrs ,  id = "pval", variable_name = 'series')
    
    # plot on same grid, each series colored differently -- 
    # good if the series have same scale
    fdrs_building<-ggplot(fdrs, aes(pval,value)) + geom_line(aes(colour = series))+coord_cartesian(ylim = c(-0.01,1.01))
    
    pdf(paste0(getwd(), "/outputs/figures/fdrs_building.pdf"))
    print(fdrs_building)
    dev.off()
  }
  print("Done.")
  return(fdrs_df)
}

#output from linear modeling
results=read.delim('BaboonLMEOutput2-1-2023.txt')

#p-value permutations
GeomPolyTable <- read.delim("BaboonLMMEmpiricalFDR.txt")

#turn table into format wanted
PermutationTable<-matrix(ncol=10,nrow= dim(results)[1])
for (i in 1:10){
  PermutationTable[,i]<-subset(GeomPolyTable,Permutation==i)$V14
}

DEGenes<-perm.fdr(results, PermutationTable, "V14")

write.table(DEGenes, "BaboonDEGenesFDR.txt", quote=F)
