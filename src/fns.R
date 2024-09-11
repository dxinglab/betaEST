############################################################################
##                            main functions                              ##
############################################################################

############################################################################
## function for calculating observed beta from a community composition    ##
## matrix                                                                 ##
############################################################################
# new unbiased beta proposed by D. Xing (unpublished):
betax=function(comp){
    #@comp S*M community matrix containing abundance data
    comp=comp[rowSums(comp)>0,]
    occ_dat=comp
    occ_dat[occ_dat>0]=1
    M=ncol(occ_dat)
    S=nrow(occ_dat)
    q=1 - rowSums(occ_dat)/M # prob. of absence
    p=rowSums(comp)/sum(comp) #relative abundance
    
    res=sum(p*q)*M/(M-1)
    
    c(betax=res)
}
# whittaker's classical beta:
betap=function(comp, type='p'){
    # @comp - S*M community matrix containing either occurrence or abundance data
    comp[comp>0]=1 # convert it to occurrence if comp is an abundance matrix
    if (any(rowSums(comp)==0)|any(colSums(comp)==0))
        stop('check comp: 0 row/col sums detected.')
    S=nrow(comp)
    M=ncol(comp)
    alphaMean=mean(colSums(comp))
    res=1-alphaMean/S
    if(type=='p')
        return(c(betap=res))
    else if(type=='w')
        return(c(betaw=1/(1-res)))
}
############################################################################
## function for calculating analytical beta deviation, from 0/1 data      ##
############################################################################
beta_dev=function(comp,N=NULL, S=NULL, version='precise'){
    # @comp - S*M community matrix containing either occurence or abundance data
    # @N - metacommunity size (total abundance), 
    #  is either specified (in cases comp is occurence data)
    #  or calculated from the community matrix (when comp is abundance data).
    # @S - number of species observed in the metacommunity. calculated from the comp data if not given.
    # @version - the method used to find the logseries parameter
    if(is.null(N)){
        if(all(comp<=1))
            stop('N not available.')
        else N=sum(comp)
    }
    M=ncol(comp)
    S=nrow(comp)
    (betap(comp) - beta_null(M,S,N,version)) / beta_null_sd(M,S,N,version)
}



############################################################################
##                other functions used in the main functions              ##
############################################################################

############################################################################
## function for finding the mete-derived logseries parameter (lambda)     ##
## based on state variables S and N (sensu Harte 2011)                    ##
############################################################################
meteSN=function(S0, N0, version='precise'){
    if(!length(S0) == length(N0)) stop("S and N must have the same length")
    if(!all(S0[] > 1,na.rm=T)) stop("S must be greater than 1")
    if(!all(N0[] > 0,na.rm=T)) stop("N must be greater than 0")
    if(!all(S0[]/N0[] < 1,na.rm=T)) stop("N must be greater than S")
    if(!version %in% c('precise', 'approx')) stop("version must be either 'precise' or 'approx'")
    
    res = rep(NA,length(N0))
    for(i in 1:length(S0)){
        S=S0[i]
        N=N0[i]
        SNid=paste(S,N)
        NSratio=N/S
        if(!is.na(S/N)){
            if(version == 'precise'& N<1e7){
                require(VGAM) # for the lerch() function 
                # eqn 7.29 of Harte (2011)
                y = function(x) S/N * x*(1-x^N)/(1-x) -
                    (-log(1-x)-x^(N+1)*lerch(x,1,N+1,iter = 1000))
                p = try(uniroot(y, lower=.5, upper=1-1e-7,tol=1e-16)$`root`)
                if(class(p)!='try-error')
                    res[i]= -log(p)
            }
            else #if(version == 'approx'){
                {# should really try to estimate alpha, instead of p:
                y= function(x) x*log(1+N/x) -S
                alpha = try(uniroot(y, lower=1, upper=S,tol=1e-15)$`root`)
                if(class(alpha)!='try-error')
                    res[i]= N/(alpha+N)
                # require(pracma) # for the fsolve function
                # # eqn 7.30 of Harte (2011)
                # y = function(x) x * log(1/x) - S / N
                # p = try(uniroot(y, lower=1e-15, upper=1-1e-13,tol=1e-15)$`root`)
                # if(class(p)!='try-error')
                #     res[i]= -log(p)
                # #res[i]=fsolve(y,1/N)
            }
        }
    }
    return (res)
}


############################################################################
##function for calculating analytical null beta (eqn 2a in Xing & He 2021)##
############################################################################
beta_null=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    log(1-p+p/M)/log(1-p)
}
############################################################################
## function for calculating standard deviaiton of null beta               ##
## (eqn 2b in Xing & He 2021)                                             ##
############################################################################
beta_null_sd=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    
    res=1/S/M/log(1-p)*((M-1)*log(1-p*(1-2/M)) +
                            log(1-p*(1-1/M)) -
                            M*log(1-p*(1-1/M)^2))
    sqrt(res)
}

############################################################################
## function for calculating standard deviaiton of beta deviation          ##
## (eqn 3 in Xing & He 2021)                                              ##
############################################################################
beta_dev_sd=function(M,S,N,version='precise'){
    # @M - num of sampling units (plots or quadrats)
    # @S - num of spp
    # @N - num of tot individuals
    # @version - the method used to find the logseries parameter
    lambda=meteSN(S,N,version=version)
    p=exp(-lambda)
    v_Phi=1/S/log(1-p)* (log(1-p*(1-1/M)^2) -
                             1/log(1-p) * log(1-p*(1-1/M))^2 )
    v_Pi=1/S/M/log(1-p)* ((M-1)*log(1-p*(1-2/M)) +
                              log(1-p*(1-1/M)) -
                              M*log(1-p*(1-1/M)^2))
    sqrt(v_Phi/v_Pi)
}
