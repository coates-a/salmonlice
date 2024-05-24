#' Farm treatment log 
#' 
#' Records each time-step and farm that delousing treatment is applied

#' @param ts Number of time steps in simulation
#' @param z Number of lice groups per farm (number stages * number genotypes)
#' @param size Number of lice groups over all farms (z * number farms)
#' @param years.start Number of spin-up years (Change from 1st (starting) to 2nd distribution after t > years.start
#' 
#' @param N Louse population matrix. Number of lice per farm, stage and genotype (rows), per time-step (cols)
#' @param l Lice limit matrix. Number of adult lice per farm that triggers treatment (l[,1]), spring treatment (l[,2]), forced harvest (l[,3])
#' @param V Vector of seasons corresponding to each time-step in simulation
#' @param X Vector of survival on farm (for each life stage and genotype) after treatment X
#' @param X Vector of survival on farm (for each life stage and genotype) after treatment Y
#' @param Zdisc Vector of survival (each life stage and genotype) after treatment Z (always discrete)
#' @param MAX Vector of survival on farm after forced harvest
#' @param treat.log Matrix in which to record when treatments are applied, for each farm (rows), and time-step (cols)
#' @param prod.cycle Matrix of farm production cycle - whether active (1) or inactive/fallow (0)
#' 
#'@return Populate treat.log with 'N' (no treatment), 'T' (treatment) or 'H' (harvest); save as data frame
#'
#'@param export

treat.log_fun <- function(N, l, V, X, Y, Zdisc, MAX, treat.log, prod.cycle) {
  r <- seq(from=1, to=size, by=z) # Create marker for the 1st row of each farm (RR larvae) in N and M
  for(tt in 2:ts){
    
    #### Each time-step build transition matrix M corresponding to the season ####
    season.t <- V[,tt-1] # Determine the name of the season at t-1
    
    d <- if (season.t == "w1" | season.t == "w2" | # Determine which dispersal matrix to use, according to season at t-1
             season.t == "s4" | season.t == "s5" | 
             season.t == "w4") dw else ds 
    
    t.yr <- tt/52 # Convert current time-step to year
    
    ## Assign either first (starting) or second Distribution depending on whether t has got past the spin-up period ('years.start')
    
    H.t <- as.matrix(if (t.yr >= years.start ) H.vector else H.vector.start)
    
    H.gen.t <- if (t.yr >= years.start ) H.gen else H.gen.start
    
    # Assign whether matrices of fecundity (by genotype) are from the starting or final distributions (according to t)
    fRRTT.t <- as.matrix(if (t.yr >= years.start ) fRRTT.h else fRRTT.h.start)
    fRSTT.t <- as.matrix(if (t.yr >= years.start ) fRSTT.h else fRSTT.h.start)
    fSSTT.t <- as.matrix(if (t.yr >= years.start ) fSSTT.h else fSSTT.h.start)
    fRRTU.t <- as.matrix(if (t.yr >= years.start ) fRRTU.h else fRRTU.h.start)
    fRSTU.t <- as.matrix(if (t.yr >= years.start ) fRSTU.h else fRSTU.h.start)
    fSSTU.t <- as.matrix(if (t.yr >= years.start ) fSSTU.h else fSSTU.h.start)
    fRRUU.t <- as.matrix(if (t.yr >= years.start ) fRRUU.h else fRRUU.h.start)
    fRSUU.t <- as.matrix(if (t.yr >= years.start ) fRSUU.h else fRSUU.h.start)
    fSSUU.t <- as.matrix(if (t.yr >= years.start ) fSSUU.h else fSSUU.h.start)
    
    # Assign whether matrices of at.yrachment (by genotype) are from the starting or final distributions (according to t)
    vRRTT.t <- as.matrix(if (t.yr >= years.start ) vRRTT.h else vRRTT.h.start)
    vRSTT.t <- as.matrix(if (t.yr >= years.start ) vRSTT.h else vRSTT.h.start)
    vSSTT.t <- as.matrix(if (t.yr >= years.start ) vSSTT.h else vSSTT.h.start)
    vRRTU.t <- as.matrix(if (t.yr >= years.start ) vRRTU.h else vRRTU.h.start)
    vRSTU.t <- as.matrix(if (t.yr >= years.start ) vRSTU.h else vRSTU.h.start)
    vSSTU.t <- as.matrix(if (t.yr >= years.start ) vSSTU.h else vSSTU.h.start)
    vRRUU.t <- as.matrix(if (t.yr >= years.start ) vRRUU.h else vRRUU.h.start)
    vRSUU.t <- as.matrix(if (t.yr >= years.start ) vRSUU.h else vRSUU.h.start)
    vSSUU.t <- as.matrix(if (t.yr >= years.start ) vSSUU.h else vSSUU.h.start)
    
    ###############
    
    M <- M_fun(season.t, d, t.yr) # Build transition matrix, M, for t
    
    #### Calculate the number of alleles produced by adults on each farm ####
    A <- (N[r+3,tt-1]+N[r+7,tt-1]+N[r+11,tt-1]+
            N[r+15,tt-1]+N[r+19,tt-1]+N[r+23,tt-1]+
            N[r+27,tt-1]+N[r+31,tt-1]+N[r+35,tt-1]) # Calculate total number of adults (across all genotypes) per farm at t-1
    
    # For each farm, calculate the number of offspring produced per adult for each genotype (factoring in temperature effects, and any effect of continuous strategies)
    f.tRRTT <- f.temp[, season.t] * # Reproductive output per RRTT adult per farm
      fRRTT.t # multiplied by proportion of f produced by RRTT according to continuous strategy and/or trade-offs
    f.tRSTT <- f.temp[, season.t] * 
      fRSTT.t 
    f.tSSTT <- f.temp[, season.t] * 
      fSSTT.t
    f.tRRTU <- f.temp[, season.t] * # Reproductive output per RRTU adult per farm
      fRRTU.t # multiplied by proportion of f produced by RRTU according to continuous strategy and/or trade-offs
    f.tRSTU <- f.temp[, season.t] * 
      fRSTU.t 
    f.tSSTU <- f.temp[, season.t] * 
      fSSTU.t
    f.tRRUU <- f.temp[, season.t] * # Reproductive output per RRUU adult per farm
      fRRUU.t # multiplied by proportion of f produced by RRUU according to continuous strategy and/or trade-offs
    f.tRSUU <- f.temp[, season.t] * 
      fRSUU.t 
    f.tSSUU <- f.temp[, season.t] * 
      fSSUU.t
    
    
    RRTTad <- N[r+3,tt-1] # Number of RRTT adults per farm at t-1
    RSTTad <- N[r+7,tt-1] # Number of RSTT adults per farm at t-1
    SSTTad <- N[r+11,tt-1] # Number of SSTT adults per farm at t-1
    RRTUad <- N[r+15,tt-1] 
    RSTUad <- N[r+19,tt-1] 
    SSTUad <- N[r+23,tt-1] 
    RRUUad <- N[r+27,tt-1] 
    RSUUad <- N[r+31,tt-1] 
    SSUUad <- N[r+35,tt-1] 
    
    # Calculate number of R alleles produced per farm at t-1
    n.R <- (2 * RRTTad * # 2 R gametes produced by adults of each RR genotype, for each offspring
              f.tRRTT) + # factoring in the reproductive output of each genotype
      (2 * RRTUad * f.tRRTU) +
      (2 * RRUUad * f.tRRUU) +
      (RSTTad * f.tRSTT) + # 1 R gamete produced by adults of each RS genotype, for each offspring
      (RSTUad * f.tRSTU) +
      (RSUUad * f.tRSUU) 
    # Calculate number of S alleles produced per farm at t-1
    n.S <- (2 * SSTTad * # 2 S gametes produced by adults of each SS genotype, for each offspring
              f.tSSTT) + # factoring in the reproductive output of each genotype
      (2 * SSTUad * f.tSSTU) +
      (2 * SSUUad * f.tSSUU) +
      (RSTTad * f.tRSTT) + # 1 S gamete produced by adults of each RS genotype, for each offspring
      (RSTUad * f.tRSTU) +
      (RSUUad * f.tRSUU)
    # Calculate number of T alleles produced per farm at t-1
    n.T <- (2 * RRTTad * # 2 T gametes produced by adults of each TT genotype, for each offspring
              f.tRRTT) + # factoring in the reproductive output of each genotype
      (2 * RSTTad * f.tRSTT) +
      (2 * SSTTad * f.tSSTT) +
      (RRTUad * f.tRRTU) + # 1 T gamete produced by adults of each TU genotype, for each offspring
      (RSTUad * f.tRSTU) +
      (SSTUad * f.tSSTU) 
    # Calculate number of U alleles produced per farm at t-1
    n.U <- (2 * RRUUad * 
              f.tRRUU) + 
      (2 * RSUUad * f.tRSUU) +
      (2 * SSUUad * f.tSSUU) +
      (RRTUad * f.tRRTU) + 
      (RSTUad * f.tRSTU) +
      (SSTUad * f.tSSTU)
    
    
    
    #### Calculate gene frequencies of R (pR), S (qS), T (pT), U (qU), in gamete pool ####
    # Create vectors for p and q for each farm at t-1
    pR <- ifelse((n.R + n.S) > 0, # Don't want to divide by 0 if n.R + n.S = 0 (i.e. no adults reproducing at farm)
                 pR <- (n.R/(n.R + n.S)), # Proportion of R alleles at each farm
                 pR <- 0.5) # If no alleles produce --> Doesn't matter what this value is, as it will be multiplied by 0.
    qS <- ifelse((n.R + n.S) > 0, # Don't want to divide by 0 if (n.R + n.S) = 0 
                 qS <- (n.S/(n.R + n.S)), # Frequency of S allele at each farm
                 qS <- 0.5) # If no alleles produce --> Doesn't matter what this value is, as it will be multiplied by 0.
    pT <- ifelse((n.T + n.U) > 0, # Don't want to divide by 0 if (n.R + n.S) = 0 
                 pT <- (n.T/(n.T + n.U)), # Frequency of T allele at each farm
                 pT <- 0.5) # If no alleles produce --> Doesn't matter what this value is, as it will be multiplied by 0.
    qU <- ifelse((n.T + n.U) > 0, # Don't want to divide by 0 if (n.R + n.S) = 0 
                 qU <- (n.U/(n.T + n.U)), # Frequency of U allele at each farm
                 qU <- 0.5)
    # p + q = 1
    
    # Multiply total reproductive output per adult per genotype per farm (values already in M) by expected gene frequencies (from Hardy-Weinburg principle)
    # Insert pR^2 * pT^2 into RRTT LARVA rows
    M[r,] <- sweep(M[r,], MARGIN=1, (pR^2)*(pT^2), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by pR^2 * pT^2 for that farm
    # RSTT Larva rows (2pRqS)*pT^2
    M[r+4,] <- sweep(M[r+4,], MARGIN=1, (2*pR*qS)*(pT^2), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by 2pq for that farm
    # SSTT rows
    M[r+8,] <- sweep(M[r+8,], MARGIN=1, (qS^2)*(pT^2), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by q^2 for that farm
    
    # RRTU
    M[r+12,] <- sweep(M[r+12,], MARGIN=1, (pR^2)*(2*pT*qU), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by pR^2 * pT^2 for that farm
    # RSTU
    M[r+16,] <- sweep(M[r+16,], MARGIN=1, (2*pR*qS)*(2*pT*qU), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by 2pq for that farm
    # SSTU
    M[r+20,] <- sweep(M[r+20,], MARGIN=1, (qS^2)*(2*pT*qU), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by q^2 for that farm
    
    # RRUU
    M[r+24,] <- sweep(M[r+24,], MARGIN=1, (pR^2)*(qU^2), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by pR^2 * pT^2 for that farm
    # RSUU
    M[r+28,] <- sweep(M[r+28,], MARGIN=1, (2*pR*qS)*(qU^2), `*`) # Sweep multiplies (*) entire row (MARGIN=1) by 2pq for that farm
    # SSUU
    M[r+32,] <- sweep(M[r+32,], MARGIN=1, (qS^2)*(qU^2), `*`) 
    
    
    #### Calculate number of lice at t (before anti-lice treatments) ####
    N[,tt] <- as.matrix(M) %*% N[,tt-1] # Multiply N(t-1) by the newly-adjusted M (need to convert sparse matrix M to full matrix)
    
    ##### APPLY DELOUSING TREATMENTS #####
    
    At <- (N[r+3,tt]+N[r+7,tt]+N[r+11,tt]+
             N[r+15,tt]+N[r+19,tt]+N[r+23,tt]+
             N[r+27,tt]+N[r+31,tt]+N[r+35,tt]) # Calcaulte total number of adults per farm at t (before treatment)
    
    # Assign lice limit, according to season at t
    ifelse (season.t == "sp", # If spring
            ll <- l[,2], # use spring limit,
            ll <- l[,1]) # otherwise use limit for rest of year
    
    ############################################################ 
    ################# CODE FOR treatments_Fun ##################
    
    # Record which farms receive treatment
    treat.log[,tt] <- ifelse(At < ll, # If total number of adults is under limit,
                             "N", # Insert 'N' (no treatment)
                             ifelse(At > l[,3], # If over maximum limit
                                    "H", # Insert 'H' (forced harvest) 
                                    ifelse(((H.t=="x" | H.t=="xy") & stratx=="disc"), # If adults in range of treatment trigger AND on farm where a discrete strategy is used,
                                           "X", # INsert 'X'
                                           ifelse(((H.t=="y" | H.t=="xy") & straty=="disc"), 
                                                  "Y", # Same for it discrete Y is used
                                                  ifelse((H.t=="z" | H.t=="xz" | H.t=="yz" | H.t=="xyz"), # If on farm using discrete Z
                                                         "Z", # Insert 'X'
                                                         "N")))))
    
    ############################################################
    
    # Create vector for treatment survival for strategy X (if X is discrete) that runs the WHOLE length of N
    Xt <- ifelse(rep(At, each=z) < rep(ll, each=z), # If number of adults is under lice limit (need to repeat each vector values so same length as z)
                 
                 1, # No treatment (survival = 1)
                 
                 ifelse(rep(At, each=z) < rep(l[,3], each=z), # If number of adults over limit but under max limit,
                        
                        ifelse((rep(H.t, each=z)=="x" | rep(H.t, each=z)=="xy" | rep(H.t, each=z)=="xz" | rep(H.t, each=z)=="xyz") & # AND if farm is using X AND strat X is discrete
                                 stratx=="disc",
                               
                               X, # Insert vector of X treatment survival (survival for each stage and genotype)
                               
                               1), # If no X or X isn't discrete, no treatment
                        
                        MAX)) # If adults OVER max limit, insert vector of forced harvest mortality
    
    # Vector treaetment survival from strategy Y
    
    Yt <- ifelse(rep(At, each=z) < rep(ll, each=z), # If number of adults is under lice limit (need to repeat each vector values so same length as z)
                 
                 1, # No treatment (survival = 1)
                 
                 ifelse(rep(At, each=z) < rep(l[,3], each=z), # If number of adults over limit but under max limit,
                        
                        ifelse((rep(H.t, each=z)=="y" | rep(H.t, each=z)=="xy" | rep(H.t, each=z)=="yz" | rep(H.t, each=z)=="xyz") & # AND farm is using Y AND strat Y is discrete
                                 straty=="disc",
                               
                               Y, # Insert vector of X treatment survival (survival for each stage and genotype)
                               
                               1), # If no X or X isn't discrete, no treatment
                        
                        MAX)) 
    # Vector treaetment survival from strategy Z
    
    Zt <- ifelse(rep(At, each=z) < rep(ll, each=z), # If number of adults is under lice limit (need to repeat each vector values so same length as z)
                 
                 1, # No treatment (survival = 1)
                 
                 ifelse(rep(At, each=z) < rep(l[,3], each=z), # If number of adults over limit but under max limit,
                        
                        ifelse((rep(H.t, each=z)=="z" | rep(H.t, each=z)=="xz" | rep(H.t, each=z)=="yz" | rep(H.t, each=z)=="xyz"), # AND farm is using Y AND strat Y is discrete
                               
                               Zdisc, # Insert vector of X treatment survival (survival for each stage and genotype)
                               
                               1), # If no Z, no treatment
                        
                        MAX)) 
    
    N[,tt] <- N[,tt] * Xt * Yt * Zt # Multiply number of lice at t by proportion treatment survival (for each strategy)
    
    N[,tt] <- floor(N[,tt]) # Round down to nearest integer (Adults with N<1 don't reproduce)
    
    # Production Cycle active vs. fallow farms
    PC <- prod.cycle[,tt] # Vector of whether farms are active (1) or inactive (0) at tt
    PC <- matrix(rep(PC, each=z), # Repeat for each genotype & stage on farm
                 nrow=size,
                 ncol=1)
    PC[c(seq(1, size, by=4)),] <- 1 # Put 1 in all Larva rows (larvae already produced at t-1 and are in environment; unaffected by inactive farms)
    
    N[,tt] <- N[,tt] * PC # Multiply farms by 0 if inactive; by 1 in active
    
  }
  ############################################################ 
  
  return(treat.log) # Return matrix of farm treatments
  
  ############################################################
}
