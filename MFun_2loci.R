#' Build stage transition matrix (M) according to season
#' 
#' Spatial variation in development, survival and reproduction, according to site-specific temperature and whether or not continuous strategy is used on farm 

#' @param season.name Name of 'season' corresponding to current time-step (in popFun function)
#' @param d Louse larval dispersal matrix corresponding to season
#' @param t.yr Current time-step of simulation, converted to years
#' 
#'@return Save completed transition matrix (M)
#'
#'@param export

M_fun <- function(season.name, d, t.yr){ # Name of season and corresponding dispersal matrix at current time-step
  
  M <- kronecker(diag(1,i), G) # Fill out diagonal of M (transition matrix for all farms) with G (transition matrix for single farm)
  
  ##### Assign Distribution based on time
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
  
  #### Fill M with development rates (s) ####
  s.t <- s.temp[, season.name] # Vector of development rates per farm, for current season
  s.t <- c(rep(s.t, each=g)) # Repeat development rates for each genotype
  
  # Adjust growth rate for any effects of continuous strategy and genotype
  s.t <- rep(c(s.t * c(H.gen.t$s.cont)), # Multiply baseline rates by proportion of rates per farm per genotype
             each = j) # Repeat values for each stage, so same length as M
  
  # Fill diagonal of M with 1-s (proportion of lice REMAINING in same stage)
  s.t1 <- s.t
  s.t1[c(seq(4,size,4), seq(1,size,4))] <- 0 # Remove values from Adult and Larva rows (no transition)
  diag(M) <- diag(M) * # Multiply whole diagonal (values won't appear in cells =0)
    (1-s.t1) # by 1-s = proportion of chalimus & pre-adults retained)
  
  # Fill off-diagonal of M with s (proportion of lice TRANSITIONING to next stage)
  s.t2 <- s.t
  s.t2[c(seq(5,size,4))] <- 1 # To fill out off-diagonal, need to put 1s BACK in the rows where f will be
  diag(M[-1,]) <- diag(M[-1,]) * # Multiply OFF-diagonal by s (proportion of chalimus & pre-adults transitioning)
    s.t2[-1] # Need to remove first value from vector, since length of off-diagonal is 1 shorter than diagonal   
  
  #### Fill M with background mortality survival (mu) ####
  
  # Multiply diagonal with rate of survival from background mortality (survival of lice remaining in same stage)
  mu.allD <- rep(muD, i) # Repeat vector of survival for all farms
  # Need to adjust vector of background survival OF CHALIMUS by any effect caused by continuous strategies and.or genotype
  mu.contD <- c(rep(1, size)) # Create vector with length of M diagonal
  mu.contD[c(seq(2,size,4))] <- # For chalimus cells (every 4th cell of diagonal, starting cell 2),
    H.gen.t$mu.cont # Insert proportion of chalimus survival, mu_, to adjust mu according to genotype and whether or not continuous strategy is used
  mu.allD <- mu.allD * mu.contD ## Adjust background-mortality survival according to effects of continuous strategy and genotype
  diag(M) <- diag(M) * mu.allD # Multiply adjusted survival vector along diagonal
  
  # Multiply off-diagonal with rate of survival from background mortality (survival of lice transitioning to next stage)
  mu.allO <- rep(muO, i)
  # Need to adjust vector of background survival by any effect caused by continuous strategies and.or genotype
  mu.contO <- c(rep(1, size))
  mu.contO[c(seq(3,size,4))] <- H.gen.t$mu.cont # Chalimus cells start 3rd cell of off-diagonal
  mu.allO <- mu.allO * mu.contO
  diag(M[-1,]) <- diag(M[-1,]) * # Multiply off-diagonal with vector of adjusted survival
    mu.allO[-1] # Remove first value of vector, since length of off-diagonal is 1 shorter than diagonal
  
  #### Fill M with larval dispersal & attachment (d * v) ####
  dRRTT <- sweep(d, MARGIN=1, vRRTT.t, `*`) # RRTT Dispersal (to each farm, from each farm) * attachment success of RRTT (according to whether continuous strategy is used on farm)
  dRSTT <- sweep(d, MARGIN=1, vRSTT.t, `*`) # RSTT Dispersal * attachment success 
  dSSTT <- sweep(d, MARGIN=1, vSSTT.t, `*`) # SSTT Dispersal * attachment success 
  dRRTU <- sweep(d, MARGIN=1, vRRTU.t, `*`)  
  dRSTU <- sweep(d, MARGIN=1, vRSTU.t, `*`) 
  dSSTU <- sweep(d, MARGIN=1, vSSTU.t, `*`) 
  dRRUU <- sweep(d, MARGIN=1, vRRUU.t, `*`)  
  dRSUU <- sweep(d, MARGIN=1, vRSUU.t, `*`) 
  dSSUU <- sweep(d, MARGIN=1, vSSUU.t, `*`)
  
  # Assign which cells of M to fill with d values for each genotype (Chalimus rows = every 4th row starting row 2; Larva columns = every 4th col starting col 1)
  drowsRRTT <- 2 # RRTT Chalimus row (in 1st farm)
  dcolsRRTT <- 1 # RRTT Larvae column (in 1st farm)
  drowsRSTT <- 6 # RS Chalimus row (in 1st farm)
  dcolsRSTT <- 5 # RS Larvae column (in 1st farm)
  drowsSSTT <- 10 # SS Chalimus row (in 1st farm)
  dcolsSSTT <- 9 # SS Larvae column (in 1st farm)
  
  drowsRRTU <- 14 # 
  dcolsRRTU <- 13 # 
  drowsRSTU <- 18 # 
  dcolsRSTU <- 17 # 
  drowsSSTU <- 22 # 
  dcolsSSTU <- 21 
  
  drowsRRUU <-  26 # 
  dcolsRRUU <-  25 # 
  drowsRSUU <-  30 # 
  dcolsRSUU <-  29 # 
  drowsSSUU <-  34 # 
  dcolsSSUU <-  33 # 
  
  M[seq(from=drowsRRTT, to=size, by=z), # RRTT dispersal TO RRTT Chalimus rows, FROM RRTT Larva columns
    seq(from=dcolsRRTT, to=size, by=z)] <- dRRTT # Dispersal & attachment of RRTT (for each farm)
  
  M[seq(from=drowsRSTT, to=size, by=z), # RSTT dispersal TO RSTT Chalimus rows, FROM RSTT Larva columns
    seq(from=dcolsRSTT, to=size, by=z)] <- dRSTT # Dispersal & attachment of RSTT (for each farm)
  
  M[seq(from=drowsSSTT, to=size, by=z), # SSTT dispersal TO SSTT Chalimus rows, FROM SSTT Larva columns
    seq(from=dcolsSSTT, to=size, by=z)] <- dSSTT # Dispersal & attachment of SSTT (for each farm)
  
  M[seq(from=drowsRRTU, to=size, by=z), 
    seq(from=dcolsRRTU, to=size, by=z)] <- dRRTU 
  
  M[seq(from=drowsRSTU, to=size, by=z), 
    seq(from=dcolsRSTU, to=size, by=z)] <- dRSTU 
  
  M[seq(from=drowsSSTU, to=size, by=z), 
    seq(from=dcolsSSTU, to=size, by=z)] <- dSSTU 
  
  M[seq(from=drowsRRUU, to=size, by=z), 
    seq(from=dcolsRRUU, to=size, by=z)] <- dRRUU 
  
  M[seq(from=drowsRSUU, to=size, by=z), 
    seq(from=dcolsRSUU, to=size, by=z)] <- dRSUU 
  
  M[seq(from=drowsSSUU, to=size, by=z), 
    seq(from=dcolsSSUU, to=size, by=z)] <- dSSUU
  
  
  #### Fill M with fecundity (f) values ####
  f.t <- f.temp[, season.name] # Vector of reproductive output per adult per farm, for current season
  f.t <- matrix(rep(c(f.t), each=g*g), # Turn into matrix with a column for each genotype 
                ncol=g, nrow=i*g, byrow=T) # repeat rows 9 times
  
  # For each row (corresponding to genotype) adjust f according to resistance/tradeoffs and whether continuous strategy used on farms
  f.t[,1] <- f.t[,1] * rep(fRRTT.t, each=g) # Multiply f by proportion of f with/without continuous strategy for RRTT adults
  f.t[,2] <- f.t[,2] * rep(fRSTT.t, each=g) # Multiply f by proportion of f with/without continuous strategy for RSTT adults
  f.t[,3] <- f.t[,3] * rep(fSSTT.t, each=g)
  f.t[,4] <- f.t[,4] * rep(fRRTU.t, each=g)
  f.t[,5] <- f.t[,5] * rep(fRSTU.t, each=g)
  f.t[,6] <- f.t[,6] * rep(fSSTU.t, each=g)
  f.t[,7] <- f.t[,7] * rep(fRRUU.t, each=g)
  f.t[,8] <- f.t[,8] * rep(fRSUU.t, each=g)
  f.t[,9] <- f.t[,9] * rep(fSSUU.t, each=g)
   
  lr <- c(seq(1, size, 4)) # Create a vector corresponding to all larva rows (every 4th row, starting row 1) in M
  
  for(hh in 1:c(i*g)){ # For each larval row, 
    # Create vector where for each farm & genotype there are zeros except for:
    f.t.g <- c(rep(c(rep(0,3), f.t[hh,1], # RRTT fecundity in RR column
                     rep(0,3), f.t[hh,2], # RSTT fecundity in RS column
                     rep(0,3), f.t[hh,3], # SSTT fecundity in SS column
                     rep(0,3), f.t[hh,4],
                     rep(0,3), f.t[hh,5],
                     rep(0,3), f.t[hh,6],
                     rep(0,3), f.t[hh,7],
                     rep(0,3), f.t[hh,8],
                     rep(0,3), f.t[hh,9]),
                   times = i)) # Run this vector along entire row of M (Values will only appear in diagonal G matrices, since rest of row is 0)
    
    M[lr[hh],] <- M[lr[hh],] * f.t.g # Add fecundity values to larva rows of M
  } # Since rows in f.t are repeated 3 times, the same vector is applied for each genotype in a farm
  
  return(M) # Return completed M
}
