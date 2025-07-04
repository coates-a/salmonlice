#### Assign parameters to model ####

rm(list=ls()) # Clear workspace

args <- commandArgs(trailingOnly = TRUE) # Create command line for interacting with job_submission.slurm script
iter <- as.numeric(args[1]) # iter corresponds to array number (iteration) of job

# Check for package dependenices and install/load in packages as required.
#.libPaths("/home/alcoates/R/lib")
#lib = .libPaths()[1]
repo<- "https://cran.ms.unimelb.edu.au/" # Set mirror to download packages.
pkgs = c("readr", "dplyr", "gdata", "ggplot2", #"ggpubr",
        # "gganimate", "gifski", 
         "mapdata", 
        #"egg", "corrplot", "Matrix", 
        "magic", "reshape", 
        #"ggstance",
         "pillar") # Define packages.
#new.pkgs <- pkgs[!(pkgs %in% installed.packages(lib=lib)[,"Package"])]
#if(length(new.pkgs)) install.packages(new.pkgs, lib=lib, repos=repo)
inst = lapply(pkgs, library, character.only = TRUE)


#### Load datasets ####
farms <- read_csv("src/2loci/farms.week1.csv") # Farms info (farm size, starting lice numbers, lat & long)
paras <- read_csv("src/2loci/paras.csv") # Model parameters for batch of simulations
Dw <- read_csv("src/2loci/Dw.csv") # Dispersal matrix (winter)
Ds <- read_csv("src/2loci/Ds.csv") # Dispersal matrix (summer)
dw <- as.matrix(Dw)
ds <- as.matrix(Ds)
H <- read_csv("src/2loci/H.csv") # Distribution of continuous strategies across farms
H <- as.matrix(H) # Need to save as matrix to avoid tibbles
s.temp <- read_csv("src/2loci/temp.sw.csv") # Weekly growth rates per farm, per season (calculated from temperature)
s.temp <- as.matrix(s.temp)
rownames(s.temp)<-s.temp[,1] # Set farm IDs as rownames
s.temp <- s.temp[,-1] # Remove farm ID column
f.temp <- read_csv("src/2loci/temp.fw.csv") # Weekly reproductive rates per farm, per season (calculated from temperature)
f.temp <- as.matrix(f.temp)
rownames(f.temp)<-f.temp[,1] # Set farm IDs as rownames
f.temp <- f.temp[,-1] # Remove farm ID column


#### Assign parameters ####
i <- length(farms$farmID)   # Number of farms
g <- 9 # Number of genotypes: 2 loci = 9 total combinations
j <- 4 # Number of life stages
z <- g * j # Number of lice groups per farm (# genotypes * # life stages)
size <- i*z # Number of lice groups over ALL farms (= Size of M and P matrices)

list2env(x = as.list(paras[iter,]), # For the #iter row of the Parameters ('paras') dataset, assign each parameter a name corresponding to its column name
         envir = .GlobalEnv)

ts <- years * 52 # Number of time-steps in simulation (with 52 weekly time-steps per year) 

######

# Vector of forced harvest survival 
MAX <- rep(c(1, # No mortality for larvae (already dispersing from farm)
             rep(xmax, 3)), g) # Proportion survival from forced harvest same for all other stages and genotypes
# MAX will be multiplied by L matrix (all lice on farm) whenever forced harvest is applied, to remove 1-xmax proportion of lice


# Vector of weekly survival rate --> DIAGONAL of S matrix (lice REMAINING in same stage (1-s))
muD <- c(rep(c(1, muC, muP, muA),# Proportion survival after weekly mortality (mu) for each stage (C, P, A)
             times=g)) # Baseline mortality same for each genotype
# muD will be multiplied across diagonal of S matrix, to remove 1-mu proportion of non-transitioning lice (lost through background mortality)

# Vector of weekly survival rate --> OFF-DIAGONAL of S matrix (lice TRANSITIONING to next stage (s))
muO <- c(rep(c(1, 1, muC, muP),
             times=g))
# muO will be multiplied across diagonal of S matrix, to remove 1-mu proportion of lice (lost through background mortality)

# Calculate proportion larval infestation success (w = winter, s = summer)
dw <- dw * # Proportion of larvae that disperse TO each farm, FROM each farm *
  v # Proportion of those larvae that attach to host
ds <- ds * 
  v


#### ###
### GENOTYPE ADJUSTMENTS
# Determine how much survival parameters are further adjusted (by adding these values to baseline adjustments) according to genotype
# Effect of genotype in response to Strat X


#Effect of R allele (where relevant)

# Order of Matrix:
#RRTT
#RSTT
#SSTT,
#RRTU,
#RSTU,
#SSTU,
#RRTU,
#RSTU,
#SSTU

# Matching this order,
# Create vector of adjustments made at R/S locus
R.effect <- c(rep(c(RR.ad, RS.ad, SS.ad), times=3))

T.effect <- c(rep(c(TT.ad, TU.ad, UU.ad), each=3))
               


# TREATMENT survival from STRATEGY X
# Values will all be =1 if Strategy X is continuous, and won't be incorporated into simulation 

# Proportion survival from Strat X (BEFORE genotype adjustment)
x.stages <- ifelse(rep(stratx == "disc", times=z), # If strat x is Discrete
                   c(rep(c(1, xC, xP, xA), times=g)), # Assign baseline treatment survival from strategy X - all stages and genotypes
                   c(rep(1, z))) # Just 1s if strat x is Continuous

X <- x.stages + # Baseline treatment X survival,
  ifelse(rep(R.sel == "x", z), # PLUS any advantage from R allele
         rep(R.effect, each=j),
         rep(0, z)) +
  ifelse(rep(T.sel == "x", z), # PLUS any advantage from T allele
         rep(T.effect, each=j), # If T allele NOT under selection by X treatment,
         rep(0, z)) # ZERo additional fitness

X <- replace(X, X > 1, 1) # Cap proportion survival at maximum 1

# X will be multiplied across L (lice of all genotypes and stages in ONE farm) if Discrete X is applied

### REPEAT FOR STRATEGY Y ####

y.stages <- ifelse(rep(straty == "disc", times=z), # If strat x is Discrete
                   c(rep(c(1, yC, yP, yA), times=g)), # Assign baseline treatment survival from strategy X - all stages and genotypes
                   c(rep(1, z))) # Just 1s if strat x is Continuous

Y <- y.stages + # Baseline treatment Y survival,
  ifelse(rep(R.sel == "y", z), # PLUS any advantage from R allele
         rep(R.effect, each=j), 
         rep(0, z)) + # If R allele NOT under selection by Y treatment, ZERO additional fitness
  ifelse(rep(T.sel == "y", z), # PLUS any advantage from T allele
         rep(T.effect, each=j),
         rep(0, z)) 

Y <- replace(Y, Y > 1, 1) # Cap proportion survival at maximum 1

# Y will be multiplied across L (lice of all genotypes and stages in ONE farm) if Discrete Y is applied

### REPEAT FOR STRATEGY Z ####
# Include Z if X and Y are both continuous and impose selection on R and T
# Z always discrete, and imposes NO selection

Zdisc <- c(rep(c(1, zC, zP, zA), times=g)) # Assign baseline treatment survival from strategy Z - all stages and genotypes

# Z will be multiplied across L (lice of all genotypes and stages in ONE farm) if Discrete Z is applied



# Vector of the distribution of continuous strategies on all farms
H.vector <- H[, Hcol] # Assign which distribution (column of 'H') to use in #iter simulation
H.gen <- data.frame(H = c(rep(H.vector, each=g)), # Expand H so there is a row for each genotype 
                    geno = c(rep(c("RRTT","RSTT","SSTT",
                                   "RRTU","RSTU","SSTU",
                                   "RRUU","RSUU","SSUU"), times=i)))

H.vector.start <- H[, Hcol.start] # Assign which STARTING distribution (column of 'H') to use in #iter simulation
H.gen.start <- data.frame(H = c(rep(H.vector.start, each=g)), # Expand H so there is a row for each genotype 
                    geno = c(rep(c("RRTT","RSTT","SSTT",
                                   "RRTU","RSTU","SSTU",
                                   "RRUU","RSUU","SSUU"), times=i)))

# Depending on whether continuous strategy is present on farm or not, the following parameters are adjusted, according to genotype:
# development rate (s), chalimus survival (mu), larval attachment (L) or fecundity (f)
# For each parameter, create a vector of how they are adjusted (proportion relative to baseline), according to genotype and farm.

# Any advantage conferred by genotype to X (IF X is continuous)
x.cont.ad <- ifelse(rep(stratx == "cont" & R.sel == "x", g), # If Strat X is CONT AND selects for R allele,
                 R.effect, # Add any advantage to growth from R
                 rep(0, g)) + # ZERO additional fitness from R if strategy X isn't Cont and/or doesn't affect R
  ifelse(rep(stratx == "cont" & T.sel == "x", g), # AND If Strat X is CONT AND selects for T allele,
         T.effect, # Add any advantage from T
         rep(0, g)) 

# ADDITIONAL s under strategy Y, according to genotype
y.cont.ad <- ifelse(rep(straty == "cont" & R.sel == "y", g), 
                 R.effect, 
                 rep(0, g)) +
  ifelse(rep(straty == "cont" & T.sel == "y", g), 
         T.effect, 
         rep(0, g)) 

# Baseline reduction in v/f/s/mu from continuous strategies will be adjusted proportionally according to genotype by ADDING x.cont.ad and y.cont.ad 

# MAKE FOLLOWING A FUNCTION? - CAN PLUG IN FOR XS, YS, YV, ETC...

#############
#### ADJUSTING GROWTH RATE (s) ACCORDING TO STRATEGIES AND GENOTYPE

# Proportion reduction in GROWTH (for susceptible lice) caused by Strat X (if X is CONTINUOUS)
x.s <- ifelse(rep(stratx == "cont", # If strat x is Continuous
                  times=g), # (Need to repeat for each genotype to return a vector with length g)
                   c(rep(xs, # Assign baseline REDUCTION in s caused by strategy X
                         times=g)), # Repeat for each genotype
              c(rep(1, times=g))) # If X is discrete, no reduction in growth (x.s=1)


x.s.g <- x.s + x.cont.ad # Adjustment of s in response to X, PLUS any effecty of genotype
x.s.g <- replace(x.s.g, x.s.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect s (xs=1), OR if outsize effect of resistant alleles)


# The same for Strategy Y
# Proportion reduction in GROWTH (for susceptible lice) caused by Strat Y (if Y is CONTINUOUS)
y.s <- ifelse(rep(straty == "cont", 
                  times=g), 
              c(rep(ys, times=g)), 
              c(rep(1, times=g))) 

y.s.g <- y.s + y.cont.ad # Adjustment of s in response to Y, PLUS any effect of genotype
y.s.g <- replace(y.s.g, y.s.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect s (xs=1), OR if outsize effect of resistant alleles)


## TRADE-OFFS
s.R.to <- c(rep(c((2*sR.to), sR.to, 0), times=3)) # Trade-offs in GROWTH to R allele (sR.to) subtracted for EACH R allele
s.T.to <- c(rep(c((2*sT.to), sT.to, 0), each=3)) # Trade-offs in GROWTH to T allele (sT.to) subtracted for EACH T allele

# Vector of trade-offs for genotypes when strategy X NOT present (IF X is CONTINUOUS)
s.to.no.contx <- ifelse(rep(stratx == "cont", times=g), # If X is CONTINUOUS, there will be trade-offs when X NOT PRESENT
               ifelse(rep(R.sel == "x" & to.R == "absent", times=g), # if R selected for by Cont X,
                      s.R.to, # Any trade-offs from R allele
                      0) + # No trade-offs if R selected for by Y
                 ifelse(rep(T.sel == "x" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont X)
                        s.T.to,
                        0),
               0) # No trade-offs in absence of X if X is Discrete (calculated below)

# Vector of trade-offs for genotypes when strategy Y NOT present (IF Y is CONTINUOUS)
s.to.no.conty <- ifelse(rep(straty == "cont", times=g), # If Y is CONTINUOUS, there will be trade-offs when Y NOT PRESENT
       ifelse(rep(R.sel == "y" & to.R == "absent", times=g), # if R selected for by Cont Y,
              s.R.to, # Any trade-offs from R allele
              0) + # No trade-offs if R selected for by X
         ifelse(rep(T.sel == "y" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont Y)
                s.T.to,
                0),
       0)

# Vector of trade-offs to resistance to DISCRETE strategies (applied at ALL TIMES, regardless of whether used or not)
s.disc.to <- ifelse(rep((R.sel=="x" & stratx=="disc") | (R.sel=="y" & straty=="disc") | to.R == "always", times=g), # If R selected for by Discrete strategy
       s.R.to, # Apply any trade-offs of R
       rep(0,g)) +
  ifelse(rep((T.sel=="x" & stratx=="disc") | (T.sel=="y" & straty=="disc") | to.T == "always", times=g), # If T selected for by Discrete strategy
         s.T.to, # Apply any trade-offs of T
         rep(0,g))

# Vector calculating all growth adjustments
# X, no Y
s.x.noy <- x.s.g - # Effect of X on s (including any advantage of genotypes),
  s.to.no.conty - # MINUS any trade-offs from NO Y (if Y is CONT)
  s.disc.to # Minus any trade-offs from Discrete resistance
s.x.noy <- replace(s.x.noy, s.x.noy < 0, 0) # Ensure survival not <0 after trade-offs
 
# no X, Y
s.nox.y <- y.s.g - # Effect of Y on s for each genotype
  s.to.no.contx - # MINUS any trade-offs from NO X (if X is CONT)
  s.disc.to # Minus any trade-offs from Discrete resistance
s.nox.y <- replace(s.nox.y, s.nox.y < 0, 0)


# X, Y
s.x.y <- x.s.g * # Adjustment to s from X on genotypes
  y.s.g - # MULTIPLIED BY adjustments from Y 
  s.disc.to # Minus any trade-offs from Discrete resistance
s.x.y <- replace(s.x.y, s.x.y < 0, 0)

s.nox.noy <- 1 - s.to.no.contx - s.to.no.conty - s.disc.to # No X or Y (just tradeoffs)

##### Match vector of s (growth rate) adjustments for genotypes to each farm (depending on which strategies are used)

H.gen$s.cont <- ifelse(H.gen$H=="x" | H.gen$H=="xz",
                       s.x.noy,
                       ifelse(H.gen$H=="y" | H.gen$H=="yz",
                              s.nox.y,
                              ifelse(H.gen$H=="xy" | H.gen$H=="xyz",
                                     s.x.y,
                                     s.nox.noy)))

# Repeat for H.gen.start (Starting distribution)
H.gen.start$s.cont <- ifelse(H.gen.start$H=="x" | H.gen.start$H=="xz",
                       s.x.noy,
                       ifelse(H.gen.start$H=="y" | H.gen.start$H=="yz",
                              s.nox.y,
                              ifelse(H.gen.start$H=="xy" | H.gen.start$H=="xyz",
                                     s.x.y,
                                     s.nox.noy)))


# s.cont will be multiplied by s (baseline growth rate on each farm) to account for any variation in growth due to genotypes (resistance & trade-offs) and/or continuous strategies

##################### REPEAT FOR mu ###################


# Proportion reduction in (Chalimus) BACKGROUND MORTALITY SURVIVAL (for susceptible lice) caused by Strat X (if X is CONTINUOUS)
x.mu <- ifelse(rep(stratx == "cont", # If strat x is Continuous
                  times=g), # (Need to repeat for each genotype to return a vector with length g)
              c(rep(xmu, # Assign baseline REDUCTION in mu caused by strategy X
                    times=g)), # Repeat for each genotype
              c(rep(1, times=g))) # If X is discrete, no reduction in survival (x.s=1)

x.mu.g <- x.mu + x.cont.ad # Adjustment of mu in response to X, for each genotype
x.mu.g <- replace(x.mu.g, x.mu.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect mu (xmu=1), OR if outsize effect of resistant alleles)

# The same for Strategy Y
# Proportion reduction in SURVIVAL (for susceptible lice) caused by Strat Y (if Y is CONTINUOUS)
y.mu <- ifelse(rep(straty == "cont", 
                  times=g), 
              c(rep(ymu, times=g)), 
              c(rep(1, times=g))) 


y.mu.g <- y.mu + y.cont.ad # Adjustment of mu in response to X, for each genotype
y.mu.g <- replace(y.mu.g, y.mu.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect mu (xmu=1), OR if outsize effect of resistant alleles)


## TRADE-OFFS
mu.R.to <- c(rep(c((2*muR.to), muR.to, 0), times=3)) # Trade-offs in GROWTH to R allele (sR.to) subtracted for EACH R allele
mu.T.to <- c(rep(c((2*muT.to), muT.to, 0), each=3)) # Trade-offs in GROWTH to T allele (sT.to) subtracted for EACH T allele

# Vector of trade-offs for genotypes when strategy X NOT present (IF X is CONTINUOUS)
mu.to.no.contx <- ifelse(rep(stratx == "cont", times=g), # If X is CONTINUOUS, there will be trade-offs when X NOT PRESENT
                         ifelse(rep(R.sel == "x" & to.R == "absent", times=g), # if R selected for by Cont X,
                                mu.R.to, # Any trade-offs from R allele
                                0) + # No trade-offs if R selected for by Y
                           ifelse(rep(T.sel == "x" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont X)
                                  mu.T.to,
                                  0),
                         0) # No trade-offs in absence of X if X is Discrete (calculated below)

# Vector of trade-offs for genotypes when strategy Y NOT present (IF Y is CONTINUOUS)
mu.to.no.conty <- ifelse(rep(straty == "cont", times=g), # If Y is CONTINUOUS, there will be trade-offs when Y NOT PRESENT
                         ifelse(rep(R.sel == "y" & to.R == "absent", times=g), # if R selected for by Cont Y,
                                mu.R.to, # Any trade-offs from R allele
                                0) + # No trade-offs if R selected for by X
                           ifelse(rep(T.sel == "y" & to.R == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont Y)
                                  mu.T.to,
                                  0),
                         0)

# Vector of trade-offs to resistance to DISCRETE strategies (applied at ALL TIMES, regardless of whether used or not)
mu.disc.to <- ifelse(rep((R.sel=="x" & stratx=="disc") | (R.sel=="y" & straty=="disc") | to.R == "always", times=g), # If R selected for by Discrete strategy
                     mu.R.to, # Apply any trade-offs of R
                     rep(0,g)) +
  ifelse(rep((T.sel=="x" & stratx=="disc") | (T.sel=="y" & straty=="disc") | to.T == "always", times=g), # If T selected for by Discrete strategy
         mu.T.to, # Apply any trade-offs of T
         rep(0,g))

# Vector calculating all mortality adjustments
# X, no Y
mu.x.noy <- x.mu.g - # Effect of strategy X on mu, plus any resistance
  mu.to.no.conty - # MINUS any trade-offs from NO Y (if Y is CONT)
  mu.disc.to # Minus any trade-offs from Discrete resistance
mu.x.noy <- replace(mu.x.noy, mu.x.noy < 0, 0) # Ensure survival after trade-offs isn't <0

# no X, Y
mu.nox.y <- y.mu.g - # PLUS any advantage of genotypes,
  mu.to.no.contx - # MINUS any trade-offs from NO X (if X is CONT)
  mu.disc.to # Minus any trade-offs from Discrete resistance
mu.nox.y <- replace(mu.nox.y, mu.nox.y < 0, 0)


# X, Y
mu.x.y <- x.mu.g * # Effect of X on mu (including resistance) 
  y.mu.g - # MULTIPLIED BY effect of strategy Y, 
  mu.disc.to # Minus any trade-offs from Discrete resistance
mu.x.y <- replace(mu.x.y, mu.x.y < 0, 0)

mu.nox.noy <- 1 - mu.to.no.contx - mu.to.no.conty - mu.disc.to

##### Match vector of s (growth rate) adjustments for genotypes to each farm (depending on which strategies are used)

H.gen$mu.cont <- ifelse(H.gen$H=="x" | H.gen$H=="xz",
                       mu.x.noy,
                       ifelse(H.gen$H=="y" | H.gen$H=="yz",
                              mu.nox.y,
                              ifelse(H.gen$H=="xy" | H.gen$H=="xyz",
                                     mu.x.y,
                                     mu.nox.noy))) # 

# mu.cont will be multiplied across M to adjust chalimus survival rate (mu) based on genotype and farm

H.gen.start$mu.cont <- ifelse(H.gen.start$H=="x" | H.gen.start$H=="xz",
                        mu.x.noy,
                        ifelse(H.gen.start$H=="y" | H.gen.start$H=="yz",
                               mu.nox.y,
                               ifelse(H.gen.start$H=="xy" | H.gen.start$H=="xyz",
                                      mu.x.y,
                                      mu.nox.noy))) # 


######################################################
# For each genotype, create vector of proportion of baseline larval attachment (v) for each farm, according to whether continuous strategy is used

# Baseline ADJUSTMENT to v (according to strategy, but NOT genotype)


# Proportion reduction in ATTACHMENT (for susceptible lice) caused by Strat X (if X is CONTINUOUS)
x.v <- ifelse(rep(stratx == "cont", # If strat x is Continuous
                   times=g), # (Need to repeat for each genotype to return a vector with length g)
               c(rep(xv, # Assign baseline REDUCTION in v caused by strategy X
                     times=g)), # Repeat for each genotype
               c(rep(1, times=g))) # If X is discrete, no reduction in survival (x.s=1)

x.v.g <- x.v + x.cont.ad # Adjustment of v in response to X, for each genotype
x.v.g <- replace(x.v.g, x.v.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect v (xv=1), OR if outsize effect of resistant alleles)

# The same for Strategy Y
# Proportion reduction in SURVIVAL (for susceptible lice) caused by Strat Y (if Y is CONTINUOUS)
y.v <- ifelse(rep(straty == "cont", 
                   times=g), 
               c(rep(yv, times=g)), 
               c(rep(1, times=g))) 

y.v.g <- y.v + y.cont.ad # Adjustment of mu in response to X, for each genotype
y.v.g <- replace(y.v.g, y.v.g > 1, 1) # Ensure survival not >1 (if X doesn't actually affect mu (xmu=1), OR if outsize effect of resistant alleles)


## TRADE-OFFS
v.R.to <- c(rep(c((2*vR.to), vR.to, 0), times=3)) # Trade-offs in GROWTH to R allele (sR.to) subtracted for EACH R allele
v.T.to <- c(rep(c((2*vT.to), vT.to, 0), each=3)) # Trade-offs in GROWTH to T allele (sT.to) subtracted for EACH T allele

# Vector of trade-offs for genotypes when strategy X NOT present (IF X is CONTINUOUS)
v.to.no.contx <- ifelse(rep(stratx == "cont", times=g), # If X is CONTINUOUS, there will be trade-offs when X NOT PRESENT
                         ifelse(rep(R.sel == "x" & to.R == "absent", times=g), # if R selected for by Cont X,
                                v.R.to, # Any trade-offs from R allele
                                0) + # No trade-offs if R selected for by Y
                           ifelse(rep(T.sel == "x" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont X)
                                  v.T.to,
                                  0),
                         0) # No trade-offs in absence of X if X is Discrete (calculated below)

# Vector of trade-offs for genotypes when strategy Y NOT present (IF Y is CONTINUOUS)
v.to.no.conty <- ifelse(rep(straty == "cont", times=g), # If Y is CONTINUOUS, there will be trade-offs when Y NOT PRESENT
                         ifelse(rep(R.sel == "y" & to.R == "absent", times=g), # if R selected for by Cont Y,
                                v.R.to, # Any trade-offs from R allele
                                0) + # No trade-offs if R selected for by X
                           ifelse(rep(T.sel == "y" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont Y)
                                  v.T.to,
                                  0),
                         0)

# Vector of trade-offs to resistance to DISCRETE strategies (applied at ALL TIMES, regardless of whether used or not)
v.disc.to <- ifelse(rep((R.sel=="x" & stratx=="disc") | (R.sel=="y" & straty=="disc") | to.R == "always", times=g), # If R selected for by Discrete strategy
                     v.R.to, # Apply any trade-offs of R
                     rep(0,g)) +
  ifelse(rep((T.sel=="x" & stratx=="disc") | (T.sel=="y" & straty=="disc") | to.T == "always", times=g), # If T selected for by Discrete strategy
         v.T.to, # Apply any trade-offs of T
         rep(0,g))

# Vector calculating all mortality adjustments
# X, no Y
v.x.noy <- x.v.g - # Effect of strategy X on v, plus any resistance
  v.to.no.conty - # MINUS any trade-offs from NO Y (if Y is CONT)
  v.disc.to # Minus any trade-offs from Discrete resistance
v.x.noy <- replace(v.x.noy, v.x.noy < 0, 0) # Ensure survival after trade-offs isn't <0

# no X, Y
v.nox.y <- y.v.g - # PLUS any advantage of genotypes,
  v.to.no.contx - # MINUS any trade-offs from NO X (if X is CONT)
  v.disc.to # Minus any trade-offs from Discrete resistance
v.nox.y <- replace(v.nox.y, v.nox.y < 0, 0)

# X, Y
v.x.y <- x.v.g * # Effect of X on v (including resistance) 
  y.v.g - # vLTIPLIED BY effect of strategy Y, 
  v.disc.to # Minus any trade-offs from Discrete resistance
v.x.y <- replace(v.x.y, v.x.y < 0, 0)

v.nox.noy <- 1 - v.to.no.contx - v.to.no.conty - v.disc.to # No X or Y (just tradeoffs)

H.gen$v.cont <- ifelse(H.gen$H=="x" | H.gen$H=="xz",
                        v.x.noy,
                        ifelse(H.gen$H=="y" | H.gen$H=="yz",
                               v.nox.y,
                               ifelse(H.gen$H=="xy" | H.gen$H=="xyz",
                                      v.x.y,
                                      v.nox.noy)))

# For starting distribution
H.gen.start$v.cont <- ifelse(H.gen.start$H=="x" | H.gen.start$H=="xz",
                       v.x.noy,
                       ifelse(H.gen.start$H=="y" | H.gen.start$H=="yz",
                              v.nox.y,
                              ifelse(H.gen.start$H=="xy" | H.gen.start$H=="xyz",
                                     v.x.y,
                                     v.nox.noy)))

       
# Create separate matrix vectors for each genotype attachment success, for EVERY FARM (according to strategies used)

vRRTT.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RRTT")),
                  ncol=1, nrow=i)
vRSTT.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RSTT")),
                  ncol=1, nrow=i)
vSSTT.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="SSTT")),
                  ncol=1, nrow=i)
vRRTU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RRTU")),
                  ncol=1, nrow=i)
vRSTU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RSTU")),
                  ncol=1, nrow=i)
vSSTU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="SSTU")),
                  ncol=1, nrow=i)
vRRUU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RRUU")),
                  ncol=1, nrow=i)
vRSUU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="RSUU")),
                  ncol=1, nrow=i)
vSSUU.h <- matrix(c(subset(H.gen$v.cont, H.gen$geno=="SSUU")),
                  ncol=1, nrow=i)
# v_.h vectors will be multiplied across M to adjust attachment by v_ based on genotype and farm

# Repeat for Starting distribution

vRRTT.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RRTT")),
                  ncol=1, nrow=i)
vRSTT.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RSTT")),
                  ncol=1, nrow=i)
vSSTT.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="SSTT")),
                  ncol=1, nrow=i)
vRRTU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RRTU")),
                  ncol=1, nrow=i)
vRSTU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RSTU")),
                  ncol=1, nrow=i)
vSSTU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="SSTU")),
                  ncol=1, nrow=i)
vRRUU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RRUU")),
                  ncol=1, nrow=i)
vRSUU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="RSUU")),
                  ncol=1, nrow=i)
vSSUU.h.start <- matrix(c(subset(H.gen.start$v.cont, H.gen.start$geno=="SSUU")),
                  ncol=1, nrow=i)

##########################

# Baseline ADJUSTMENT to f (according to strategy, but NOT genotype)


# Proportion reduction in FECUNDITY (for susceptible lice) caused by Strat X (if X is CONTINUOUS)
x.f <- ifelse(rep(stratx == "cont", # If strat x is Continuous
                  times=g), # (Need to repeat for each genotype to return a vector with length g)
              c(rep(xf, # Assign baseline REDUCTION in f caused by strategy X
                    times=g)), # Repeat for each genotype
              c(rep(1, times=g))) # If X is discrete, no reduction in f (x.f=1)

x.f.g <- x.f + x.cont.ad # Adjustment of f in response to X, for each genotype
x.f.g <- replace(x.f.g, x.f.g > 1, 1) # Ensure not >1 (if X doesn't actually affect f (xf=1), OR if outsize effect of resistant alleles)

# The same for Strategy Y
# Proportion reduction in FECUNDITY (for susceptible lice) caused by Strat Y (if Y is CONTINUOUS)
y.f <- ifelse(rep(straty == "cont", 
                  times=g), 
              c(rep(yf, times=g)), 
              c(rep(1, times=g))) 

y.f.g <- y.f + y.cont.ad # Adjustment of f in response to X, for each genotype
y.f.g <- replace(y.f.g, y.f.g > 1, 1) # Ensure f adjustment not >1


## TRADE-OFFS
f.R.to <- c(rep(c((2*fR.to), fR.to, 0), times=3)) # Trade-offs in GROWTH to R allele (sR.to) subtracted for EACH R allele
f.T.to <- c(rep(c((2*fT.to), fT.to, 0), each=3)) # Trade-offs in GROWTH to T allele (sT.to) subtracted for EACH T allele

# Vector of trade-offs for genotypes when strategy X NOT present (IF X is CONTINUOUS)
f.to.no.contx <- ifelse(rep(stratx == "cont", times=g), # If X is CONTINUOUS, there will be trade-offs when X NOT PRESENT
                        ifelse(rep(R.sel == "x" & to.R == "absent", times=g), # if R selected for by Cont X,
                               f.R.to, # Any trade-offs from R allele
                               0) + # No trade-offs if R selected for by Y
                          ifelse(rep(T.sel == "x" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont X)
                                 f.T.to,
                                 0),
                        0) # No trade-offs in absence of X if X is Discrete (calculated below)

# Vector of trade-offs for genotypes when strategy Y NOT present (IF Y is CONTINUOUS)
f.to.no.conty <- ifelse(rep(straty == "cont", times=g), # If Y is CONTINUOUS, there will be trade-offs when Y NOT PRESENT
                        ifelse(rep(R.sel == "y" & to.R == "absent", times=g), # if R selected for by Cont Y,
                               f.R.to, # Any trade-offs from R allele
                               0) + # No trade-offs if R selected for by X
                          ifelse(rep(T.sel == "y" & to.T == "absent", times=g), # PLUS any trade-offs from T allele (if T selected for by Cont Y)
                                 f.T.to,
                                 0),
                        0)

# Vector of trade-offs to resistance to DISCRETE strategies (applied at ALL TIMES, regardless of whether used or not)
f.disc.to <- ifelse(rep((R.sel=="x" & stratx=="disc") | (R.sel=="y" & straty=="disc") | to.R == "always", times=g), # If R selected for by Discrete strategy
                    f.R.to, # Apply any trade-offs of R
                    rep(0,g)) +
  ifelse(rep((T.sel=="x" & stratx=="disc") | (T.sel=="y" & straty=="disc") | to.T == "always", times=g), # If T selected for by Discrete strategy
         f.T.to, # Apply any trade-offs of T
         rep(0,g))

# Vector calculating all mortality adjustments
# X, no Y
f.x.noy <- x.f.g - # Effect of strategy X on f, plus any resistance
  f.to.no.conty - # MINUS any trade-offs from NO Y (if Y is CONT)
  f.disc.to # Minus any trade-offs from Discrete resistance
f.x.noy <- replace(f.x.noy, f.x.noy < 0, 0) # Ensure surfifal after trade-offs isn't <0

# no X, Y
f.nox.y <- y.f.g - # PLUS any adfantage of genotypes,
  f.to.no.contx - # MINUS any trade-offs from NO X (if X is CONT)
  f.disc.to # Minus any trade-offs from Discrete resistance
f.nox.y <- replace(f.nox.y, f.nox.y < 0, 0)

# X, Y
f.x.y <- x.f.g * # Effect of X on f (including resistance) 
  y.f.g - # MULTIPLIED BY effect of strategy Y, 
  f.disc.to # Minus any trade-offs from Discrete resistance
f.x.y <- replace(f.x.y, f.x.y < 0, 0)

f.nox.noy <- 1 - f.to.no.contx - f.to.no.conty - f.disc.to # No X or Y (just tradeoffs)

H.gen$f.cont <- ifelse(H.gen$H=="x" | H.gen$H=="xz",
                       f.x.noy,
                       ifelse(H.gen$H=="y" | H.gen$H=="yz",
                              f.nox.y,
                              ifelse(H.gen$H=="xy" | H.gen$H=="xyz",
                                     f.x.y,
                                     f.nox.noy)))
# For starting distribution
H.gen.start$f.cont <- ifelse(H.gen.start$H=="x" | H.gen.start$H=="xz",
                       f.x.noy,
                       ifelse(H.gen.start$H=="y" | H.gen.start$H=="yz",
                              f.nox.y,
                              ifelse(H.gen.start$H=="xy" | H.gen.start$H=="xyz",
                                     f.x.y,
                                     f.nox.noy)))


# Create separate matrix vectors for each genotype attachment success, for EVERY FARM (according to strategies used)

fRRTT.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RRTT")),
                  ncol=1, nrow=i)
fRSTT.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RSTT")),
                  ncol=1, nrow=i)
fSSTT.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="SSTT")),
                  ncol=1, nrow=i)
fRRTU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RRTU")),
                  ncol=1, nrow=i)
fRSTU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RSTU")),
                  ncol=1, nrow=i)
fSSTU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="SSTU")),
                  ncol=1, nrow=i)
fRRUU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RRUU")),
                  ncol=1, nrow=i)
fRSUU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="RSUU")),
                  ncol=1, nrow=i)
fSSUU.h <- matrix(c(subset(H.gen$f.cont, H.gen$geno=="SSUU")),
                  ncol=1, nrow=i)

# f_.h vectors will be multiplied across M to adjust fecundity by f_ based on genotype and farm

# Repeat for STARTING distribution
fRRTT.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RRTT")),
                  ncol=1, nrow=i)
fRSTT.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RSTT")),
                  ncol=1, nrow=i)
fSSTT.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="SSTT")),
                  ncol=1, nrow=i)
fRRTU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RRTU")),
                  ncol=1, nrow=i)
fRSTU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RSTU")),
                  ncol=1, nrow=i)
fSSTU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="SSTU")),
                  ncol=1, nrow=i)
fRRUU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RRUU")),
                  ncol=1, nrow=i)
fRSUU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="RSUU")),
                  ncol=1, nrow=i)
fSSUU.h.start <- matrix(c(subset(H.gen.start$f.cont, H.gen.start$geno=="SSUU")),
                  ncol=1, nrow=i)

# Matrix of production cycles of each farm
prod.cycle <- read_csv("src/2loci/prod.cycle.csv")
prod.cycle <- as.matrix(prod.cycle)
# Subset prod.cycle so same dimensions as N
prod.cycle <- prod.cycle[, c(105: # Start at beginning of THIRD year (time for all farms to start up)
                       (104+ts))] # Run for the number of time-steps


#### Call script to run Simulation ####

source("src/2loci/popSims_2loci.R") # Pull parameters & variables to use


