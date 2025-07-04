# Create data-frame template for model parameters

paras_gm <- data.frame(scenarios = NA, # Name of batch of simulations
                    sim = NA, # Name of simulation
                    years = NA, # Number of years to run
                    Hcol = NA, # Which column of 'H.csv' to use. In each column, every farm is assigned one or more anti-lice treatments (given as 'x', 'y' and 'z'))
                    years.start = NA, # Number of years 'spin-up' (before Hcol used)
                    Hcol.start = NA, # Which column of 'H.csv' to use in the time BEFORE years.start
                    
                    # Lice limits (abundance of adults per fish) that trigger:
                    kws = NA, # treatment (winter & summer)
                    ksp = NA, # treatment (spring)
                    kmax = NA, # forced harvest
                    xmax = NA, # Proportion survival during forced harvest
        
                    # Proportion background mortality survival (1 - weekly mortality rate)
                    muA = NA, # Adults
                    muP = NA, # Pre-adults
                    muC = NA, # Chalimus
                    
                    v = NA, # Proportion attachment success of Larvae (after successful dispersal)
                    
                    # Parameters for Control strategy 'X'
                    stratx = NA, # either 'cont' (continuous strategy, e.g. lice skirt) or 'disc' (discrete strategy, e.g. chemical delousing)
                    # These values are equivalent to a_mu and b_mu in Coate et al., (2025)
                    # They can take the value between 0 and 1 and are multiplied against the parameter of interest
                    # E.g. if larval attachment with Treatment X is 25% that of normal levels, then xv = 0.25 (attachment will be adjusted v*xv)
                    xA = NA, # Proportion adults surviving x (if discrete)
                    xP = NA, # Proportion pre-adults surviving x (if discrete)
                    xC = NA, # Proportion chalimus surviving x (if discrete)
                    xv = NA, # Proportion of baseline larval attachment under x (if continuous)
                    xf = NA, # Proportion of baseline fecundity under x (if continuous)
                    xs = NA, # Proportion of baseline growth rate under x (if continuous)
                    xmu = NA, # Proportion of baseline chalimus survival under x (if continuous)
                    
                    # Same as above, for control strategy 'Y'
                    straty = NA, #
                    yA = NA, #
                    yP = NA, # 
                    yC = NA, # 
                    yv = NA, # 
                    yf = NA, # 
                    ys = NA, # 
                    ymu = NA, #
                    
                    # Fitness advantages conferred by R and T under selection:
                    R.sel = NA, # Strategy ('x' or 'y') selecting for R allele
                    T.sel = NA, # Strategy ('x' or 'y') selecting for T allele
                    # These values are equivalent to alpha_R and alpha_T in Coates et al. (2025)
                    # Take the value between 0 and 1 and are ADDED to the parameter affected by the Treatment 
                    # E.g. If survival of SUSCEPTIBLE adults from Treat Y is 10% (yA = 0.1) but survival of RR adults is 60%, then:
                    # RR.ad = 0.5 (as yA + RR.ad = 0.1 + 0.5 = 0.6)
                    RR.ad = NA, # Advantage conferred to RR genotypes under treatment 
                    RS.ad = NA, # Advantage conferred to RS genotypes
                    SS.ad = NA, # Advantage conferred to SS genotypes
                    TT.ad = NA, # Advantage conferred to TT genotypes under treatment 
                    TU.ad = NA, # Advantage conferred to TU genotypes
                    UU.ad = NA, # Advantage conferred to UU genotypes
                    
                    # Fitness trade-offs associated with R and T alleles
                    # Take the value between 0 and 1 and are SUBTRACTED (for each copy of the allele) from the parameter of interest (thus reducing attachment, fecundity, etc)
                    vR.to = NA, # Trade-off to larval attachment from R
                    fR.to = NA, # Trade-off to fecundity from R
                    sR.to = NA, # Trade-off to development from R
                    muR.to = NA, # Trade-off to chalimus survival from R
                    vT.to = NA, # Same as above for T allele
                    fT.to = NA,
                    sT.to = NA,
                    muT.to = NA,
                    to.R = NA, # Assign whether R trade-offs are applied at all times ('always') or only when main selection pressure is absent ('absent')
                    to.T = NA, # Assign whether T trade-offs are applied at all times ('always') or only when main selection pressure is absent ('absent')
              

                    prop.RSUU = NA, # Proportion of RSUU lice on farms at start of simulation
                    prop.SSTU = NA, # Proportion of SSTU lice on farms at start of simulation
                    prop.RSTU = NA, # Proportion of RSTU lice on farms at start of simulation
                    
                    # Can assign a lower and upper latitude limit to where the R and T alleles are first present at the start of the sim
                    lwr.latR = NA, # 
                    upr.latR = NA, # 
                    lwr.latT = NA, # 
                    upr.latT = NA, # 
                    
                    zA = NA, # Proportion adult survival after Treatment Z (which is always discrete)
                    zA = NA, # Proportion pre-adult survival after Treatment Z (which is always discrete)
                    zA = NA, # Proportion chalimus survival after Treatment Z (which is always discrete)
                    
                    # Explore different starting lice levels
                    lice.start = NA # Takes value 1 - 4 (see popSims)
                    # 1 = levels from Barentswatch, 2 = levels following 2-year spin-up,
                    # 3 = 0.2 lice/stage/fish on all farms, 4 = 0.1 lice/stage/fish on all farms
)

# Save as CSV file to be filled parameter values
write.csv(paras, file = "src/aza_temp/paras_aza_temp.csv", row.names = F)
