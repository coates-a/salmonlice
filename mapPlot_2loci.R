################## PLOTS ######################
############ MAPS #############


########## NEED TO DO MAPS ######

# Adding map
norwaymap <- map_data("worldHires", "Norway", xlim = c(4,11), ylim = c(58,66.5)) # Save map of Norway (within lat & long of sim)

### ANIMATION OF P (FREQUENCY OF R) THROUGH TIME & SPACE 

# Looking just at adults
adrr <- subset(popdyn, stage=="ad" & geno=="RRTT" & abund.stage >=0.01) # p is repeated for each genotype (so only need to look at 1), ignore farms with <0.01 abund (avoid NAs)
options(bitmapType='cairo') # To save png correctly

adrr.p <- melt(adrr, # Melt p for each locus (pT, pR) into single column
               id.vars = colnames(select(adrr, -c(pT, pR)))) # Keep all other columns the same  
data.table::setnames(adrr.p,'variable','p.locus') # New column names
data.table::setnames(adrr.p,'value','p')


map <- ggplot()+
  geom_map(data = norwaymap, aes(map_id = region), # Draw map
           map = norwaymap, colour = NA, fill = "grey85")+
  geom_point(data=adrr.p, # Add farm sites
             aes(x=lon, y=lat, 
                 size=abund.stage, 
                 col=p,
                 shape = H))+ # Shape of point indicates host type at farm
  
  ##  ** For now, just using H AFTER new distribution is added (ie ignoring H.start). Need to change if I want to plot new distributions through time
  
  facet_wrap(. ~ p.locus, # Facet plots by R/S and T/U loci
             ncol=2,
             labeller = as_labeller(c(`pR` = "R frequency", `pT` = "T frequency")))+
  scale_color_gradientn(name="Allele frequency", # Colour gradient based on p
                        limits=c(0,1), 
                        colors = c("#00CCFF","#FF9933","#FF0000"), # Set a gradient along these colours (blue to red)
                        breaks=c(0,0.5,1))+
  scale_size(name="Abundance\n(adult lice/fish)", # Scale point size according to number of lice
             limits=c(0.01, 10), # Only plot farms with at least 1 adult per 100 fish
             trans="log10", # Log scale
             range=c(0.75,2.75))+ # max and min size of points
  scale_shape_manual(name="Strategies",
                     values=c(17, 16, 15, 18, 2, 1, 0, 6))+ # Up to 7 treatment combos
  scale_x_continuous(name = expression(Longitude~(degree~E)))+
  scale_y_continuous(name = expression(Latitude~(degree~N)))+
  coord_fixed(1.8)+ # Fixed aspect ratio (so map not distorted)
  geom_hline(yintercept = c(lwr.latR, upr.latR),
             alpha = ifelse(lwr.latR > 58 | upr.latR < 67, # If it's a simulation where R allele starts within subset of study area, 
                            1, # add lines marking the upper and lower starting range
                            0), # Otherwise, lines are transparent
             lty=2)+
  geom_hline(yintercept = c(lwr.latT, upr.latT),
             alpha = ifelse(lwr.latT > 58 | upr.latT < 67, # If it's a simulation where T allele starts within subset of study area, 
                            1, # add lines marking the upper and lower starting range
                            0), # Otherwise, lines are transparent
             lty=3)+ # Different line than for R allele
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key= element_rect(fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.1))

adrr.p$tt <- round(adrr.p$tt,2) # Round tt to 3 decimal places (clearer labels)

anim <- map + 
  transition_states(adrr.p$tt) + # Different frame for each time point
  ggtitle("Simulation: '{sim}' \n Year = {closest_state}") # Label by simulation (parameter combo) and nearest tt of that frame

anim <- animate(anim, # Animate frames
                renderer  = gifski_renderer(), 
                nframes = (2*ts)+30, # Animate seems to use 2 frames for each 1 tt...
                end_pause = 30, # Freeze last frame for 30 frames
                fps = 20,
                width=800, height=720) # try and keep 1:1.8 aspect

anim_save(file = sprintf("outputs/2loci/plots/animRxT_%s_%s.gif", scenarios, sim), anim) # Save


######## SUMMARY PLOTS ######
options(bitmapType='cairo') # To save png correctly

### CALCULATE N and GENE FREQUENCY (p) IN WHOLE METAPOPULATION ####

# Look just at adults for now
ads <- subset(popdyn, stage=="ad") 

gf <- aggregate(N ~ geno + genoR + genoT + tt, data=ads, FUN = sum) # Sum all adults of each genotype across all farms
Ntot <- aggregate(N ~ tt, data=ads, FUN=sum) # Total number of adults in sim
data.table::setnames(Ntot,'N','Ntot')
gf <- merge(x=gf, y=Ntot, by="tt", all.y=F)
gf$geno.prop <- gf$N / gf$Ntot # Proportion of genotypes in whole metapopulation
# Number of R alleles in whole adult metapopulation
gf$nR <- gf$N * ifelse(gf$genoR == "RR", # For each louse,
                           2, # 2x R alleles in RR
                           ifelse(gf$genoR == "RS", 
                                  1, # 1x in RS
                                  0)) # 0 in SS
nRtot <- aggregate(nR ~ tt, data=gf, FUN=sum)
data.table::setnames(nRtot,'nR','nRtot')
gf <- merge(x=gf, y=nRtot, by="tt")
gf$pR <- gf$nRtot / (2* gf$Ntot) # p = number R alleles / 2*number adults
# Number of T alleles in whole adult metapopulation
gf$nT <- gf$N * ifelse(gf$genoT == "TT", # For each louse,
                       2, # 2x T alleles in TT
                       ifelse(gf$genoT == "TU", 
                              1, # 1x in TU
                              0)) # 0 in UU
nTtot <- aggregate(nT ~ tt, data=gf, FUN=sum)
data.table::setnames(nTtot,'nT','nTtot')
gf <- merge(x=gf, y=nTtot, by="tt")
gf$pT <- gf$nTtot / (2* gf$Ntot) # p = number T alleles / 2*number adults


gf.p <- subset(gf, geno=="SSUU") # Just need to plot for 1 genotype, since p values are repeated for each genotype
gf.p <- gf.p[,c("tt","pR","pT")]
gf.p <- melt(gf.p, # Melt data into single column
                    id.vars = "tt")
colnames(gf.p)[c(2,3)]<-c('locus', 'p')

#### Plot OVERALL p through time ####
p.plot <- ggplot(data = gf.p, 
 aes(x=tt, y=p, col=locus))+
  geom_line(lwd=2)+
  scale_color_manual(name="Allele",
                     label = c('R','T'),
                     values=c("#FF0000", "#00CCFF"))+
  scale_x_continuous(name="Year", limits=c(0,years), breaks=seq(0,years,by=1))+
  scale_y_continuous(name=expression("Frequency of allele ("*italic(p)*")"), 
                     limits=c(0,1), breaks=seq(0,1,by=0.25))+
  ggtitle("A)")+
  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=20, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "gray50", linewidth = 0.6, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA,linewidth=1.5))

p.plot


### Plot OVERALL N of each genotype ####
genoNplot <- ggplot(gf, 
                    aes(x=tt, y=N, col=geno))+
  geom_line(lwd=1)+
  scale_x_continuous(name="Year", limits=c(0,years+1), breaks=seq(0,years,by=1))+
  scale_y_continuous(name="Total number of adults")+
  #scale_color_manual(name="Genotype",
  #                   values=c("#FF0000","#FF9933","#00CCFF"))+
  ggtitle("B)")+
  geom_text(data=subset(gf,tt==max(gf$tt)), # Only label final point
            aes(x=tt, y=N, label=geno, col=geno),
            #position=position_dodgev(height=1000),
            vjust=-0.5,
            hjust=-0.01)+
  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=20, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 14),
        legend.key= element_rect(fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "gray50", linewidth = 0.6, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA,linewidth=1.5))

genoNplot


##### Calculate MEAN ABUNDACE across all farms #####

m.ab <- subset(popdyn, geno=="SSUU" & stage=="ad" &
                 prodcy==1) # Only average ACTIVE farms
m.ab <- aggregate(abund.stage ~ tt, data=m.ab, FUN = mean)

# Plot OVERALL mean abundance over time
ab.plot <- ggplot(m.ab, 
                  aes(x=tt, y=abund.stage))+
  geom_line(lwd=2, colour="red")+
  scale_x_continuous(name="Year", limits=c(0,years), breaks=seq(0,years,by=1))+
  scale_y_continuous(name=expression("Mean adult abundance (lice/fish)"))+
  ggtitle("C)")+ # Want to be 3rd plot
  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=20, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "gray50", linewidth = 0.6, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA,linewidth=1.5))


##### Total treatments PER YEAR ####

treat.t$year <- floor(treat.t$tt) # treat.t calculated in treatments_2loci 
treat.y <- aggregate (freq ~ year + treat, data=treat.t, FUN=sum)

treat.y.plot <- ggplot(treat.y, 
                     aes(x=year, y=freq, col=treat, group=treat))+
  geom_line(lwd=0.9)+
  geom_point()+
  #scale_color_manual(name="Treatment",
   #                    labels = c("X","Y"),
    #                   values=c("#FF0000",'#00CCFF', "#FF9933"))+
  scale_y_continuous(name="Total treatments per year")+
  scale_x_continuous(name="Year", breaks=seq(0,15,1))+
  ggtitle("B)")+
  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=18, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        strip.text = element_text(size=15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.key= element_rect(fill = NA),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.line = element_line(colour = "gray50", linewidth = 0.6, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA,linewidth=1.5))


############ Treatments PER WEEK
#treat.plot <- ggplot(treat.t, 
#                     aes(x=tt, y=freq, col=treat))+
#  geom_line(lwd=1)+
#  scale_x_continuous(name="Year", limits=c(0,years), breaks=seq(0,years,by=1))+
#  scale_y_continuous(name="Number of treatments")+
#  #scale_color_manual(name="Treatment",
#   #                  labels = c("X","Y"),
#    #                 values=c("#FF0000","#FF9933"))+
#  ggtitle("D)")+
#  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
#        axis.text.y = element_text(size = 14),
#        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
#        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
 #       plot.title = element_text(size=20, margin = margin(t = 0, r = 0, b = 10, l = 0)),
  #      legend.title = element_text(size = 18), 
#        legend.text = element_text(size = 14),
 #       legend.key= element_rect(fill = NA),
  #      panel.grid.major.x = element_blank(),
   #     panel.grid.major.y = element_blank(),
#        panel.background = element_rect(fill = "white"),
 #       axis.line = element_line(colour = "gray50", size = 0.6, linetype = "solid"),
  #      panel.border = element_rect(colour = "black", fill=NA,size=1.5))



pop.plots <- ggarrange(p.plot, genoNplot, ab.plot, treat.y.plot, nrow=4)

ggsave(file = sprintf("outputs/2loci/plots/pop.plots_%s_%s.png", scenarios, sim), pop.plots, # Save as png
       width = 10, height = 20, dpi = 1000, units = "in", device='png')


#### Save summary files for WHOLE METAPOPULATION ####

write.csv(gf, file = sprintf("outputs/2loci/pop_summary_%s_%s.csv", scenarios, sim), row.names = F) # Save whole-population summary



# Save p 
colnames(gf.p)[3] <- sim
write.csv(gf.p, file = sprintf('outputs/2loci/p/p_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data

geno.pp <- gf[,c("tt","geno","geno.prop")]
colnames(geno.pp)[3]<-sim
write.csv(geno.pp, file = sprintf('outputs/2loci/geno.prop/geno.prop_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data


# Save mean abundance (across farms) in a separate folder - will combine all files in folder to compare abundance across sims
colnames(m.ab)[2] <- sim
write.csv(m.ab, file = sprintf('outputs/2loci/abund/mab_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data

# Mean abundance by Zone
m.abz <- subset(popdyn, geno=="SSUU" & stage=="ad")
m.abz <- aggregate(abund.stage ~ tt + zone, data=m.abz, FUN = mean)
colnames(m.abz)[3] <- sim
write.csv(m.abz, file = sprintf('outputs/2loci/abundz/mabz_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data


# Mean abundance of MOTILES (adults & pre-adults)
m.ma <- subset(popdyn, geno=="SSUU" & (stage=="ad" | stage=="pa"))
m.ma <- aggregate(abund.stage ~ tt + farm, data=m.ma, FUN = sum) # Add ad & pa abundances on each farm
m.ma <- aggregate(abund.stage ~ tt, data=m.ma, FUN = mean) # Average across farms
colnames(m.ma)[2] <- sim
write.csv(m.ma, file = sprintf('outputs/2loci/abund.mot/mma_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data

# Mean motile abundance by Zone
m.maz <- subset(popdyn, geno=="SSUU" & (stage=="ad" | stage=="pa"))
m.maz <- aggregate(abund.stage ~ tt + farm + zone, data=m.maz, FUN = sum) # Add ad & pa abundances on each farm
m.maz <- aggregate(abund.stage ~ tt + zone, data=m.maz, FUN = mean) # Average across farms
colnames(m.maz)[3] <- sim
write.csv(m.maz, file = sprintf('outputs/2loci/abund.motz/mmaz_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data



########## MAP AT 4 POINTS SIMULATION ############

# Subset data at 4 time-points
ad.t1 <- subset(adrr.p, tt == 1)
ad.t2 <- subset(adrr.p, tt == round(years*0.4))
ad.t3 <- subset(adrr.p, tt == round(years*0.7))
ad.t4 <- subset(adrr.p, tt == max(tt))

####### Scale points by ABUNDANCE ############

map.a <- function(data, norwaymap, label, legend, ytext){
ggplot()+
  geom_map(data = norwaymap, aes(map_id = region), # Draw map
           map = norwaymap, colour = NA, fill = "grey85")+
  geom_point(data= data, ## Insert data subset
             aes(x=lon, y=lat, 
                 size=abund.stage, 
                 col=p,
                 shape = H))+ # Shape of point indicates host type at farm
  facet_wrap(. ~ p.locus, # Facet plots by R/S and T/U loci
             nrow=2,
             labeller = as_labeller(c(`pR` = "R frequency", `pT` = "T frequency")))+
  scale_color_gradientn(name="Allele frequency", # Colour gradient based on p
                        limits=c(0,1), 
                        colors = c("#00CCFF","#FF9933","#FF0000"), # Set a gradient along these colours (blue to red)
                        breaks=c(0,0.5,1))+
  scale_size(name="Abundance\n(adult lice/fish)", # Scale point size according to number of lice
             limits=c(0.01, 10), # Only plot farms with at least 1 adult per 100 fish
             trans="log10", # Log scale
             range=c(0.5,2.75))+ # max and min size of points
  scale_shape_manual(name="Strategies",
                     values=c(17, 16, 15, 18, 2, 1, 0, 6))+
  scale_x_continuous(limits=c(4.5,13),
                     name = expression(Longitude~(degree~E)))+
  scale_y_continuous(limits=c(58,66.6),
                     name = expression(Latitude~(degree~N)))+
  coord_fixed(1.8)+ # Fixed aspect ratio (so map not distorted)
  geom_hline(yintercept = c(lwr.latR, upr.latR),
             alpha = ifelse(lwr.latR > 58 | upr.latR < 67, # If it's a simulation where R allele starts within subset of study area, 
                            1, # add lines marking the upper and lower starting range
                            0), # Otherwise, lines are transparent
             lty=2)+
  geom_hline(yintercept = c(lwr.latT, upr.latT),
             alpha = ifelse(lwr.latT > 58 | upr.latT < 67, # If it's a simulation where T allele starts within subset of study area, 
                            1, # add lines marking the upper and lower starting range
                            0), # Otherwise, lines are transparent
             lty=3)+ # Different line than for R allele
  ggtitle(label)+ # Label according to function
  theme(axis.text = element_text(size=10),
        axis.title.y = ytext, # Only label y-axis on leftmost plot (function)
        axis.title.x = element_text(size=10,margin = margin(t = 8, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=20),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key= element_rect(fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.1),
        legend.position = legend) # Only add legend for final panel in series 
}

map.a.t1 <- map.a(ad.t1, norwaymap, "A)", "none", element_text(size=10, margin = margin(t = 0, r = 8, b = 0, l = 0))) # Run plots for each subset
map.a.t2 <- map.a(ad.t2, norwaymap, "B)","none", element_blank()) # label = A), B), etc.
map.a.t3 <- map.a(ad.t3, norwaymap, "C)","none", element_blank()) # Only add legend for last plot
map.a.t4 <- map.a(ad.t4, norwaymap, "D)","right", element_blank()) # Only provide y-axis label for left-most plot

maps.ab <- ggarrange(map.a.t1, map.a.t2, map.a.t3, map.a.t4, nrow=1)

ggsave(file = sprintf("outputs/2loci/plots/mapsAb_%s_%s.png", scenarios, sim), maps.ab, # Save as png
       width=15, height=15, dpi = 500, units = "in", device='png')

#png(file = sprintf("outputs/plots/mapend.abund_%s_%s.png", iter, sim))
#print(map.end.ab)
#dev.off()

