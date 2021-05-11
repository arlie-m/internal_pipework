#Josh plays with some stuff
#Apologies for not doing in RMarkdown
#I am a noob 

#PERMANOVA - to get statistical significance for NMDS results
#(are the clusterings significant?)

#Need to make sure have community matrix and site data with matching numbers of rows
community_matrix_sites <- community_matrix %>% 
  rownames_to_column(var = "id") %>% 
  left_join(site_info)
#Re-extract community matrix and site info data now they are aligned
community_matrix_AB <- select(community_matrix_sites, c_f_corophium_na:conchoderma_auritum)
community_site_info <- select(community_matrix_sites, id, area, stage, quadrat_number)

#Make a presence-absence version of community data
community_matrix_PA <- decostand(community_matrix_AB, method = "pa")

###########################################################################################

#Perform PERMANOVA - right now PERMANOVAs use bray-curtis distances, as I think from NMDS 
#output it also uses the bray-curtis distance, so these results are directly comparable to 
#the ordination visualisation

#PERMANOVA for abundance data
ABpermanova <- adonis2(community_matrix_AB ~ stage, community_site_info, permutations=1000,
                       method="bray")
ABpermanova
##         Df SumOfSqs      R2      F Pr(>F)
##stage     1   0.4827 0.03656 1.4798 0.1249
##Residual 39  12.7221 0.96344              
##Total    40  13.2048 1.00000
#Cannot reject null hypothesis that centroids and/or dispersion of the groups are the same 

#PERMANOVA for presence-absence data
PApermanova <- adonis2(community_matrix_PA ~ stage, community_site_info, permutations=1000,
                       method="bray")
PApermanova
##         Df SumOfSqs      R2      F   Pr(>F)    
##stage     1   1.3602 0.15336 7.0643 0.000999 ***
##Residual 39   7.5093 0.84664                    
##Total    40   8.8695 1.00000 
#CAN reject null hypothesis.

#So basically: heaps of variation in abundances WITHIN each stage (cannot detect overall
#differences between stages), BUT just looking at the presence of species, the four stages ARE
#significantly different i.e. different communities of organisms are in the four stages. 
#Suggests multiple things going on - there are broad differences between the four stages in 
#terms of what species are there, but there is also an area-specific effect leading to very 
#variable abundances within each stage.

###########################################################################################

#Analysis of nestedness

#First, look at nestedness of presence-absence data:

#Use least-constrained null model (only total count maintained)

oecosimu(community_matrix_PA, nestednodf, "r00", nsimul=1000)
##N columns  : 42.16523 
##N rows     : 60.64611 
##NODF       : 51.87331 
##Matrix fill: 0.1651032 
##
##          statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    42.165 28.000 17.811 16.148 17.794 19.608 0.000999 ***
##N.rows       60.646 48.669 17.836 16.199 17.790 19.582 0.000999 ***
##NODF         51.873 41.577 17.824 16.283 17.798 19.442 0.000999 ***

#Use null model where individual site richness is constrained

oecosimu(community_matrix_PA, nestednodf, "r0", nsimul=1000)
##N columns  : 42.16523 
##N rows     : 60.64611 
##NODF       : 51.87331 
##Matrix fill: 0.1651032 
##
##          statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    42.165 17.995 23.324 21.299 23.321 25.460 0.000999 ***
##N.rows       60.646 45.033 21.584 20.011 21.537 23.582 0.000999 ***
##NODF         51.873 34.233 22.410 20.750 22.379 24.143 0.000999 ***

#Use null model where species prevalences are constrained

oecosimu(community_matrix_PA, nestednodf, "c0", nsimul=1000)
##N columns  : 42.16523 
##N rows     : 60.64611 
##NODF       : 51.87331 
##Matrix fill: 0.1651032 
##
##          statistic     SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    42.165 15.4292 26.032 24.124 26.010 28.125 0.000999 ***
##N.rows       60.646  8.2024 43.598 39.589 43.577 47.738 0.000999 ***
##NODF         51.873 11.3218 35.259 32.480 35.240 38.271 0.000999 ***

#Use null model where both site richness and species prevalences are constrained 
oecosimu(community_matrix_PA, nestednodf, "quasiswap", nsimul=1000)
##N columns  : 42.16523 
##N rows     : 60.64611 
##NODF       : 51.87331 
##Matrix fill: 0.1651032 
##
##          statistic      SES   mean   2.5%    50%  97.5% Pr(sim.)  
##N.columns    42.165  1.64850 40.140 37.600 40.238 42.292  0.07093 .
##N.rows       60.646 -1.48999 63.236 59.253 63.420 66.074  0.15485  
##NODF         51.873 -0.33869 52.273 49.825 52.353 54.333  0.68831


#So basically: we observe that species-poor sites tend to be nested subsets of species-rich
#ones (r00 model). Allowing for the fact that some sites will be more favourable than others (r0),
#we still see this pattern. Allowing for the fact that some species are inherently more common than
#others (c0), we still see this pattern. However, allowing for both these things (quasiswap) we do 
#NOT see this pattern: nestedness is a combined product of the suitability of certain sites and the
#differential prevalence of species.


#Next, look at nestedness of abundance data

#Algorithms require integer values
community_matrix_AB_rounded <- ceiling(community_matrix_AB) #round up because otherwise some
#columns would be entirely 0

oecosimu(community_matrix_AB_rounded, nestednodf, "r00_ind", nsimul=1000)
#Don't think this code is right - need to have the weighted version of nestednodf 
#but can't get the nested version to work with oecosimu right now









