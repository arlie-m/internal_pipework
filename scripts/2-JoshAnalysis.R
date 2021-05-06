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
community_matrix_AB <- select(community_matrix_sites, c_f_conchoderma_na:other_bivalve)
community_site_info <- select(community_matrix_sites, id, area, stage, quadrat_number)

#Make a presence-absence version of community data
community_matrix_PA <- decostand(community_matrix_AB, method = "pa")

#Perform PERMANOVA - right now PERMANOVAs use bray-curtis distances, as I think from NMDS 
#output it also uses the bray-curtis distance, so these results are directly comparable to 
#the ordination visualisation

#PERMANOVA for abundance data
ABpermanova <- adonis2(community_matrix_AB ~ stage, community_site_info, permutations=1000,
                       method="bray")
ABpermanova
##         Df SumOfSqs      R2      F Pr(>F)
##stage     1   0.4517 0.03426 1.3834 0.1628
##Residual 39  12.7328 0.96574              
##Total    40  13.1844 1.00000 
#Cannot reject null hypothesis that centroids and/or dispersion of the groups are the same 

#PERMANOVA for presence-absence data
PApermanova <- adonis2(community_matrix_PA ~ stage, community_site_info, permutations=1000,
                       method="bray")
PApermanova
##          Df SumOfSqs      R2      F   Pr(>F)    
##stage     1   1.3204 0.14654 6.6961 0.000999 ***
##Residual 39   7.6903 0.85346                    
##Total    40   9.0107 1.00000 
#CAN reject null hypothesis.

#So basically: heaps of variation in abundances WITHIN each stage (cannot detect overall
#differences between stages), BUT just looking at the presence of species, the four stages ARE
#significantly different i.e. different communities of organisms are in the four stages. 
#Suggests multiple things going on - there are broad differences between the four stages in 
#terms of what species are there, but there is also an area-specific effect leading to very 
#variable abundances within each stage.

#Analysis of nestedness
#For now, just use least constrained null model (we need to talk about this) - 
#the nestedness algorithms we use need serious consideration 
#Will also increase simulation number 

#Look at nestedness of presence-absence data
oecosimu(community_matrix_PA, nestednodf, "r00", nsimul=500)
##nullmodel method ‘r00’ with 500 simulations
##
##alternative hypothesis: statistic is less or greater than simulated values
##
##N columns  : 43.11143 
##N rows     : 60.68859 
##NODF       : 51.90001 
##Matrix fill: 0.1659726 
##
##          statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)   
##N.columns    43.111 28.761 17.939 16.201 17.918 19.781 0.001996 **
##N.rows       60.689 49.952 17.921 16.377 17.931 19.668 0.001996 **
##NODF         51.900 41.510 17.930 16.450 17.919 19.664 0.001996 **

#i.e. species presences in the community matrix are signifcantly more nested than expected 
#by chance - i.e. species-poor samples tend to be nested subsets of species-rich ones.
#This is good - especially if we can show a linear decrease in species richness through the
#ship's system!

#Look at nestedness of abundnace data
#Algorithms require integer values
community_matrix_AB_rounded <- ceiling(community_matrix_AB) #round up because otherwise some
#columns would be entirely 0
oecosimu(community_matrix_AB_rounded, nestednodf, "r00_ind", nsimul=500)
##nullmodel method ‘r00_ind’ with 500 simulations
##
##alternative hypothesis: statistic is less or greater than simulated values
##
##N columns  : 43.11143 
##N rows     : 60.68859 
##NODF       : 51.90001 
##Matrix fill: 0.1659726 
##
##          statistic      SES   mean   2.5%    50%  97.5% Pr(sim.)  
##N.columns    43.111 -2.69421 57.028 45.279 57.465 65.748  0.02196 *
##N.rows       60.689  0.72595 56.984 44.892 57.497 65.268  0.51697  
##NODF         51.900 -1.01389 57.006 45.303 57.534 65.118  0.30539 

#In terms of species abundances, they are NOT more nested than expected by chance!
#This agrees thematically with the PERMANOVA results 










