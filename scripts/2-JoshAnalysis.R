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
community_matrix_AB <- select(community_matrix_sites, -id)
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
#         Df SumOfSqs      R2      F  Pr(>F)  
#stage     1   0.8247 0.08526 3.6351 0.01099 *
#Residual 39   8.8477 0.91474                 
#Total    40   9.6724 1.00000                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#CAN reject null hypothesis? 

#So, unlike previous time the analysis was run (on data without verified species), 
#we can detect overall differences between stages on the abundance matrix.


#PERMANOVA for presence-absence data
PApermanova <- adonis2(community_matrix_PA ~ stage, community_site_info, permutations=1000,
                       method="bray")
PApermanova
##         Df SumOfSqs      R2      F   Pr(>F)    
##stage     1   0.4618 0.15706 7.2668 0.000999 ***
##Residual 39   2.4784 0.84294                    
##Total    40   2.9402 1.00000  
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
#N columns  : 51.34796 
#N rows     : 80.44558 
#NODF       : 65.89677 
#Matrix fill: 0.2355741 
#
#statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
#N.columns    51.348 28.868 24.819 23.065 24.816 26.638 0.000999 ***
#N.rows       80.446 60.828 24.828 23.000 24.850 26.632 0.000999 ***
#NODF         65.897 47.811 24.823 23.191 24.812 26.545 0.000999 ***

#Use null model where individual site richness is constrained

oecosimu(community_matrix_PA, nestednodf, "r0", nsimul=1000)
#N columns  : 51.34796 
#N rows     : 80.44558 
#NODF       : 65.89677 
#Matrix fill: 0.2355741 
#
#statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
#N.columns    51.348 26.839 27.561 25.865 27.559 29.280 0.000999 ***
#N.rows       80.446 95.989 27.882 26.930 27.839 29.011 0.000999 ***
#NODF         65.897 55.755 27.721 26.489 27.713 29.066 0.000999 ***

#Use null model where species prevalences are constrained

oecosimu(community_matrix_PA, nestednodf, "c0", nsimul=1000)
#N columns  : 51.34796 
#N rows     : 80.44558 
#NODF       : 65.89677 
#Matrix fill: 0.2355741 
#
#statistic     SES   mean   2.5%    50%  97.5% Pr(sim.)    
#N.columns    51.348 14.9415 37.895 36.257 37.858 39.842 0.000999 ***
#N.rows       80.446  7.6746 62.071 57.127 62.233 66.400 0.000999 ***
#NODF         65.897 10.4494 49.983 47.015 50.038 52.768 0.000999 ***


#Use null model where both site richness and species prevalences are constrained 
oecosimu(community_matrix_PA, nestednodf, "quasiswap", nsimul=1000)
#N columns  : 51.34796 
#N rows     : 80.44558 
#NODF       : 65.89677 
#Matrix fill: 0.2355741 

#statistic      SES   mean   2.5%    50%  97.5% Pr(sim.)  
#N.columns    51.348  1.34289 50.143 48.245 50.204 51.757  0.16484  
#N.rows       80.446 -1.86316 80.961 80.342 80.990 81.404  0.08891 .
#NODF         65.897  0.69391 65.552 64.467 65.581 66.402  0.51049  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


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

#Stuff to do double dendrogram heat map:
##Heatmap:
#image(rectangular_matrix,col=your_colourscale,xaxt="n",yaxt="n") ##rectangular_matrix is a numeric matrix with the colour levels you want to represent.
#
##Horizontal colour values (levels of invasion):
#image(as.matrix(numeric_vector),col=brewer.pal(5, "OrRd"),xaxt="n",yaxt="n") ##the colour scale was set to 5 levels for the "OrRd" set of colours (see help file for brewer.pal). The numeric_vector has the same number of elements as there are columns in the rectangular_matrix above.
#
##Vertical colour values (driver categories):
#image(t(as.matrix(another_numeric_vector)),col=brewer.pal(length(levels(typesVariables)),"Set1"),xaxt="n",yaxt="n") ##the colour scale was set to the number of levels corresponding to the number of categories for the "Set1" set of colours (see help file for brewer.pal). The numeric_vector has the same number of elements as there are rows in the rectangular_matrix above. Note that the matrix was transposed to get it vertically, but you can do this by hand in illustrator too.
#
##Dendrogram:
#dendrogram <- hclust(distance_matrix)
#plot(dendrogram ,hang=-1)
#
#For the rotated dendrogram, you need to rotate it in Illustrator.









