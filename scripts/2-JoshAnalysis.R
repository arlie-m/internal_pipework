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
##N columns  : 52.12222 
##N rows     : 80.47978 
##NODF       : 66.65547 
##Matrix fill: 0.2414634 
##
##          statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    52.122 28.871 25.428 23.626 25.432 27.272 0.000999 ***
##N.rows       80.480 60.098 25.463 23.684 25.494 27.258 0.000999 ***
##NODF         66.655 48.315 25.446 23.857 25.429 27.189 0.000999 ***

#Use null model where individual site richness is constrained

oecosimu(community_matrix_PA, nestednodf, "r0", nsimul=1000)
##N columns  : 52.12222 
##N rows     : 80.47978 
##NODF       : 66.65547 
##Matrix fill: 0.2414634 
##
##          statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    52.122 25.050 28.210 26.278 28.248 30.018 0.000999 ***
##N.rows       80.480 92.183 28.578 27.534 28.533 29.784 0.000999 ***
##NODF         66.655 53.135 28.399 27.013 28.398 29.922 0.000999 ***

#Use null model where species prevalences are constrained

oecosimu(community_matrix_PA, nestednodf, "c0", nsimul=1000)
##N columns  : 52.12222 
##N rows     : 80.47978 
##NODF       : 66.65547 
##Matrix fill: 0.2414634 
##
##          statistic     SES   mean   2.5%    50%  97.5% Pr(sim.)    
##N.columns    52.122 15.0916 38.680 37.082 38.617 40.583 0.000999 ***
##N.rows       80.480  7.7312 61.892 56.961 61.982 66.585 0.000999 ***
##NODF         66.655 10.4793 50.576 47.514 50.540 53.676 0.000999 ***


#Use null model where both site richness and species prevalences are constrained 

oecosimu(community_matrix_PA, nestednodf, "quasiswap", nsimul=1000)
##N columns  : 52.12222 
##N rows     : 80.47978 
##NODF       : 66.65547 
##Matrix fill: 0.2414634 
##
##statistic     SES   mean   2.5%    50%  97.5% Pr(sim.)  
##N.columns    52.122  0.8833 51.259 49.102 51.354 52.917  0.39261  
##N.rows       80.480 -1.8819 81.010 80.411 81.049 81.457  0.08492 .
##NODF         66.655  0.2791 66.507 65.386 66.549 67.456  0.84016


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









