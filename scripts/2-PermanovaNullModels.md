2-PermanovaNullModels
================
Arlie McCarthy
25/06/2021

# Statistical analysis using PERMANOVA and Null Models

\#PERMANOVA - to get statistical significance for NMDS results

(are the clusterings significant?)

``` r
#Need to make sure have community matrix and site data with matching numbers of rows
community_matrix_sites <- community_matrix %>% 
  left_join(site_info)
```

    ## Joining with `by = join_by(id)`

``` r
#Re-extract community matrix and site info data now they are aligned
community_matrix_AB <- select(community_matrix_sites, -id)
community_site_info <- select(community_matrix_sites, id, area, stage, quadrat_number)
#Make a presence-absence version of community data
community_matrix_PA <- decostand(community_matrix_AB, method = "pa")
```

# Perform PERMANOVA

Right now PERMANOVAs use bray-curtis distances, as the NMDS also uses
the bray-curtis distance, these results are directly comparable to the
ordination visualisation

## PERMANOVA for abundance data

``` r
ABpermanova <- adonis2(community_matrix_AB ~ stage, community_site_info, permutations=1000,
                       method="bray")
ABpermanova
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## adonis2(formula = community_matrix_AB ~ stage, data = community_site_info, permutations = 1000, method = "bray")
    ##          Df SumOfSqs      R2      F   Pr(>F)   
    ## Model     1   0.9094 0.09027 3.8697 0.004995 **
    ## Residual 39   9.1653 0.90973                   
    ## Total    40  10.0747 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

CAN reject null hypothesis

So, we can detect overall differences between stages on the abundance
matrix.

## PERMANOVA for presence-absence data

``` r
PApermanova <- adonis2(community_matrix_PA ~ stage, community_site_info, permutations=1000,
                       method="bray")
PApermanova
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## adonis2(formula = community_matrix_PA ~ stage, data = community_site_info, permutations = 1000, method = "bray")
    ##          Df SumOfSqs      R2      F   Pr(>F)    
    ## Model     1  0.44868 0.16069 7.4668 0.000999 ***
    ## Residual 39  2.34350 0.83931                    
    ## Total    40  2.79218 1.00000                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

CAN reject null hypothesis.

So basically: some variation in abundances WITHIN each stage (less
significant differences between stages), BUT just looking at the
presence of species, the four stages are still significantly different
i.e. different communities of organisms are in the four stages. Suggests
multiple things going on - there are broad differences between the four
stages in terms of what species are there, but there is also an
area-specific effect leading to very variable abundances within each
stage.

# Analysis of nestedness

Look at nestedness of presence-absence data:

``` r
#Use least-constrained null model (only total count maintained)
oecosimu(community_matrix_PA, nestednodf, "r00", nsimul=1000)
```

    ## oecosimu object
    ## 
    ## Call: oecosimu(comm = community_matrix_PA, nestfun = nestednodf, method
    ## = "r00", nsimul = 1000)
    ## 
    ## nullmodel method 'r00' with 1000 simulations
    ## 
    ## alternative hypothesis: statistic is less or greater than simulated values
    ## 
    ## N columns  : 53.8666 
    ## N rows     : 80.99604 
    ## NODF       : 68.83709 
    ## Matrix fill: 0.2570864 
    ## 
    ##           statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
    ## N.columns    53.867 25.683 26.950 24.827 26.979 28.886 0.000999 ***
    ## N.rows       80.996 54.820 27.021 25.100 27.038 28.881 0.000999 ***
    ## NODF         68.837 44.453 26.989 25.182 27.000 28.813 0.000999 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Use null model where individual site richness is constrained
oecosimu(community_matrix_PA, nestednodf, "r0", nsimul=1000)
```

    ## oecosimu object
    ## 
    ## Call: oecosimu(comm = community_matrix_PA, nestfun = nestednodf, method
    ## = "r0", nsimul = 1000)
    ## 
    ## nullmodel method 'r0' with 1000 simulations
    ## 
    ## alternative hypothesis: statistic is less or greater than simulated values
    ## 
    ## N columns  : 53.8666 
    ## N rows     : 80.99604 
    ## NODF       : 68.83709 
    ## Matrix fill: 0.2570864 
    ## 
    ##           statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
    ## N.columns    53.867 25.191 29.380 27.491 29.402 31.185 0.000999 ***
    ## N.rows       80.996 92.260 30.037 29.081 29.996 31.224 0.000999 ***
    ## NODF         68.837 55.663 29.742 28.442 29.731 31.132 0.000999 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Use null model where species prevalences are constrained
oecosimu(community_matrix_PA, nestednodf, "c0", nsimul=1000)
```

    ## oecosimu object
    ## 
    ## Call: oecosimu(comm = community_matrix_PA, nestfun = nestednodf, method
    ## = "c0", nsimul = 1000)
    ## 
    ## nullmodel method 'c0' with 1000 simulations
    ## 
    ## alternative hypothesis: statistic is less or greater than simulated values
    ## 
    ## N columns  : 53.8666 
    ## N rows     : 80.99604 
    ## NODF       : 68.83709 
    ## Matrix fill: 0.2570864 
    ## 
    ##           statistic     SES   mean   2.5%    50%  97.5% Pr(sim.)    
    ## N.columns    53.867 12.8634 41.055 39.313 41.036 43.118 0.000999 ***
    ## N.rows       80.996  7.3801 62.666 57.515 62.868 67.123 0.000999 ***
    ## NODF         68.837  9.5104 52.980 49.593 53.042 56.104 0.000999 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Use null model where both site richness and species prevalences are constrained 
oecosimu(community_matrix_PA, nestednodf, "quasiswap", nsimul=1000)
```

    ## oecosimu object
    ## 
    ## Call: oecosimu(comm = community_matrix_PA, nestfun = nestednodf, method
    ## = "quasiswap", nsimul = 1000)
    ## 
    ## nullmodel method 'quasiswap' with 1000 simulations
    ## 
    ## alternative hypothesis: statistic is less or greater than simulated values
    ## 
    ## N columns  : 53.8666 
    ## N rows     : 80.99604 
    ## NODF       : 68.83709 
    ## Matrix fill: 0.2570864 
    ## 
    ##           statistic      SES   mean   2.5%    50%  97.5% Pr(sim.)
    ## N.columns    53.867  1.05936 52.746 50.444 52.859 54.459   0.2647
    ## N.rows       80.996 -1.52577 81.415 80.820 81.452 81.853   0.1828
    ## NODF         68.837  0.50129 68.566 67.411 68.618 69.425   0.6563

So basically: we observe that species-poor sites tend to be nested
subsets of species-rich ones (r00 model). Allowing for the fact that
some sites will be more favourable than others (r0), we still see this
pattern. Allowing for the fact that some species are inherently more
common than others (c0), we still see this pattern. However, allowing
for both these things (quasiswap) we do NOT see this pattern: nestedness
is a combined product of the suitability of certain sites and the
differential prevalence of species.
