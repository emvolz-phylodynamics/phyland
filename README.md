# phyland: Phylogeography using structured coalescent models

This package implements classic ['island' models](https://en.wikipedia.org/wiki/Isolation_by_distance) for migration of lineages between discrete demes. Typical usage would involve using a dated phylogeny and location of sampling for each lineage to infer migration rates between different spatial locations by maximum likelihood. The effective population size in each deme is also estimated. A design matrix can be provided to estimate a subset of rates or to test the hypothesis that rates are equal. 

To cite:
* SDW Frost, Dearlove, B, and Volz, EM, Phylodynamic modeling of the number of host jumps in a zoonosis, in preparation

# Installing the package

This package is an extension of the [phydynR](https://github.com/emvolz-phylodynamics/phydynR) package and requires the `bbmle` package. 
To install the devel version of the package, type:

```{r install, eval=FALSE}
devtools::install_github("emvolz-phylodynamics/phyland")
```
Note that this requires the package *devtools* installed. Or clone the repository and use `R CMD INSTALL <path to repository>`. 

Note if using a Mac, you will need lgfortran and lquadmath libraries installed before hand in order to compile `phydynR`, one of the dependencies.

### Contributors:
- Erik Volz (@emvolz-phylodynamics)
- Simon Frost (@sdwfrost)

Maintainer: Erik Volz (erik.volz@gmail.com)

# Middle East Respiratory Syndrome 

This demonstration is based on a [recent analysis](http://www.biorxiv.org/content/early/2017/08/10/173211) by Gytis et al. of migration of MERS lineages between camels and humans. They did a Bayesian analysis in BEAST using the [structured coalescent](https://github.com/blab/mers-structure) implemented by [Tim Vaughan](https://github.com/tgvaughan). We will check that we get similar results with maximum likelihood. 

We load the maximum clade credibility tree: 
```
require(phyland)
require(ape)
mcc <- read.nexus(system.file( 'inst/MERS_274_sCoal.combinedTyped.mcc.tree', package='phyland') )
```

Next we fit the model by maximum likelihood. Note that `delimiter` and `index` options tell us where to look in phylogeny tip labels for the deme of sampling, which in this case will be the host species (camel or human). Alternatively, we could provide a regular expression to infer the deme from tip labels (`regex` option). 
```
fit <- phylandml(mcc, delimiter="|", index=3)
fit

Summary of log transformed parameters:
              Estimate Std. Error    z value        Pr(z)
human        -1.092695 0.13102339 -8.3396972 7.448103e-17
camel         1.292827 0.08786038 14.7145640 5.198017e-49
human->camel -7.491869 8.07100469 -0.9282449 3.532806e-01
camel->human  1.699715 0.16106809 10.5527744 4.931886e-26

Design matrix
      human camel
human    NA     2
camel     4    NA

Estimated effective sizes and rates:
                         
human        0.3353114669
camel        3.6430716096
camel->human 5.4723887166
human->camel 0.0005575999
```
The 'human' and 'camel' parameters provide the effective population size in both demes. 

Finally, we can derive confidence intervals using likelihood profiles. The `ncpu` option tells the profiler to spread the job among multiple cpus. 
```
confint( fit, whichparms= c('camel->human', 'human->camel', 'human' ) , ncpu = 3)

       camel..human
2.5 %      3.912770
97.5 %     7.418204
       human..camel     human
2.5 %            NA 0.2598891
97.5 %   0.07553906 0.4355089
```
Note the `NA` value indicates the human to camel rate was indistinguishable from zero. 
