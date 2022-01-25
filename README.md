# Code from "Asychrony, density dependence, and persistence in an amphibian" (Manuscript ID ECY21-0446)
## Published in Ecology 2022

## Authors
- Freya E. Rowland<sup>1,6</sup>*
- Elizabeth Schyling<sup>1</sup>
- L. Kealoha Freidenburg<sup>1</sup>
- Mark C. Urban<sup>2</sup>
- Jonathan L Richardson<sup>3</sup>
- A.Z. Andis Arietta<sup>1</sup>
- Susan B. Rodrigues<sup>4</sup>
- Adriana D. Rubinstein<sup>1</sup>
- Michael F. Benard<sup>5</sup>
- David K. Skelly<sup>1,4</sup>

<sup>1</sup>School of the Environment, Yale University, 370 Prospect St., New Haven, CT 06511 USA
<sup>2</sup>Department of Ecology and Evolutionary Biology and Center of Biological Risk, University of Connecticut, Storrs, CT 06269 USA
<sup>3</sup>Department of Biology, University of Richmond, VA 23173 USA
<sup>4</sup>Yale Peabody Museum of Natural History, Yale University, New Haven, CT 06511 USA
<sup>5</sup>Department of Biology, Case Western Reserve University, 10900 Euclid Ave., Cleveland, OH 44106 USA
<sup>6</sup>Corresponding author: frowland@usgs.gov

*Current address Columbia Environmental Research Center, US Geological Survey, Columbia MO 65201 USA

## Code

There are five scripts associated with this manuscript all linked to [WoodFrogCode.Rproj](code/WoodFrogCode.Rproj):

1) [NeighborhoodCompCode.R](<code/NeighborhoodCompCode.R>) shows how we computed the weighted effect of competition from neighboring ponds.
2) [RegressionModelsFigures.R](<code/RegressionModelsFigures.R>) shows the Bayesian hierarchical code using rstanarm to allow slope and intercept to vary by pond.
3) [SynchronyModelsFigures.R](<code/SynchronyModelsFigures.R>) shows the calculations for spatial synchrony within years and across years.
4) [YMFmap.R](<code/YMFmap.R>) is the code for making the map figure.
5) [BayesianVariableSelection.R](<code/BayesianVariableSelection.R>) is the script for Bayesian variable selection of variables and models.

All scripts have been archived via Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5903388.svg)](https://doi.org/10.5281/zenodo.5903388) permanent DOI [10.5281/zenodo.5893576](https://doi.org/10.5281/zenodo.5893576).

## Data

The data are published on Dryad [doi:10.5061/dryad.0cfxpnw3r](https://doi.org/10.5061/dryad.0cfxpnw3r). Note: no pond coordinates are included in the released data to keep populations safe from intrusion. 
