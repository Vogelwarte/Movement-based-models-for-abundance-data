# Instructions

Information for reproducing the results in Carrizo Vergara et al. (2024).

## R functions

The files "Sim_Capture.R", "Sim_Snapshot.R" and "Sim_EcoDiff.R" contain R functions to simulate the proposed models.
The files "Moments_Capture.R", "Moments_Snapshot.R" and "Moments_EcoDiff.R" contain functions to compute mean vectors and covariance matrices for the proposed models.
The files "EcoDiff_ll.R" and "Gauss_ll.R" contain custom functions to compute the exact log-likelihood of an EcoDiff model and of a Gaussian random vector respectively. The last one can be regularized and continuity-corrected.

## Fly data and fitting

The R file "RoverFlyData.R" initializes the data set used. This data set consists only on the information about Rover fruit flies during the first two time takes (0.5 and 1.5 hours after first release) on the experiment carried out by Edelsparre et al. (2021). The data is saved in a matrix called `Q.obs`. The spatial regions are saved in the matrix `Regions.B.L` which contains the bottom-left coordinates of the squares surrounding the traps. The double `dx` indicates the length of each side of the squares. The times of takes are saved as `t`. The initial time and the spatial dimension are saved respectively in `t0` and `d`.




## References

Edelsparre, A. H., Hefley, T. J., Rodríguez, M. A., Fitzpatrick, M. J., & Sokolowski, M. B. (2021). Scaling up: understanding movement from individual differences to population-level dispersal. bioRxiv, 2021-01.


