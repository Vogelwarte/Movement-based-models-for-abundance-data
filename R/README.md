# Instructions

Information for reproducing the results in Carrizo Vergara et al. (2024).

## R functions

The files "Sim_Capture.R", "Sim_Snapshot.R" and "Sim_EcoDiff.R" contain R functions to simulate the proposed models.
The files "Moments_Capture.R", "Moments_Snapshot.R" and "Moments_EcoDiff.R" contain functions to compute mean vectors and covariance matrices for the proposed models.
The files "EcoDiff_ll.R" and "Gauss_ll.R" contain custom functions to compute the exact log-likelihood of an EcoDiff model and of a Gaussian random vector respectively. The last one can be regularized and continuity-corrected.
All the functions have detailed comments for their use, as well as examples.

## Fly data and fitting

The R file "RoverFlyData.R" initializes the data set used. This data set consists only on the information about Rover fruit flies during the first two time takes (0.5 and 1.5 hours after first release) on the experiment carried out by Edelsparre et al. (2021). The data is saved in a matrix called `Q.obs`. The spatial regions are saved in the matrix `Regions.B.L` which contains the bottom-left coordinates of the squares surrounding the traps. The double `dx` indicates the length of each side of the squares. The times of takes are saved as `t`. The initial time and the spatial dimension are saved respectively in `t0` and `d`.


For fitting the proposed models to the data, one must launch the routines in the files "FitCaptureDataMGLE.R", "FitSnapshotDataMGLE.R", "FitEcoDiffDataMGLE.R", or "FitEcoDiffDataMLE.R", depending on the desired model and fitting method. For the Capture and Snapshot models with MGLE fitting method, the computing time on a 1.50 GHz and 16 GB RAM laptot was about 30 min. For the EcoDiff model with MGLE, the computing time was about 15s, while with MLE fitting method was less than 1s. The computation of the continuity corrected version of the Gaussian likelihood, done with the `TruncatedNormal` package, took around 10 min.

## Simulation studies

To perform simulation studies, one must use the routines "SimuStudyCaptureMGLE.R", "SimuStudySnapshotMGLE.R", "SimuStudyEcoDiffMGLE.R", or "SimuStudyEcoDiffMLE.R" depending on the desired model and fitting method. For each routine launch, one must initialize a parameter's setting theoretical value: $(\sigma , v_{x} , v_{y} )$ for the movement model (assumed Brownian with advection), $p$ for Snapshot or EcoDiff, and $\alpha$ for Capture. One must also provide the desired number of individuals $N$.


**IMPORTANT:** Since the simulation-fitting methods take some time, these routines are parallelaized.


## References

Carrizo Vergara, R., Kéry, M., & Hefley, T. (2024). Movement-based models for abundance data. arXiv preprint arXiv:2407.13384.


Edelsparre, A. H., Hefley, T. J., Rodríguez, M. A., Fitzpatrick, M. J., & Sokolowski, M. B. (2021). Scaling up: understanding movement from individual differences to population-level dispersal. bioRxiv, 2021-01.


