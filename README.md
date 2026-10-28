The MATLAB code bayes3parICrot_S22LOD_public.m performs the Bayesian inversion using the gravity time-series (file cs22_ts_errors.nc) 
and the length-of-day time-series corrected from atmospheric, oceanic, hydrological and sea-level angular momentum (file Chiz1976-2026_AOHSlcorr.txt). 
The results are the joint probability density functions between:
* the gravitational coupling constant \Gamma (its log10) and the degree-2 topographic amplitude at the inner-core boundary (\delta h),
* the gravitational coupling constant \Gamma (its log10) and the characteristic time of viscous deformation (\tau),
* the degree-2 topographic amplitude at the inner-core boundary (\delta h) and the characteristic time of viscous deformation (\tau).
The marginal laws as well as the maximum a posteriori probability values (MAP), the mean values and their standard deviations are also computed for these parameters. 
