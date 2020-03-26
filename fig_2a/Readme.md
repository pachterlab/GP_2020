# fig_2a
The directory `fig2a` contains codes to reproduce Figure 2a, which demonstrates the effect of approximation order on a sample probability mass function. The following files are provided: 

`fig2a.m`: MATLAB code to compute distributions. The computation uses a simulation (gray histogram)  and special-function approximations over a user-defined grid (red line). 
`fig_2a.pdf`: Published MATLAB output of the visualization code.

`gg_200110_gillespie_geom_1.m`: Gillespie algorithm implementation used to generate realizations of the stochastic system.
`gg_200325_analyt_geom_tdep_vec_31.m`: Implementation of the routine for calculating approximations. 
`mhygfx.m`: Implementation of the Gaussian hypergeometric function by Ben Barrowes ([https://www.mathworks.com/matlabcentral/fileexchange/6218-computation-of-special-functions](https://www.mathworks.com/matlabcentral/fileexchange/6218-computation-of-special-functions)), based on FORTRAN code by Zhang and Jin.
`me1z_gg200131_comb_4.m`: Implementation of the exponential integral routine described in the **Supplementary Note**.
All others: Padé, Chebyshev, and Taylor approximations to the exponential integral. 
