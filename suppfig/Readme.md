# suppfig
The directory `suppfig` contains codes to reproduce Supplementary Figure 1, which demonstrates the combination of approximations. In the algorithm, these are used to compute the exponential integral for the Laurent expansion of the degenerate system solution. The following files are provided: 

`supp_fig_1.m`: MATLAB code to compute and visualize each approximation schema. The function compares the results with the native MATLAB function `expint`, outputting both log relative error and runtime improvement. The visualization conventions are described in comments and the **Supplementary Note**.

`supp_fig_1.pdf`: Published MATLAB output of the visualization code.

`me1z_gg200131_comb_4.m`: Implementation of the combined exponential integral routine described in the **Supplementary Note**.

All others: Padé, Chebyshev, and Taylor approximations to the exponential integral. 
