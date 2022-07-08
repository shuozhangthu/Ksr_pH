Model and data analysis scripts for "A Model for pH Dependent Strontium Partitioning during Calcite Precipitation from Aqueous Solutions".
 
This directory includes:
-Model input data from this study.
   (1) experimental data of Tang et al. (2008), Lorens (1981) and Tesoriero and Pankow (1996).   Files: lorens.mat, tang.mat, tang_all.mat, tp.mat
   
-Estimation of the optimal parameters of fitting carves of experimental data.
   (1) Script for eatimating the optimal parameters of data of Tang et al. (2008), Lorens (1981) and Tesoriero and Pankow (1996) at 25C. File: optimize_1D_jia.m
   (2) Script for eatimating the optimal parameters of data of Tang et al. (2008) at 5C, 25C, and 40C. File: optimize_5.m, optimize_25.m, optimize_40.m              

-fitting of the experimental data.
   (1) Script for fitting the data of Tang et al. (2008), Lorens (1981) and Tesoriero and Pankow (1996) at 25C. Files: figure_1D_jia.m
   (2) Script for fitting the data of Tang et al. (2008) at 5C, 25C, and 40C. Files: run_5.m, run_25.m, run_40.m, figure_3T_Tang.m
 
-pH_T_Ksr_model.
   (1) A function for calculating calcite precipitation rate and partitioning coefficient based on known seawater chemistry. The inputs are pH, temperature, saturation and growth mechanism. The outputs are Ksr and Rp. Files: function_ion_by_ion.m
   (2) A simple case . Files: case1.m
 

Note: File [solveP.m] is used to calculate the probability of each kink site in run files and function files.

Reference:
Lorens, R.B. (1981) Sr, Cd, Mn and Co distribution coefficients in calcite as a function of calcite precipitation rate. Geochim. Cosmochim. Acta 45, 553–561.
Tesoriero, A.J. and Pankow, J.F. (1996) Solid solution partitioning of Sr2+, Ba2+, and Cd2+ to calcite. Geochimica et Cosmochimica Acta, 60, 1053-1063.
Tang, J., Köhler, S.J. and Dietzel, M. (2008) Sr2+/Ca2+ and 44Ca/40Ca fractionation during inorganic calcite formation: I. Sr incorporation. Geochimica et Cosmochimica Acta 72, 3718-3732.