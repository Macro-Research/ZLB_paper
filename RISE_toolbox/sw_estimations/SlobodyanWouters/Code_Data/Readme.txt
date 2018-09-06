This file contains description of the data and code supplied together with the paper 

"Learning in a medium-scale DSGE model with expectations based on small forecasting models"
by SERGEY SLOBODYAN AND RAF WOUTERS

accepted for publication at the American Economic Journal: Macroeconomics

The archive contains 3 directories. 

Directory DATA contains file SW_data_Apr09_AEJM.xls which performs derivation of the data used in the paper from sources described in the Appendix to the paper. File usmodel_data_SW_data_april2009.mat is a MATLAB data file containing this calculated and transformed data; it is this data which is used as an input to all estimations described in the paper.

Directory FIG_3 contains Fig_3_plot.m and Fig_3_data.mat, a MATLAB file and the data which are needed to produce Fig 3 in the paper. See the commented section in the Fig_3_plot.m for derivation of the data contained in Fig_3_data.mat from our adaptive learning estimation output and publicly available SPF and real-time inflation data. References to the data are given in the header of Fig_3_plot.m.

Directory CODE contains several directions with the adaptive learning code needed to replicate all the estimates listed in Table 1, plus five-model estimates described in Section V of the paper. Make a particular directory home directory in MATLAB and run the *.mod file contained in the directory. The code runs under DYNARE 3.64 for MATLAB. It WILL NOT run under DYNARE 4. For simplicity, no MCMC is performed; modify the corresponding .mod file if you need to obtain MCMC output for a particular estimation. See Code_Description.pdf for a detailed explanation of what the code does.