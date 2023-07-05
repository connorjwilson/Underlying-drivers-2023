README for materials accompanying Wilson et al. (2023) “Unveiling the underlying drivers of Phanerozoic marine diversification”. This file explains the uses of supplied code and data. The material is formatted as an R project for ease of use; after opening 'Underlying-drivers.Rproj', alterations to directory structure should not be needed.



1. Input Data

The file 'occ_data.csv' contains the raw fossil observation data as downloaded from the Paleobiology Database (PBDB), including the URL needed to access the data directly from the PBDB.

'fragmentation_data.csv' contains continental fragmentation data in the form of the fragmentation index from Zaffos et al. (2021).

'stages.csv' contains the temporal data for the time bins used in the analyses.




2. Analyses and results

'normalize_frag_data.R' takes the raw fragmentation index data and normalizes it for the calculation of the fragmentation rate.

'fragmentation-rate.R' takes the normalized fragmentation index ('lfrag_norm.txt') as an input and outputs the fragmentation rate as the file 'smoothed_frag_rate.csv'

'CMR_analysis.R' uses fossil observation data from the PBDB ('occ_data.csv') to estimate time series of sampling, origination, and extinction (fossil time series). The estimates are saved as 'pradel_reulsts.RData' and standard deviations as 'pradel_stdev.RData'. Due to its larger file size, 'occ_data.csv' can be found in the supplementary materials. To run the code, move 'occ_data.csv' into the 'data' folder.

The script 'indiv-process-analysis.R' performs the standalone analyses, testing the possibility of hidden layers underlying the fossil time series of sampling, origination, and extinction. The outputs of these analyses are saved as: 'single_samp.RData', 'single_orig.RData', and 'single_ext.RData' respectively.

'fragindex-div-analysis.R' performs the link analyses comparing the fossil time series with the fragmentation index, and 'fragrate-div-analysis.R' does the link analyses for the fossil time series and fragmentation rate. The outputs of these scripts are in the folder 'results'.

The results of the fragmentation index analyses are saved as 'index_samp.RData', 'index_orig.RData', and 'index_ext.RData'. The results of the fragmentation rate analyses are in 'rate_samp.RData', 'rate_orig.RData', and 'rate_ext.RData'. The folder 'layer_reconstructions' contains the estimated top and hidden layers of the fossil time series: 'orig.mod.RData', 'ext.mod.RData', and 'samp.mod.RData'.

'figures.R' creates the plots seen in the paper.




