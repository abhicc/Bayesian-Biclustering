# Bayesian-Biclustering
Supplementary material for the paper 'Some Bayesian biclustering Methods: Modeling and inference'

# Description of the files

- **BayesianBiclustering_NonMissing_v1.R** - first version of the MCMC sampler for non-missing data
- **BayesianBiclustering_NonMissing_v2.R** - second version of the MCMC sampler for non-missing data
- **BayesianBiclustering_WithMissing.R** - MCMC sampler for datasets with missing entries
- **cppFunc.cpp** - C++ implementation to compute average Rand index (for rows, columns, and elements)
- **SimStudies.R** - additional code for simulation studies 6.1.3, 6.2, and 6.3
- **SimStudies_X.rds** - dataset for simulation studies 6.1.3, 6.2, and 6.3
- **MultBiclust.R** - additional code for simulation study 6.4
- **MultBiclust_X.rds** - dataset for simulation study 6.4
- **OverlapMultBiclust.R** - additional code for simulation study 6.5
- **OverlapMultBiclust_X.rds** - dataset for simulation study 6.5
- **LungCancer.rds** - lung cancer gene expression dataset
- **trueColumnCluster.rds** - "true" column clustering for lung cancer gene expression dataset
- **AgYield.rds** - agricultural yield dataset
- **idr0.rds** - initial row clustering used for agricultural dataset analysis
- **idc0.rds** - initial column clustering used for agricultural dataset analysis
- **AgYield_MAR.R** - MAR version of Bayesian biclustering with agricultural yield dataset
- **Supplementary_Document_for_Paper_on_Bayesian_Biclustering.pdf** - tables and plots supporting the paper
