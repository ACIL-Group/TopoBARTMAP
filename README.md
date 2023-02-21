# TopoBARTMAP


# About

# Quick Start

## Installation
Clone the code to needed folder and add to the path in MATLAB,
```
> addpath("Classes","Evaluation","Functions");
```
## Requirements
This work assumes that [CVAP: Cluster Validity Analysis Platform](https://www.mathworks.com/matlabcentral/fileexchange/14620-cvap-cluster-validity-analysis-platform-cluster-analysis-and-validation-tool) is assumed to be installed.

## Usage
*Input Data*: The data is assumed to be a matrix of shape #observationsX#features, ie., features/genes are assumed to be the columns, whilst observations are assumed to be rows. This is representation of data is different from how most gene expression data sets are stored in Gene Expression Omnibus (GEO), where data is a matrix of shape #genesX#observations.

*Classes*: If ground truth is available for the observations, then the clustering performance can be evaluated using external cluter validity index Adujsted Rand Index. Please refer to CVAP for further details.

### With classes
An example run file, *run_TopoBARTMAP.m* is provided in **Functions** to use with gene expression data from GEO. This run file assumes data is of shape #genesX#observations. *run_TopoBARTMAP* returns an instance of TopoBARTMAP and the indices value.
```
> [y,TBM] = run_TopoBARTMAP(data,classes);
```

### Without classes
run 


# Citing

If you use TopoBARTMAP in your research, please cite the [TopoBARTMAP paper](https://doi.org/10.1016/j.neunet.2022.12.010).

*Yelugam, R., da Silva, L.E.B. and Wunsch II, D.C., 2023. Topological biclustering ARTMAP for identifying within bicluster relationships. Neural Networks, 160, pp.34-49.*

In BibTeX format:

```tex
@article{YELUGAM202334,
title = {Topological biclustering ARTMAP for identifying within bicluster relationships},
journal = {Neural Networks},
volume = {160},
pages = {34-49},
year = {2023},
issn = {0893-6080},
doi = {https://doi.org/10.1016/j.neunet.2022.12.010},
url = {https://www.sciencedirect.com/science/article/pii/S0893608022005020},
author = {Raghu Yelugam and Leonardo Enzo {Brito da Silva} and Donald C. {Wunsch II}},
keywords = {Biclustering, Topological data analysis, Adaptive resonance theory (ART), Gene expression, Gene Co-expression}
}
```
