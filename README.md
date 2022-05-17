# Cooperative Learning for Multi-view Analysis

This repository hosts the code for <b>cooperative learning</b>, a supervised learning method for multiple sets of features ("views"), as described in ["Cooperative Learning for Multi-view Analysis"](https://arxiv.org/abs/2112.12337).

## Overview
* Cooperative learning combines the usual squared error loss of predictions with an <b>"agreement" penalty </b> to encourage the predictions from different data views to agree. 
* By varying the weight of the agreement penalty, we get <b>a continuum of solutions that include the well-known early and late fusion approaches</b>. Cooperative learning chooses the degree of agreement (or fusion) in an adaptive manner, using a validation set or cross-validation to estimate test set prediction error.
* One version of our fitting procedure is <b>modular</b>, where one can choose different fitting mechanisms (e.g. lasso, random forests, boosting, neural networks) appropriate for different data views. 
* In the setting of cooperative regularized linear regression, the method combines the lasso penalty with the agreement penalty, yielding <b>feature sparsity</b>.
* The method can be especially powerful when the different data views <b>share some underlying relationship</b> in their signals that can be exploited to strengthen signal, while each view has its idiosyncratic noise that needs to be reduced.
    
## Usages
* The "cooperative_learning/" folder contains the implementation of cooperative learning. The subfolder contains the code to reproduce the results of our simulation and real data example studies. 
  *  Specifically, the subfolder "regularized_cooperative_regression/" contains code for the simulation experiments based on the factor model that generates Figure 2A and Figure 3 in the paper and Figure S1-6 in SI appendix; the subfolder "solution_sparsity/" contains code for the simulation experiments for exploring solution sparsity and generates Figure 2B in the paper; the subfolder "general_cooperative_learning/" contains code for the simulation studies with more distinct data modalities (e.g. imaging and omics data), and corresponds to Figure 4 and Figure 5 in the paper; the subfolder "more_than_two_data_views/" contains code for the simulation experiments for more than 2 data views, which corresponds to Figure 7-8 in the SI appendix; at last, the subfolder "real_data_examples/" contains code for the real data experiments of labor onset prediction and breast ductal carcinoma in situ and invasive breast cancer classification, which corresponds to Table 1 and 2 in the paper.

## Reference 
Daisy Yi Ding, Shuangning Li, Balasubramanian Narasimhan, Robert Tibshirani. <b>Cooperative Learning for Multi-view Analysis</b>. 2022. [[arXiv](https://arxiv.org/abs/2112.12337)]
