Adaptive Reverse Constrained IVA Gauss and Threshold Free Constrained IVA Gauss, were implemented. About 15 hyperparameters per algorithm may be set. Works for both old and new GUI and batch mode.

Implementation of adaptive-reverse constrained IVA with multivariate Gaussian distribution (ar-cIVA-G) that incorporates prior informationabout the sources into the IVA cost function. By alternating between a conservative scheme and an assertive scheme, ar-cIVA-G optimally controls the effect of each reference on the corresponding estimated source. There is no need for users to specify the degree of similarity between the estimate and the reference signal. For a general description of the algorithm and its relationship with others, see http://mlsp.umbc.edu/jointBSS_introduction.html.

Implementation of threshold-free constrained IVA with multivariate Gaussian distribution (tf-cIVA-G) that utilizes the references as regularization for the IVA cost function. Both the similarity between the reference and the corresponding source and the (dis)similarity between that reference and the other sources are exploited. tf-cIVA-G eliminates the need for constraint-threshold selection.

Implementation of threshold-free constrained IVA and adaptive-reverse constrained IVA is in accordance with following publication:
Vu, Trung, Francisco Laport, Hanlu Yang, Vince D. Calhoun, and Tulay Adali. 
"Constrained independent vector analysis with reference for multi-subject fMRI analysis," 
IEEE Transactions on Biomedical Engineering, 12 pages, July, 2024, DOI: 10.1109/TBME.2024.3432273.
