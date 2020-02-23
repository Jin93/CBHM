# CBHM
This repository provides the R code for Bayesian hierarchical model with a correlated prior (CBHM) for the analysis of early-phase oncology basket trials. The details of the method is described in Jin et al. (arXiv, 2020) https://arxiv.org/abs/2002.03007


1. CBHM_Prior_Specification.R:   specifying the prior distributions for CBHM in a general trial setting
2. Calibration.R:                calibrating tuning parameters for CBHM, EXNEX, Liu's two-stage method, BHM and independent analysis
3. Simulation.R:                 an example code for simulation.
4. CBHM_h.R:                     an example code for CBHM using H distance.
5. CBHM_kl.R:                    an example code for CBHM using KL distance.


## Citation
please cite:
Jin, J., Riviere, M. K., Luo, X., & Dong, Y. (2020). Bayesian Methods for the Analysis of Early-Phase Oncology Basket Trials with Information Borrowing across Cancer Types. arXiv preprint arXiv:2002.03007.
