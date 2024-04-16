# Sparse Non-Uniform Optical Phased Array (OPA) Design
 This repository is for my thesis work [*Sparse Non-Uniform Optical Phased Array Design*](https://repository.tudelft.nl/islandora/object/uuid%3A4a67e143-0004-48cb-963f-c6c874282a9d). The general sparse OPA design problem can be classified into several cases according to the architecture of the OPA (the geometry and the controllability of the excitations to each element). Depending the corresponding mathematical description for each case, the formulation of the design problem can be very different. 
 My work starts from the background and attempt to address one of the cases, namely the 1D linear amplitude control structure. This case is the simplest to handle from the optimization perspective. The problem is solved 

 The problem formulation is still improvable, especially regarding the constraints on the beampattern: 
 1. Hard (equality) constraint can force the peak of beampattern to be unit. 
 2. Beampattern can be bounded to control the peak sidelobe


## Explanation of this repository
The explanation of this repository is as following:
Only the folder "**src/After_graduation**" contains meaningful scripts, which are used to generate plots included in my thesis. I implemented three algorithms solving the problem formulation, LASSO, thresholding and reweighted l1 norm minimization. 
The implementations of algorithms are in *Scripts*. The involving functions are *Functions*. The results and the parameters that are used to generate these results are in *Results and Parameters*. 
### Scripts
Three scripts contain the code for minimization problem and the code to generate the plots. 
- LASSO.m
- thresholding.m
- reweighted_l1norm.m

Loading the stored results and skipping the section running optimization, the script can generate plots. 
Adjusting parameters and solving the optimization problem, the script can be adapted to different requriements. 

### Functions
- desiredbeam.m (To generate an ideal beampattern according to the user-defined requirements)
- threshold.m (To threshold an excitation vector to different level of sparsity)

### Results and Parameters
- LASSO_parameters.mat (The results and the paramters used for LASSO)
- reweighted_results.mat (The results and the parameters used for reweigthed l1 norm minimization)

### Animation
Generated by script "thresholding.m", demonstrate the evolution of the beampattern with thresholding the excitation vector. 
- thresholding_animation.avi

Good luck with exploring!