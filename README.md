
# sMSNF-EM: Switching Multiscale Numerical Integration Filter - Expectation Maximization <br/>


Learn parameters for a switching multiscale dynamical system with Gaussian and Poisson observations using a switch EM framework with a multiscale numerical integration filter (MSNF) embedded.


# Publication

Kim D., Song C.Y., Hsieh H-L., Shanechi M.M., "Unsupervised learning of multiscale switching dynamical system models from multimodal neural data", Journal of Neural Engineering, Apr. 2026

Link to the paper: https://doi.org/10.1088/1741-2552/ae5688

# Usage guide
## Dependencies
Code was developed and tested on MATLAB R2021b and the following packages are used:
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

## Initialization
Add the source directory and its subdirectories to the path. Adding the experimental directory is optional unless you want to run the example script. You can run init.m to do this. 

## Main learning function
The main learning function [source/fitSwitchMultiscaleEM.m](source/fitSwitchMultiscaleEM.m). A complete usage guide is available in the function. Also see [source/scriptMSNF.m](source/example%20scripts/scriptMSNF.m). The following shows an example case: 
```
trainSet = struct('fs',fs,'obsNt',nt, 'obsYt', yt);
paramSet = struct('dimXtEst',nx,'dimStEst',ns);
mdlS = fitSwitchEM(trainSet,paramSet,nIter);
```

Inputs:
- nt is the input to trainSet.obsNt which is a dimension x time matrix with spiking activity observation data. It can also be a cell array of dimension x varying time length matrices for trial based data.
- yt is the input to trainSet.obsYt which is in the same format as yt but for field observation data.
- fs is the input to trainSet.fs which is the sampling frequency in Hz.
- nx is the input to prmS.dimXtEst which is the assumed latent state xt dimension.
- ns is the input to prmS.dimStEst which is the assumed number of regimes.
- nIter is the number of iterations to run EM.

Output:
- mdlS: a structure containing learned model parameters and any tracked metrics. The final parameters can be obtained by mdlS.thetaCell{end}. See code for details.

## Estimating latent and regime states in test data
Once a model is learned, you can apply the model to new data to estimate latent and regime states. See [source/filterStates.m](source/filterStates.m) for more details on outputs.
```
theta = mdlS.thetaCell{end};
mngr = filterStates(theta,ytTest,ntTest);
PStDec = mngr.PStEst;
```
Input:
- theta: model parameters
- ntTest: spike observations dimension x time.
- ytTest: field observations dimenxion x time.

Outputs:
- mngr: structure containing filtered states and estimation covariances
- mngr.PStDec: decoded probability of regimes. Choose regime with highest probability to get decoded regimes.
- mngr.xDec: decoded latent states xt.

# Licence
Copyright (c) 2026 University of Southern California  
See full notice in [LICENSE.md](LICENSE.md)
DongKyu Kim and Maryam M. Shanechi  
Shanechi Lab, University of Southern California