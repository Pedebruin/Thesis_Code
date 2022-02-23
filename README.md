# Thesis_Code
The code necessary for graduation

This is a 1D Euler bernoulli simulation with piëzo patches

The beam model is generated based on Euler Bernoulli beam theory. The
piezoelec coupling is implemented as in the paper from 

Aktas, K. G.
Esen, I.
Doi:10.48084/etasr.3949

The author of this code is Pim de Bruin. 

The main model parameters are found at the top of the script.

The patches can be placed through modelSettings.patches. This is a vector
with the bottom location of each patch in natural coordinates. the length
of this vector determines the amount of patches. The size of the patches is 
determined through patchL. Same goes for the accelerometers and
modelSettings.Acc. 

The rest of the parameters are relatively self explanatory. If you have any
questions or find an error (very very probable) don't hesitate to send me a
message! @ P.E.deBruin@student.tudelft.nl or 0615331950

XXX Pim

This script requires:
    Symbolic math toolbox 
    Control system toolbox (For feedback command)
    Robust Control toolbox (For Hinf command)
    Statistics and machine learning toolbox (mvnrnd command)
    signal processing toolbox 
    System identification toolbox (goodnessOfFit command)
    Parallel computing toolbox (for parfor)

# How to run:
