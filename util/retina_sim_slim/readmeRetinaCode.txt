

Greetings, 


The folder 'Simulation' contains the few m-files I wrote to simulate a piece of retina.

genFilters: generates filters for DoG-type RGC models. It currently supports only DoG 
center-surround cells, but as I got it we won't need much more.
The placement of these is currently either random (within the central 80% along X- and Y-axis, 
such that not too much gets clipped off towards the retina patch borders) or on some lattice. 
        As I had started out more fancy than might be needed now, I handed over parameters for 
Wishart-priors over the DoG-component covariance matrices to allow some variability of RF form. 
Other inputs of interest control filter balancedness (integral over the filter values) by 
changing relative DoG-component hights. Feel free to cut as much of this out as desired. 
I enforce 'sparsity' of the filter by setting all filter values below a threshold to zero. Idea was
to save memory/computation for very large retina patches by using sparse(W). 

retSim: simulates the retina given by filters W and settings for a sigmoidal functions to map filter
responses to RGC outputs. I didn't fiddle too much yet with introducing noise correlations. The current
code as is would generate them from spatially correlated noise (another pink noise image added on a 
trial-by-trial-varying basis on top of the actual stimulus) that are yet not implemented
.

spatialPattern: a stack-exchange file for generating white, pink or brown noise images. Mostly chosen 
because it was available and seemingly ready for use! Slight changes from my side to allow returning 

entire image stacks at once.  

extractKeyStats: some file for plotting full/signal/noise correlations from spiking data. Out-commented
part plots correlations of means over trials. Signal correlations computed by shuffling trials individually
for each of the n neurons and then simply using Matlab's corr(x). 

demo_retina: Shows how one would use this code. Should run right away if the other files are in the path. 
 
The code is part of a larger repository also with functions for fitting MaxEnt models and sampling
from the fits, so there might be confusing parts left in the file version I sent you that most likely have 
something to do with these other project parts. Overall, especially the retina simulation part that you now
hold in your hands was mostly written without thinking that anybody beyond me and Jakob would ever get to see 
this thing. Not sure if all decisions/conventions make sense outside of my head. 


The folder 'Data' contains figures (generated at heart with extractKeyStats.m) for the 'real' and desired
full/signal/noise correlations and their histograms. SalamanderDataStatistics.mat contains the values that
go along with the plots. I tried to give it sensible names. The distribution of population activity counts
starts at K = 0 and goes up to K = 40, even though these high counts never occured in the data. 



Anything missing or unclear? -> marcel.nonnenmacher@tuebingen.mpg.de ! 


Regards, Marcel