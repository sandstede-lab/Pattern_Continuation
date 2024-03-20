# Pattern_Continuation

Code for data-driven continuation of patterns and their bifurcations (link XXX)

Authors: Wenjun Zhao, Sam Maffa, Bjorn Sandstede

For questions/comments please contact XXX

## Description

This software computes pattern statistics in spatially extended systems and uses them for bifurcation tracing. Depending on the task, typically it requires (1) an initial value problem solver which generates desired patterns, e.g. snapshots of Brusselator in spatial dimension 2 at a certain terminal time; (2) a pattern statistic to be used when performing the bifurcation tracing task, usually computed via $\alpha$-shapes of level sets.

For simulating spiral waves, this repository provides the [EZ-Spiral](http://homepages.warwick.ac.uk/~masax/Software/ez_software.html) code written by Dwight Barkley and [Matlab code written by Stephanie Dodson and Bjorn Sandstede](https://github.com/sandstede-lab/Spiral-Waves-Boundary-Sinks-and-Spectra) that were used in [Sandstede and Scheel. Spiral waves: linear and nonlinear theory](http://bjornsandstede.com/publications.html).

The implementation for computing the Wasserstein distance can be found on: 
Niklas Kolbe (2024). Wasserstein distance (https://github.com/nklb/wasserstein-distance), GitHub. Retrieved March 20, 2024.

## Requirement

MATLAB

## Source

```
git clone https://github.com/WenjunZHAOwO/Pattern_Continuation.git
```

## Example usage in MATLAB

# To perform bifurcation tracing

The following file contains examples for (1) snaking in 1D Swift-Hohenberg model, and (2) stripes/spots in various 2D reaction-diffusion systems (Brusselator, Swift-Hohenberg, Gray-Scott, Schnakenberg).

```
cd BifurcationTracing
Main_adapt_bisection
```

The following file contains examples for spiral waves:
```
cd BifurcationTracing
Main_Barkley_bisection
```
Note that users may need to create their own Mex files by typing `make' in command line within the corresponding folder, such as DataGenerator/Spiral_Wave/Barkley.


# To reproduce figures in paper

The following file shows how to display results together with solutions of IVP problems in parameter space for Brusselator:

```
cd DisplayResults
overlay_simulations_Brusselator
```
Visualization for Barkley:
```
cd DisplayResults
overlay_simulations_Barkley
```





