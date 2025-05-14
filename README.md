# Pattern_Continuation

Code for the manuscript [Data-driven continuation of patterns and their bifurcations](http://www.bjornsandstede.com/papers/Data_Driven_Continuation.pdf)

Authors: Wenjun Zhao, Samuel Maffa, Bjorn Sandstede

## Description

This software computes pattern statistics in spatially extended systems and uses them for bifurcation tracing. Depending on the task, typically it requires (1) an initial value problem solver which generates desired patterns, e.g. snapshots of Brusselator in spatial dimension 2 at a certain terminal time; (2) a pattern statistic to be used when performing the bifurcation tracing task, usually computed via $\alpha$-shapes of level sets.

For simulating spiral waves, this repository provides the [EZ-Spiral](http://homepages.warwick.ac.uk/~masax/Software/ez_software.html) code written by Dwight Barkley and [Matlab code written by Stephanie Dodson and Bjorn Sandstede](https://github.com/sandstede-lab/Spiral-Waves-Boundary-Sinks-and-Spectra) that were used in [Sandstede and Scheel. Spiral waves: linear and nonlinear theory](http://bjornsandstede.com/publications.html).

Node: the link for EZ-Spiral may be deprecated -- a copy of it is hosted [here](https://github.com/sandstede-lab/Pattern_Continuation/blob/main/Models/ezspiral-3.2.zip).

## Computation of the Wasserstein distance

This work uses the following implementation:

Wasserstein distance between 1D histograms: Niklas Kolbe (2024) https://github.com/nklb/wasserstein-distance

Transportation cost given a cost matrix: Antoine Rolet (2014) https://github.com/gpeyre/2017-ot-beginners/blob/master/matlab/mexEMD/mexEMD.cpp (based on an implementation by Nicolas Bonneeel).


## Requirement

MATLAB

## Source

```
git clone https://github.com/sandstede-lab/Pattern_Continuation.git
```

## Structure of scripts

The scripts are organized (tentatively) as the following:

* ```model.m```: takes a struct ```modelpar``` as input, with all model settings as its attributes, and outputs a pattern, which is either (1) positive/negative sets for spot/stripe related applications, or (2) tip point trajectory for spiral wave related applications.

* ```feature_evaluation.m```: takes two structs ```modelpar``` for evaluation of patterns, and ```featpar``` for mapping patterns to feature functions. The output is pattern statistics associated with the models.

* ```objective_evaluation.m```: takes two sets of pattern statistics as inputs, and output the distance metric between them, depending on the data type of inputs. Note that this function has dependency on ```ws_distance.m``` and ```medEMD/``` for evaluating Wasserstein distance.

* ```continuation.m```: this is the main function for continuation. It takes ```modelpar```, ```featpar```, and two function handles for feature evaluation and objective evaluation as inputs. In addition, it requries ```contpar```, which is a struct specifying details on running the continuation, and a starting point for the algorithm ```start```.

## Example usage in MATLAB

### To perform bifurcation tracing

The following file contains examples for (1) snaking in 1D Swift-Hohenberg model, and (2) stripes/spots in various 2D reaction-diffusion systems (Brusselator, Swift-Hohenberg, Gray-Scott, Schnakenberg).

Example usage: tracing out stripe/spot interface for Brusselator:

```
cd ReproduceCurves
Main_Brusselator
```
The file has starting points/directions already. An alternative is to use the automated starting point search tool:
```
cd ReproduceCurves
Test_IC_Brusselator
```

Spiral wave (Barkley):
```
cd ReproduceCurves
Main_Barkley
```
The file has starting points/directions already. An alternative is to use the automated starting point search tool:
```
cd ReproduceCurves
Test_IC_Barkley
```


Note that the spiral wave simulations may require compilation locally. Users may need to create their own Mex files by typing `make' in command line within the corresponding folder, such as DataGenerator/Spiral_Wave/Barkley.

### To use freezing method:

```
cd FreezingMethod
freezing_method
```

### To reproduce specific figures in paper
Scripts to reproduce each figure in paper are hosted under ReproduceFigures and labeled according to the index in paper. For example, to generate Figure 6 (simulations displayed together with our results):

```
Cd ReproduceFigures
Figure6
```






