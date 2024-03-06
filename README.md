# Pattern_Continuation

Code for data-driven continuation of patterns and their bifurcations (link XXX)

Authors: Wenjun Zhao, Sam Maffa, Bjorn Sandstede

For questions/comments please contact XXX

## Description

This software computes pattern statistics in spatially extended systems and uses them for bifurcation tracing. Depending on the task, typically it requires (1) an initial value problem solver which generates desired patterns, e.g. snapshots of Brusselator in spatial dimension 2 at a certain terminal time; (2) a pattern statistic to be used when performing the bifurcation tracing task, usually computed via $\alpha$-shapes of level sets.

## Requirement

MATLAB

## Source

```
git clone https://github.com/WenjunZHAOwO/Pattern_Continuation.git
```

## Example usage in MATLAB

To perform bifurcation tracing, please see the example below for Brusselator:

```
cd BifurcationTracing
Main_adapt_bisection
```

To reproduce figures in paper (display results together with solutions of IVP problems in parameter space), please see the example below for Brusselator:

```
cd DisplayResults
overlay_simulations_Brusselator
```





