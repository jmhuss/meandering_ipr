# meandering_ipr
#### Quantify and categorize meandering intensity of wind direction

As a submeso-sclae phenomenon, meandering winds are a potential turbulence-generating mechanism especially in the weak-wind stable boundary layer.
Quantifying the intensity of meandering for a series of given time windows may allow to draw meaningful conclusions about the state of the boundary layer.

The here provided method basically computes the interpercentile range (IPR) of the directional changes in any given window.

## Contents

This repository holds two R functions in [meand_ipr.R](/meand_ipr.R):
#### meand.ipr()
...quantifies meandering as a continuous measure
#### ipr.categ()
...is applied to the output of `meand.ipr()` and categorizes this measure, allowing to e.g. integrate the metric into boundary-layer regimes

Each function provides a basic description as well as a detailed explanation of input variables and output at the top of the respective code.
