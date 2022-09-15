# meandering
Quantify and categorize meandering intensity of wind direction

As a submeso-sclae phenomenon, meandering winds are a potential turbulence-generating mechanism especially in the weak-wind stable boundary layer.
Quantifying the intensity of meandering for a series of given time windows may allow to draw meaningful conclusions about the state of the boundary layer.

The here provided method basically computes the interpercentile range (IPR) of the directional changes in any given window.

## Contents

This repository holds two R functions, the first of which (`meand.ipr()`) quantifies meandering as a continuous measure.
The second function (`ipr.categ()`), applied to the output of the first, categorizes this measure, allowing for e.g. application to generate boundary-layer regimes.
