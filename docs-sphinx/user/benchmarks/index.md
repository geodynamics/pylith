# Benchmarks

:::{warning}
None of the benchmark input files in the PyLith benchmarks repository on GitHub have been updated for v3.0.
:::

## Overview

The Crustal Deformation Modeling and Earthquake Source Physics Focus Groups within the Southern California Earthquake Center and the Short-Term Tectonics Working Group within CIG have developed a suite of benchmarks to test the accuracy and performance of 3D numerical codes for quasistatic crustal deformation and earthquake rupture dynamics.
The benchmark definitions for the quasistatic crustal deformation benchmarks are posted on the CIG website at Short-Term Tectonics Benchmarks <https://geodynamics.org/cig/workinggroups/short/workarea/benchmarks/> and the definitions for the earthquake rupture benchmarks are posted on the SCEC website <https://strike.scec.org/cvws/cgi-bin/cvws.cgi>.
This suite of benchmarks permits evaluating the relative performance of different types of basis functions, quadrature schemes, and discretizations for geophysical applications.
The files needed to run the 3D benchmarks are in the CIG GitHub Repository <https://github.com/geodynamics/pylith_benchmarks>.
In addition to evaluating the efficiency and accuracy of numerical codes, the benchmarks also make good test problems, where users can perform simulations based on actual geophysical problems.
The benchmarks are performed at various resolutions and using different element types.
By comparing the runtime and accuracy for different resolutions and element types, users can evaluate which combination will be best for their problems of interest.

:::{toctree}
quasistatic-strikeslip.md
savage-prescott.md
scec-dynrup.md
:::
