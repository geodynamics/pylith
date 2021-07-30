# Preface

## About This Document

This document is organized into two parts.
The first part begins with an introduction to PyLith and discusses the types of problems that PyLith can solve and how to run the software; the second part provides appendices and references.

## Who Will Use This Documentation

This documentation is aimed at two categories of users: scientists who prefer to use prepackaged and specialized analysis tools, and experienced computational Earth scientists.
Of the latter, there are likely to be two classes of users: those who just run models, and those who modify the source code.
Users who modify the source are likely to have familiarity with scripting, software installation, and programming, but are not necessarily professional programmers.

## Conventions

:::{warning}
This is a warning.
:::

:::{important}
This is something important.
:::

:::{tip}
This is a tip, helpful hint, or suggestion.
:::

For features recently added to PyLith, we show the version number when they
were added.
*New in v3.0.0dev*

### Command Line Arguments

Example of a command line argument: `--help`.

### Filenames and Directories

Example of filenames and directories: `pylith`, `/usr/local`.

### Unix Shell Commands

Commands entered into a Unix shell (i.e., terminal) are shown in a box.
Comments are delimited by the # character. We use `$` to indicate the bash shell prompt.

```{code-block} bash
#This is a comment.
$ ls -l
```

### Excerpts of cfg Files

Example of an excerpt from a `.cfg` file:

```{code-block} cfg
# This is a comment.
[pylithapp.problem]
timestep = 2.0*s
bc = [x_pos, x_neg]
```

## Citation

The Computational Infrastructure for Geodynamics (CIG) ([geodynamics.org](https://geodynamics.org/)) is making this source code available to you at no cost in hopes that the software will enhance your research in geophysics.
A number of individuals have contributed a significant portion of their careers toward the development of this software.
It is essential that you recognize these individuals in the normal scientific practice by citing the appropriate peer-reviewed papers and making appropriate acknowledgments in talks and publications.
The preferred way to generate the list of publications (in BibTeX format) to cite is to run your simulations with the `--include-citations` command line argument, or equivalently, the `--petsc.citations` command line argument.
The `--help-citations` command line argument will generate the BibTeX entries for the references mentioned below.

The following peer-reviewed paper discussed the development of PyLith:

-   Aagaard, B. T., M. G. Knepley, and C. A. Williams (2013). A domain
    decomposition approach to implementing fault slip in finite-element models
    of quasistatic and dynamic crustal deformation, *Journal of Geophysical
    Research: Solid Earth*, 118, doi: 10.1002/jgrb.50217.

To cite the software and manual, use:

-   Aagaard, B., M. Knepley, C. Williams (2017), *PyLith v3.0.0dev.* Davis,
    CA: Computational Infrastructure of Geodynamics. DOI:
    10.5281/zenodo.XXXXXX.

-   Aagaard, B., M. Knepley, C. Williams (2017), *PyLith User Manual, Version
    3.0.0dev.* Davis, CA: Computational Infrastructure of Geodynamics. URL:
    geodynamics.org/cig/software/github/pylith/v3.0.0dev/pylith-3.0.0dev_manual.pdf

## Support

Current PyLith development is supported by the CIG, and internal GNS Science <https://www.gns.cri.nz/> and U.S. Geological Survey <https://www.usgs.gov/> funding.
Pyre development was funded by the Department of Energy's <https://www.energy.gov/energygov> Advanced Simulation and Computing program and the National Science Foundation's Information Technology Research (ITR) program.

This material is based upon work supported by the National Science Foundation under Grants No. 0313238, 0745391, 1150901, and EAR-1550901.
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

## Acknowledgments

Many members of the community contribute to PyLith through reporting bugs, suggesting new features and improvements, running benchmarks, and asking questions about the software.
In particular, we thank Surendra Somala for contributing to the development of the fault friction implementation.

## Request for Comments

Your suggestions and corrections can only improve this documentation. Please
report any errors, inaccuracies, or typos to the CIG Short-Term Tectonics
email list <cig-short@geodynamics.org> or create a GitHub pull request.
