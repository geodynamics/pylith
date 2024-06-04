# Preface

This documentation is aimed at two categories of users: (1) scientists who prefer to use prepackaged and specialized analysis tools, and (2) experienced computational Earth scientists.
If you want to modify the source code, you will likely need to be familiar with object-oriented programming, Python and C++, and finite-element analysis.

## Citation

The Computational Infrastructure for Geodynamics (CIG) ([geodynamics.org](https://geodynamics.org/)) is making this software and source code available to you at no cost in hopes that the software will enhance your research in geophysics.
A number of individuals have contributed a significant portion of their careers toward the development of this software.
It is essential that you recognize these individuals in the normal scientific practice by citing the appropriate software release and peer-reviewed papers as well as making appropriate acknowledgments in talks and publications.
The preferred way to generate the list of publications (in BibTeX format) to cite is to run your simulations with the `--include-citations` command line argument, or equivalently, the `--petsc.citations` command line argument.
The `--version` command line argument will generate the BibTeX entries for the references mentioned below.

The following peer-reviewed paper discusses the development of PyLith:

- Aagaard, B. T., M. G. Knepley, and C. A. Williams (2013). A domain decomposition approach to implementing fault slip in finite-element models of quasistatic and dynamic crustal deformation, *Journal of Geophysical Research: Solid Earth*, 118, doi: 10.1002/jgrb.50217.

To cite the software and manual, use:

- Aagaard, B., M. Knepley, C. Williams (2023a), *PyLith v4.1.0.* Davis, CA: Computational Infrastructure of Geodynamics. DOI: 10.5281/zenodo.11479383.

- Aagaard, B., M. Knepley, C. Williams (2023b), *PyLith Manual, Version 4.1.0.* Davis, CA: Computational Infrastructure of Geodynamics. https://pylith.readthedocs.io/en/v4.1.0

## Publishing Models

Open research statements are now a common requirement when publishing research.
These support reuse, validation, and citation and often take the form of *Data availability*, *Data access*, *Code availability*, *Open Research*, and *Software* availability statements.
We recommend depositing input files that allow your published research to be reproduced and output model data in support of your research outcomes and figures.
In addition, consider depositing model files that may be reused by others.

Remember to cite software and data in your text as well as in your *Data availability* or similar statement.

Files should be deposited in an approved repository.

Additional information on [*Publishing*](https://geodynamics.org/software/software-bp/software-publishing) and repositories are available on the CIG website.

### Files

Common input files include `.cfg`, , `.spatialdb`, and mesh files `.exo`, `.msh`, or `mesh`.
Common output files include the solution fields and state variables as `.h5`, `.xmf`, `vtk`, and `vtu`..

### Example Statement

The configuration files, parameters of the simulation, and solution field for the models in this study are available at DOI (Authors X, Y, Z) under CC BY-NC-SA 4.0.

PyLith version 4.1.0 (Aagaard et al., 2013; Aagaard et al., 2023a; Aagaard et al., 2023b) used in these models is freely available under the MIT license for download through its software landing page https://geodynamics.org/resources/pylith or Zenodo (10.5281/zenodo.11479383).
The project is being actively developed on GitHub and can be accessed via https://github.com/geodynamics/pylith.

## Support

Current PyLith development is supported by the CIG, and internal GNS Science <https://www.gns.cri.nz/> and U.S. Geological Survey <https://www.usgs.gov/> funding.
Pyre development was funded by the Department of Energy's <https://www.energy.gov/energygov> Advanced Simulation and Computing program and the National Science Foundation's Information Technology Research (ITR) program.

This material is based upon work supported by the National Science Foundation under Grants No. 0313238, 0745391, and EAR-1550901.
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

## Acknowledgments

Many members of the community contribute to PyLith through reporting bugs, suggesting new features and improvements, running benchmarks, and asking questions about the software.
See the contributors list for each release for specific contributions.

## Request for Comments

Your suggestions and corrections can only improve this documentation.
Please report any errors, inaccuracies, or typos to the PyLith section of the CIG Community Forum <https://community.geodynamics.org/c/pylith> or create a GitHub pull request.
