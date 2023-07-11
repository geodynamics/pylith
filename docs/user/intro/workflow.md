# PyLith Workflow

PyLith is one component in the process of investigating problems in tectonics ({numref}`fig:workflow:summary`).
Given a geological problem of interest, you must first provide a geometrical representation of the desired domain.
Once you define the geometry, ou need to discretize it as a finite-element mesh.
PyLith presently provides three mesh importing options: CUBIT Exodus format, Gmsh ASCII and binary files, and PyLith mesh ASCII format.
The modeling of the physical processes of interest is performed by a code such as PyLith.
Present output consists of VTK or HDF5/Xdmf files which can be used by a number of visualization codes (e.g., ParaView, Visit, and Matlab).

:::{figure-md} fig:workflow:summary
<img src="figs/workflow.*" alt="Workflow involved in going from geologic structure to problem analysis." width="100%" />

Workflow involved in going from geologic structure to problem analysis.
:::

% End of file
