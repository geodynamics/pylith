(sec-user-run-pylith-define-simulation)=
# Defining the Simulation

The parameters for PyLith are specified as a hierarchy or tree of components.
The application assembles the hierarchy of components from user input and then calls the `main` function in the top-level module in the same manner as a C or C++ program.
The behavior of the application is determined by the components included in the hierarchy as specified by the user.
The Pyre framework provides the interface for defining this hierarchy.
{ref}`sec-user-components` includes detailed descriptions of the components provided by PyLith.
Pyre properties correspond to simple settings in the form of strings, integers, and real numbers. 
Pyre facilities correspond to software modules. Facilities may have their own facilities (branches in the tree) and any number of properties.
See {numref}`fig:pyre:architecture` for the general concept of Pyre facilities and properties.

## Organization of Simulation Components

The components in a PyLith simulation generally fall into four main categories:

Topology
: Components associated with the spatial discretization of the domain, such as the finite-element mesh;

Physics
: Components specifying the physics to be solved, such as materials associated with a governing equation, bulk rheologies, boundary conditions, and fault interface conditions;

Physics Implementation
: Components that perform the finite-element operations, such as integration of the residual and system Jacobian; and

Observers
: Components that get notified of updates to the solution and state variables, such as writers for saving the solution to a file.

The physics components provide the point-wise functions (kernels) used by the physics implementation components, the auxiliary field, and the layout of the derived field (subfields computed from the auxiliary field and solution, such as stress and strain).

## Simulation Input and Output

{numref}`fig:pylith:simulation` shows the inputs and outputs for a PyLith simulation.
The user supplies:

1. Mesh information. This includes the topology of the finite-element mesh (coordinates of vertices and how the vertices are connected into cells), a material identifier for each cell, and sets of vertices associated with boundary conditions, faults, and output (for subsets of the mesh). This information can be provided using the PyLith mesh ASCII format (see {ref}`sec-user-file-formats-meshio-ascii` for the format specification) or by importing the information from Cubit, Gmsh, or LaGriT mesh generation software (see {ref}`sec-user-mesh` for more information).
2. A set of parameters describing the problem. These parameters describe the type of problem to be run, solver information, time-stepping information, boundary conditions, materials, etc. This information can be provided from the command-line or by using a `cfg` file.
3. Spatial databases specifying the values for the material properties and boundary conditions. Arbitrarily complex spatial variations in boundary and fault conditions and material properties may be given in the spatial database (see {ref}`sec-examples` and [Spatialdata Documentation](https://spatialdata.readthedocs.io).

PyLith writes solution information, such as solution fields and state variables, to either VTK files or HDF5/Xdmf files using the observer components.
ParaView and Visit as well as several other visualization tools can read both types of files.
Post-processing of output is generally performed using HDF5 files accessed via a Python script and the h5py package or a Matlab script.

:::{figure-md} fig:pylith:simulation
<img src="figs/pylith_simulation.*" alt="Diagram of a PyLith simulation" width="100%"/> 

PyLith requires a finite-element mesh (three different mechanisms for generating a mesh are currently supported), simulation parameters, and spatial databases (defining the spatial variation of various parameters).
PyLith writes the solution output to either VTK or HDF5/Xdmf files, which can be visualized with ParaView or Visit. Post-processing is generally done using the HDF5 files with Python or Matlab scripts.
:::

## Finite-Element Implementation User Interface

In specifying simulation parameters, some details of the finite-element implementation using the PETSc `DMPlex` is exposed to the user.
In this section we describe the data structures to give the user greater context for understanding what the parameters mean.

:::{tip}
See {ref}`sec-developer-code-layout` for a detailed discussion of the we organize the PyLith code.
:::

### Fields and Subfields

Finite-element coefficients for the finite-element basis functions (sometimes thought of as the values at vertices, on edges and faces, or in cells) are stored in a `Field`.
A `Field` is composed of a `Section`, which associates the points (vertices, edges, faces, and cells) with the finite-element coefficients, and a `Vec`, which is a vector storing the finite-element coefficients.
A `Field` may hold a single subfield, such as displacement, or it may hold several subfields, such as the density, shear modulus, and bulk modulus for an isotropic, linear elastic material.

Spatial discretization is specified for each subfield.
That is, each subfield within a `Field` can have a different discretization.
For example, a displacement field may use a second order discretization while a pressure field may use a first order discretization.
If we have uniform material properties, we use a zero order discretization (uniform values within a cell) to reduce the storage requirements.

The two main types of fields are the solution field and auxiliary fields.

#### Solution Field

The solution field contains all of the finite-element coefficients corresponding to the problem solution.
As discussed in the multiphysics finite-element formulation in {ref}`sec-user-petsc-fe-formulation`, if the governing equations have multiple unknowns, such as displacement and fluid pressure for poroelasticity, then the solution field will have multiple subfields.
See [Solution component](../components/problems/Solution.md) for details of the user interface and predefined containers for common subfield collections.

#### Auxiliary Field

We specify parameters for materials, boundary conditions, and fault interfaces using fields we refer to as the "auxiliary" field.
Each parameter (scalar, vector, tensor, or other) is held in a separate subfield.
We also store state variables in the auxiliary field, with each state variable as a different subfield.
This provides a single container for the collection of spatially varying parameters while maintaining the flexibility to specify the discretization of each parameter separately.

#### Discretization

The discretization of a field is given in terms of the topology (vertices, edges, faces, and cells) associated with the field and the basis order and quadrature order.
The basis order refers to the highest order in the basis functions.
For example, a basis order of 0 has just a constant and a basis order of 2 for a polynomial basis has constant, linear, and quadratic terms.

:::{warning}
Currently, the quadrature order **MUST** be the same for all subfields in a simulation.
This restriction may be relaxed in the future.
PyLith verifies that the quadrature order is the same for all subfields, and it will indicate if a subfield has a quadrature order that does not match the quadrature order of the first solution subfield.
:::

(sec-user-run-pylith-setting-parameters)=
## Setting PyLith Parameters

There are several methods for setting input parameters for the `pylith` executable: via the command line or by using a text file in `cfg` or `pml` format.
Both facilities and properties have default values provided, so you only need to set values when you want to deviate from the default behavior.

### Units

All dimensional parameters require units.
The units are specified using Python syntax, so square meters is `m**2`.
Whitespace is not allowed in the string, for units and dimensioned quantities are multiplied by the units string; for example, two meters per second is `2.0*m/s`.
Available units are shown in {numref}`tab-pyre-units`.

```{table} Pyre supported units. Aliases are in parentheses.
:name: tab-pyre-units
| Scale | Available Units |
| :------| :--------------- |
| length | meter (m), micrometer (um, micron), millimeter (mm), centimeter (cm), kilometer (km), inch, foot, yard, mile |
| time | second (s), nanosecond (ns), microsecond (us), millisecond (ms), minute, hour, day, year |
| mass | kilogram (kg), gram (g), centigram (cg), milligram (mg), ounce, pound, ton |
| pressure | pascal (Pa), kPa, MPa, GPa, bar, millibar, atmosphere (atm) |
```

### Using the Command Line

The `--help` command line argument displays links to useful resources for learning PyLith.

Pyre uses the following syntax to change properties from the command line.
To change the value of a property of a component, use `--COMPONENT.PROPERTY=VALUE`.
Each component is attached to a facility, so the option above can also be written as `--FACILITY.PROPERTY=VALUE`.
Each facility has a default component attached to it.
A different component can be attached to a facility by `--FACILITY=NEW_COMPONENT`.

PyLith's command-line arguments can control Pyre and PyLith properties and facilities, MPI settings, and PETSc settings.
You can get a list of all of these top-level properties along with a description using the `--help-properties` command-line argument.
To get information on user-configurable facilities and components, use the `--help-components` command-line argument.
To find out about the properties associated with a given component, use the `--COMPONENT.help-properties`} command-line argument.

```{code-block} console
$ pylith --problem.help-properties

# Show problem components.
$ pylith --problem.help-components

# Show bc components (bc is a component of problem).
$ pylith --problem.bc.help-components

# Show bc properties.
$ pylith --problem.bc.help-properties
```

### Using a `.cfg` File

Entering more than a few parameters via the command line is cumbersome.
You will generally find it easier to collect parameters into a `cfg` file.
The file is composed of one or more sections which are formatted as follows:

```{code-block} cfg
[pylithapp.COMPONENT1.COMPONENT2]
# This is a comment.

FACILITY3 = COMPONENT3
PROPERTY1 = VALUE1
PROPERTY2 = VALUE2
```

:::{tip}
We strongly recommend that you use `cfg` files for your work.
The files are syntax-highlighted in most editors, such as vim, Emacs, Atom, and VS Code.
:::

### Using a `.pml` File

A `pml` file is an XML file that specifies parameter values in a highly structured format.
It is composed of nested sections which are formatted as follows:

```{code-block} xml
<component name="COMPONENT1">
    <component name="COMPONENT2">
        <property name="PROPERTY1">VALUE1</property>
        <property name="PROPERTY2">VALUE2</property>
    </component>
</component>
```

XML files are intended to be read and written by machines, not edited manually by humans.
The `pml` file format is intended for applications in which PyLith input files are generated by another program, e.g., a GUI, web application, or a high-level structured editor.
This file format will not be discussed further here, but if you are interested in using `pml` files, note that `pml` files and `cfg` files can be used interchangeably; in the following discussion, a file with a `pml` extension can be substituted anywhere a `cfg` file can be used.

### Specification and Placement of Configuration Files

Configuration files may be specified on the command line:

```{code-block} console
$ pylith example.cfg
```

In addition, the Pyre framework searches for configuration files named `pylithapp.cfg` in several predefined locations.
You may put settings in any or all of these locations, depending on the scope you want the settings to have:

1. `$PREFIX/etc/pylithapp.cfg` for system-wide settings;
2. `$HOME/.pyre/pylithapp/pylithapp.cfg` for user settings and preferences; and
3. `./pylithapp.cfg`, for local settings.

All of the example problems are set up using configuration files and specific problems are defined by including the appropriate configuration file on the command-line.

:::{important}
The Pyre framework will search these directories for `cfg` files matching the names of components (for example, `timedependent.cfg`, `faultcohesivekin.cfg`, `greensfns.cfg`, etc) and will attempt to assign all parameters in those files to the respective component.
:::

:::{important}
Parameters given directly on the command line will override any input contained in a configuration file.
Configuration files given on the command line override all others.
The `pylithapp.cfg` files placed in (3) will override those in (2), (2) overrides (1), and (1) overrides only the built-in defaults.
:::

:::{tip}
See {ref}`sec-user-run-pylith-utilities` and {ref}`sec-user-run-pylith-parameter-gui` for several helpful utilities for viewing PyLith parameters and finding examples using specific features.
:::

% End of file