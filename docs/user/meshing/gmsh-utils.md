(sec-user-meshing-gmsh-utils)=
# `gmsh_utils`

In Gmsh we use physical groups to associate cells with materials and mark entities for boundary conditions and faults.
The names of physical groups for materials must follow the syntax `material-id:TAG`, where TAG is the tag of the physical group.
PyLith includes a Python module `pylith.meshio.gmsh_utils` to make it easy to generate a PyLith compatible Gmsh file.

The function `create_material()` generates physical groups following the required naming convention of `material-id:TAG` given the tag and names of entities.
Similarly, the function `create_group()` will construct physical groups for boundary conditions and faults compatible with PyLith.

The physical groups for boundary conditions and faults must include entities at the topological dimension of the boundary condition as well as all lower dimensions.
For example, for a boundary condition on curves a physical group must include the entities on the curves as well as the vertices defining the curves.
For a boundary condition on surfaces a physical group must include the entities on the surfaces as well as the curves and vertices defining the surfaces.

## `GenerateMesh` Application Template

The `gmsh_utils` module also includes a application template object (Python abstract base class) called `GenerateMesh` for writing Python scripts that generate meshes using Gmsh.
The application template defines the steps for generating the mesh with a separate function (to be implemented by the user) for each step:

1. `initialize()`: Initialize Gmsh;
2. `create_geometry()`: Create the geometry (implemented in user application);
3. `mark()`: Create physical groups for boundary conditions, faults, and materials (implemented in user application);
4. `generate_mesh()`: Generate the finite-element mesh (implemented in user application);
5. `write()`: Save the mesh to a file; and
6. `finalize()`: Start the Gmsh graphical user interface, if requested, and then finalize Gmsh.

The command line arguments specify which step(s) to run, the output filename, and whether to invoke the Gmsh graphical user interface upon completing the steps:

:`--geometry`: Generate the geometry by calling `create_geometry()`.
:`--mark`: Create physical groups by calling `mark()`.
:`--generate`: Generate the mesh by calling `generate_mesh()`.
:`--write`: Save the mesh by calling `write()`.
:`--name`: Name of the mesh in Gmsh (default="mesh").
:`--filename=FILENAME`: Name of output mesh file (default="mesh.msh").
:`--ascii`: Write mesh to ASCII file (default is binary).
:`--cell=[tri,quad,tet,hex]`: Generate mesh with specified cell type.
:`--gui`: Start the Gmsh graphical user interface after running steps.

The application template always calls the `initialize()` and `finalize()` methods.
Additionally, the application will run any prerequisite steps.
For example, specifying `--generate` will trigger creating the geometry and physical groups before generating the mesh.

The application is discussed in more detail in the examples.

## MaterialGroup

`MaterialGroup` is a Python data class that holds information about a physical group associated with a material.
The data members include:

:tag (int): Integer tag for the physical group.
:entities (list): List (array) of entities for the material.

The `MaterialGroup` data class include a method `create_physical_group()` that will create a physical group from the information in the `MaterialGroup`.

## VertexGroup

`VertexGroup` is a Python data class that holds information about a physical group associated with a boundary or fault.
The data members include:

:name (str): Name for the physical group.
:tag (int): Integer tag for the physical group.
:dim (int): Dimension of the entities (0=points, 1=curves, 2=surfaces)
:entities: List (array) of entities for the boundary condition or fault.

The `VertexGroup` data class include a method `create_physical_group()` that will create a physical group from the information in the `VertexGroup`.

