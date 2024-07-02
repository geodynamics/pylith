(sec-user-meshing-cubit)=
# Cubit

For many years we constructed Cubit meshes using journal files while leverating the Cubit APREPRO scripting language.
The Cubit Python interface has improved over this time, and with the use of Python scripting with Gmsh, we have replaced some of the Cubit journal files with Python scripts.
In a Python script, you can run any Cubit command using the `cubit.cmd()` function; however, for a more genuine Python experience, we recommend using the available Python functions.
Using Python scripts provides access to the complete programming language, making it easier to read in information from external files (for example, fault traces).

## Complex fault surfaces

Cubit has limited support for embedding fault surfaces into a volume domain.
We have found that Cubit will sometimes generate a suitable mesh with a single embedded fault, but usually it fails to generate a mesh when we have two embedded faults.
The workaround is to extend each fault surface to the edges of a volume (refer to `examples/crustal-strikeslip-3d` for an example). 

## Python scripting

In the PyLith examples, when we generate a Cubit mesh we create a Python script `generate_cubit.py`.
The script is setup to run from within Cubit or externally.
We usually write the script in an external editor with Python language support and then copy and paste commands into the Cubit Python command panel.

### Useful Cubit functions

The Cubit Python API is available in an appendix of the manual.
For example, the Cubit Python API for version 16.16 is available in [CubitInterface Namespace](https://cubit.sandia.gov/files/cubit/16.16/help_manual/WebHelp/cubithelp.htm).
s
`cubit.cmd(CMD: str)`
: Run Cubit command `CMD`.

`cubit.get_last_id(entity_type: str)`
: Get id of last entity create of type `entity_type`.

`cubit.get_id_from_name(name: str)`
: Get id of enty with name `name`.

`cubit.get_id_string(entity_ids: list)`
: Get string with list of ids for list of entity ids `entity_ids`.

`cubit.get_idless_signature(entity_type: str, id: int)`
: Get the idless signatude of a Cubit entity type `entity_type` with id `id`.

`cubit.vertex(id: int)`
: Get vertex given entity id `id`.

`cubit.curve(id: int)`
: Get curve given entity id `id`.

`cubit.surface(id: int)`
: Get surface given entity id `id`.

`cubit.volume(id: int)`
: Get volume given entity id `id`.

`cubit.body(id: int)`
: Get body given entity id `id`.
