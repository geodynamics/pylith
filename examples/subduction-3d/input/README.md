# Directory for simulation input data that must be downloaded

This directory is simulation input data that is too large to include in the PyLith GitHub repository.
Follow the instructions below the download the necessary files.

## Mesh files

```bash
# Gmsh
curl -L -O https://github.com/geodynamics/pylith/releases/download/v5.0.0/examples-subduction-3d-mesh_tet.msh.gz
mv examples-subduction-3d-mesh_tet.msh.gz mesh_tet.msh.gz 
gunzip mesh_tet.msh.gz 

# Cubit
curl -L -O https://github.com/geodynamics/pylith/releases/download/v5.0.0/examples-subduction-3d-mesh_tet.exo.gz
mv examples-subduction-3d-mesh_tet.exo.gz mesh_tet.exo.gz 
gunzip mesh_tet.exo.gz 

# Topography file used in generate_localdem.py
curl -L -O https://github.com/geodynamics/pylith/releases/download/v5.0.0/examples-subduction-3d-etopo2020_bedrock_local.tiff
mv examples-subduction-3d-etopo2020_bedrock_local.tiff etopo2020_bedrock_local.tiff 
```
