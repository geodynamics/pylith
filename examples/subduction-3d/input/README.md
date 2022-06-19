# Directory for simulation input data that must be downloaded

This directory is simulation input data that is too large to include in the PyLith GitHub repository.
Follow the instructions below the download the necessary files.

## Cubit Mesh

```
curl -L -O https://github.com/geodynamics/pylith/releases/download/v3.0.1/examples-subduction-3d-mesh_tet.exo.gz
mv examples-subduction-3d-mesh_tet.exo.gz mesh_tet.exo.gz 
gunzip mesh_tet.exo.gz 
```
