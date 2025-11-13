# Directory for simulation input data that must be downloaded

This directory is simulation input data that is too large to include in the PyLith GitHub repository.
Follow the instructions below the download the necessary files.

## Cubit Mesh

```
curl -L -O https://github.com/geodynamics/pylith/releases/download/v5.0.0/examples-subduction-3d-mesh_tet.gmsh.gz
mv examples-subduction-3d-mesh_tet.gmsh.gz mesh_tet.gmsh.gz 
gunzip mesh_tet.gmsh.gz 

curl -L -O https://github.com/geodynamics/pylith/releases/download/v5.0.0/examples-subduction-3d-etopo2020_bedrock_local.tiff
mv examples-subduction-3d-etopo2020_bedrock_local.tiff etopo2020_bedrock_local.tiff 
```
